module analysis;

import std.algorithm;
import std.complex;
import std.json;
import std.math;
import std.range;
import std.stdio;
import std.typecons;

import deint;
import utils;

alias R = double;
alias C = Complex!R;


/** Struct for transceiver model
*/
struct Transceiver
{
    R g;                    /// g_{n}
    C delegate(C)   alpha,  /// alpha(x)
                    beta;   /// beta(x)
    R rho11;                /// rho
    R noisePower;           /// N_{thermal}
    uint P;                 /// Order of canceller
    R OS_RATE;              /// Ratio of bandwidth of all active subcarriers and sampling rate
    uint M_QAM;             /// use M-ary QAM
}


// Parameters of analysis
struct Parameter
{
    Transceiver T1,     /// terminal#1
                T2;     /// terminal#2
    R rho21;            /// rho_{21}
}


/** Struct for analysis result
*/
struct Result
{
    R[2] SICR;                  /// performance of digital canceller
    R[2] SER;                   /// symbol error rate
    R[2] achRate;               /// achievable rate
    R[][2] nlpowers;            /// powers of nonlinearity
    R[2] txpower;               /// transmit power (nu_n^2)
    R[2] residual;              /// residual self-interference powers
    R[2] N_tot;                 /// total noise power
    R[2] N_NL;                  /// power of nonlinear noise caused by the LNA

    R[2] SER_HD;                /// SER on half duplex
    R[2] achRate_HD;            /// achievable rate on half duplex


    JSONValue toJSONValue() const
    {
        JSONValue res;
        res["SICR1"] = 10*log10(SICR[0]);
        res["SICR2"] = 10*log10(SICR[1]);
        res["SER1"] = SER[0];
        res["SER2"] = SER[1];
        res["Rate1"] = achRate[0];
        res["Rate2"] = achRate[1];
        res["NLPowers1"] = JSONValue(nlpowers[0][].map!(a => JSONValue(a)).array());
        res["NLPowers2"] = JSONValue(nlpowers[1][].map!(a => JSONValue(a)).array());
        res["SER1_HD"] = SER_HD[0];
        res["SER2_HD"] = SER_HD[1];
        res["Rate1_HD"] = achRate_HD[0];
        res["Rate2_HD"] = achRate_HD[1];
        res["TXPower1"] = txpower[0];
        res["TXPower2"] = txpower[1];
        res["ResPower1"] = residual[0];
        res["ResPower2"] = residual[1];
        res["N_tot1"] = N_tot[0];
        res["N_tot2"] = N_tot[1];
        res["N_NL1"] = N_NL[0];
        res["N_NL2"] = N_NL[1];
        return res;
    }
}


/** Analyze the performance of full-duplex terminals
*/
Result analyze(Parameter param)
{
    immutable R[2]
        rho11 = [param.T1.rho11, param.T2.rho11],
        g = [param.T1.g, param.T2.g],
        N_T = [param.T1.noisePower, param.T2.noisePower],
        OS_RATE = [param.T1.OS_RATE, param.T2.OS_RATE];

    immutable uint[2] P = [param.T1.P, param.T2.P],
                      M_QAM = [param.T1.M_QAM, param.T2.M_QAM];

    auto alpha = [param.T1.alpha, param.T2.alpha],
         beta = [param.T1.beta, param.T2.beta];

    immutable MAX_P = max(P[0], P[1]) + 10;


    // helper for make both T1 -> T2 and T2 -> T1 case
    X[2] make12_21(X)(X delegate(uint i, uint j) dg) {
        X[2] arr = [dg(0, 1), dg(1, 0)];
        return arr;
    }


    // Define numerical integration object
    auto deintRayleigh = makeDEInt!R(0, R.infinity, Yes.isExpDecay, 401, -4, 4)
                        .withWeight((R r) => 2 * r * exp(-r^^2));
    
    // Define expectation operator int_0^{infty} func(r) 2*r * e^{-r^2} dr = 1/pi * int_{C} func(x) e^{-|x|^2} dx 
    T expectOnRayleigh(T)(T delegate(C) func)
    {
        return deintRayleigh.integrate((R r) => func(C(r)));
    }


    // Step2: Compute the output power nu_1^2 from the n-th terminal,
    //        and nu1 = sqrt(nu_1^2)
    immutable R[2] nu = make12_21((i, j) => sqrt(expectOnRayleigh(x => alpha[i](g[i] * x).sqAbs )));

    // Step3: Compute the Fourier coefficients alpha_{n,p}
    immutable C[][2] alpha_np = make12_21((i, j) => iota(1, MAX_P+2, 2).map!(p => expectOnRayleigh(x => alpha[i](g[i] * x) * psi(p, x).conj)).array.idup);

    // Step4: Compute the linear gain and distortion power
    immutable C[2] beta_11 = make12_21((i, j) =>  1/(nu[i] * rho11[i]) * expectOnRayleigh(x => beta[i](nu[i] * rho11[i] * x) * x.conj));
    immutable R[2] N_NL = make12_21((i, j) => expectOnRayleigh(x => param.T1.beta(nu[i] * rho11[i] * x).sqAbs) - sqAbs(beta_11[i] * rho11[i] * nu[i]));

    // Step5: Compute SICR
    immutable R[2] I_11 = make12_21((i, j) => sqAbs(beta_11[i] * rho11[i] * nu[i]));
    immutable R[2] I_11_R = make12_21((i, j) => sqAbs(beta_11[i] * rho11[i]) * (nu[i]^^2 - alpha_np[i][0 .. P[i]/2+1].map!sqAbs.sum));
    immutable R[2] N_tot = make12_21((i, j) => beta_11[i].sqAbs * N_T[i] + N_NL[i]);

    Result result;
    result.SICR = (I_11[0] + N_tot[0]) / (I_11_R[0] + N_tot[0]);
    result.txpower = nu;
    result.residual = I_11_R;
    result.N_tot = N_tot;
    result.N_NL = N_NL;


    // Step6: Compute SER and achievable rate
    auto deintFreq = makeDEInt!R(-0.5, 0.5, No.isExpDecay, 51);
    // auto deintLambdaSINR = makeDEInt!R(0, R.infinity, No.isExpDecay, 201);
    auto deint0IExp = makeDEInt!R(0, R.infinity, Yes.isExpDecay, 101, -8, 8);

    R expectOnSIDNR(R delegate(R) func, size_t i, size_t j, Flag!"isHD" isHD = No.isHD)
    {
        return deintFreq.parallelIntegrate((R f){
            // In this implementation, we apply transformation to the variable $x$ for accurate integration of $\int_{0}^{\infty} f(x) p_{\lambda_{sinr}}(x) dx$.
            // The applied transformation is $t = N_{tot,1}(f) / ( x U_{21} )$   <=>   $x = N_{tot,1}(f) / ( t U_{21} )$.
            return deint0IExp.integrate((R t) {
                R LambdaSDR = 0;
                foreach(p; iota(3, MAX_P+2, 2))
                    LambdaSDR += sqAbs(alpha_np[j][p/2] / alpha_np[j][0]) * nTimesConvRectSpectrum(f, p);

                R I_11_R_f = 0;
                if(! isHD) {    // for full-duplex
                    foreach(p; iota(P[i] + 2, MAX_P+2, 2))
                        I_11_R_f += sqAbs(beta_11[i] * rho11[i] * alpha_np[i][p/2]) * nTimesConvRectSpectrum(f, p);
                }

                // LNA is not saturated on HD system
                immutable b11 = isHD ? C(1) : beta_11[i];
                immutable nnl = isHD ? 0 : N_NL[i];

                immutable R N_tot_f = sqAbs(b11) * N_T[i] / OS_RATE[i] + nnl * nTimesConvRectSpectrum(f, 3);
                immutable R U_21 = sqAbs(b11 * param.rho21 * alpha_np[j][0]);

                immutable NtI = N_tot_f + t * I_11_R_f;
                immutable pdf = N_tot_f * (NtI + I_11_R_f) / NtI^^2 * exp(-t);
                return func(U_21 * t / (N_tot_f + U_21 * t * LambdaSDR)) * pdf;
            });
        });
    }

    result.SER = make12_21((i, j) => expectOnSIDNR(x => qamSER(x, M_QAM[i]), i, j));
    result.achRate = make12_21((i, j) => expectOnSIDNR(x => log2(1 + x), i, j));
    result.nlpowers = make12_21((i, j) => alpha_np[i].map!sqAbs.array ~ (nu[i] - alpha_np[i].map!sqAbs.sum) );

    result.SER_HD = make12_21((i, j) => expectOnSIDNR(x => qamSER(x, M_QAM[i]), i, j, Yes.isHD));
    result.achRate_HD = make12_21((i, j) => expectOnSIDNR(x => log2(1 + x), i, j, Yes.isHD));

    return result;
}

