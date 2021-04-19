import analysis;
import utils;

import std.algorithm;
import std.complex;
import std.conv;
import std.format;
import std.json;
import std.math;
import std.range;
import std.stdio;


auto makeRapp(R vsat, R sf)
{
    return delegate(C x){
        return rapp(x, vsat, sf);
    };
}


auto makeSaleh(R vsat, R phisat)
{
    //return delegate(C x) {
    //    immutable r = std.complex.abs(x),
    //              u = x / r;    // unit vector

    //    return vsat * normalized_saleh(r / vsat, phisat) * u;
    //};
    immutable R a1 = 1,
                a2 = 2*phisat/4/vsat/vsat,
                b1 = 1.0/4/vsat/vsat,
                b2 = b1;

    return delegate(C x) {
        immutable r = std.complex.abs(x),
                  u = x / r;    // unit vector

        return cast(C)(a1 * x / (1 + b1 * r^^2) * std.complex.expi(a2*r^^2/(1 + b2*r^^2)));
    };
}


void main()
{
    immutable R OS_RATE = 8;                          // Oversampling rate of simulation

    immutable R vsat_alpha = 30.0.dBmToVolt;          // Output saturation level of rapp is 30 dBm
    immutable R vsat_beta = (-6.0).dBmToVolt;         // Input saturation level of LNA is -6 dBm

    static auto makeRapp3(R v) { return makeRapp(v, 3); }
    static auto makeSalehPI6(R v) { return makeSaleh(v, PI/6); }
    static auto makeSaleh0(R v) { return makeSaleh(v, 0); }

    typeof(&makeRapp3)[string] genNLFuncAA;
    genNLFuncAA["Rapp"] = &makeRapp3;
    genNLFuncAA["Saleh"] = &makeSalehPI6;
    genNLFuncAA["Saleh_noPM"] = &makeSaleh0;

    Parameter defaultSetting;
    with(defaultSetting) {
        T1.g = vsat_alpha / 7.0L.dBToVGain;             // 7-dB IBO
        T1.alpha = makeRapp(vsat_alpha, 3);              // AM-AM characteristic
        T1.beta = makeRapp(vsat_beta, 3);               // ditto
        T1.rho11 = (-50.0L).dBToVGain;
        T1.noisePower = (-174 + 10*log10(20*10^^6 * OS_RATE) + 4).dBmToVolt^^2; // -174 dBm/Hz * (20MHz * 8 Oversampling) * 4-dB Noise Figure of LNA
        T1.P = 7;
        T1.OS_RATE = OS_RATE * 64 / 52;
        T1.M_QAM = 16;
        T1.IRR = real.infinity;

        T2 = T1;        // symmetric (T1 == T2)
        rho21 = (-70.0).dBToVGain;
    }

    foreach(iden, makeNLFunc; genNLFuncAA)
    {
        Parameter setting = defaultSetting;
        setting.T1.alpha = makeNLFunc(vsat_alpha);
        setting.T1.beta = makeNLFunc(vsat_beta);
        setting.T2 = setting.T1;

        // sweep rho_11
        {
            JSONValue[string] results; 
            foreach(P; iota(1, 9, 2)) {
                JSONValue[] list;
                foreach(rho_11_dB; linspace(-110.0L, -20.0L, 300)) {
                    Parameter param = setting;
                    param.T1.rho11 = rho_11_dB.dBToVGain;
                    param.T1.P = P;
                    param.T2 = param.T1;

                    auto result = analyze(param);
                    JSONValue res = result.toJSONValue();
                    res["rho_11_dB"] = rho_11_dB;
                    list ~= res;
                }

                results[P.to!string] = JSONValue(list);
            }

            import std.file : write;
            write("sweep_rho_11_%s.json".format(iden), JSONValue(results).toString());
        }

        // sweep back-off
        {
            JSONValue[string] results; 
            foreach(P; iota(1, 9, 2)) {
                JSONValue[] list;
                foreach(ibo_dB; linspace(-5.0L, 25.0L, 300)) {
                    Parameter param = setting;
                    param.T1.g = vsat_alpha / ibo_dB.dBToVGain;
                    param.T1.P = P;
                    param.T2 = param.T1;

                    auto result = analyze(param);
                    JSONValue res = result.toJSONValue();
                    res["ibo_dB"] = ibo_dB;
                    list ~= res;
                }
                results[P.to!string] = JSONValue(list);
            }

            import std.file : write;
            write("sweep_backoff_%s.json".format(iden), JSONValue(results).toString());
        }

        // sweep back-off with IRR = 25 dB
        {
            JSONValue[string] results; 
            foreach(P; iota(1, 9, 2)) {
                JSONValue[] list;
                foreach(ibo_dB; linspace(-5.0L, 25.0L, 300)) {
                    Parameter param = setting;
                    param.T1.g = vsat_alpha / ibo_dB.dBToVGain;
                    param.T1.P = P;
                    param.T1.IRR = 10.0^^(25.0/10);
                    param.T2 = param.T1;

                    auto result = analyze(param);
                    JSONValue res = result.toJSONValue();
                    res["ibo_dB"] = ibo_dB;
                    list ~= res;
                }
                results[P.to!string] = JSONValue(list);
            }

            import std.file : write;
            write("sweep_backoff_%s_withIRR25dB.json".format(iden), JSONValue(results).toString());
        }
    }

    // sweep smoothness factor of Rapp model
    {
        JSONValue[string] results; 
        foreach(P; iota(1, 9, 2)) {
            JSONValue[] list;
            foreach(sf; linspace(0.1L, 6.0L, 300).chain(only(real.infinity))) {
                Parameter param = defaultSetting;
                param.T1.alpha = makeRapp(vsat_alpha, sf);
                param.T1.P = P;
                param.T2 = param.T1;

                auto result = analyze(param);
                JSONValue res = result.toJSONValue();
                res["sf"] = sf;
                list ~= res;
            }
            results[P.to!string] = JSONValue(list);
        }

        import std.file : write;
        write("sweep_sf_Rapp.json", JSONValue(results).toString(JSONOptions.specialFloatLiterals));
    }

    // sweep smoothness factor of Rapp model with IRR = 25dB
    {
        JSONValue[string] results; 
        foreach(P; iota(1, 9, 2)) {
            JSONValue[] list;
            foreach(sf; linspace(0.1L, 6.0L, 300).chain(only(real.infinity))) {
                Parameter param = defaultSetting;
                param.T1.alpha = makeRapp(vsat_alpha, sf);
                param.T1.P = P;
                param.T1.IRR = 10.0^^(25.0/10);
                param.T2 = param.T1;

                auto result = analyze(param);
                JSONValue res = result.toJSONValue();
                res["sf"] = sf;
                list ~= res;
            }
            results[P.to!string] = JSONValue(list);
        }

        import std.file : write;
        write("sweep_sf_Rapp_withIRR25dB.json", JSONValue(results).toString(JSONOptions.specialFloatLiterals));
    }
}
