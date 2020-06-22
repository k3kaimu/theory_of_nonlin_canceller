module utils;

import std.complex;
import std.math;
import std.range;
import std.traits;

import deint;


// factorial
ulong fact(uint n)
{
    ulong dst = 1;
    foreach(i; 2 .. n + 1)
        dst *= i;

    return dst;
}


/// binomial coefficient
C binom(C)(C x, size_t k)
{
    C num = 1;
    C den = 1;
    foreach(i; 0 .. k) {
        num *= x - i;
        den *= k - i;
    }

    return num / den;
}

//
unittest
{
    assert(binom(5, 0) == 1);
    assert(binom(5, 1) == 5);
    assert(binom(5, 2) == 10);
    assert(binom(5, 3) == 10);
    assert(binom(5, 4) == 5);
    assert(binom(5, 5) == 1);
}


/// generalized Laguerre polynomial
C laguerre(C)(C x, int n, C a = 0)
{
    C sum = 0;
    foreach(i; 0 .. n + 1)
        sum += (-1) ^^ i * binom(n + a, n - i) * x ^^ i / fact(i);

    return sum;
}

//
unittest
{
    import std.math : approxEqual;

    // L_0^a(x) = 1
    auto l0a = delegate(real a, real x) { return 1.0L; };

    // L_1^a(x) = -x + a + 1
    auto l1a = delegate(real a, real x) { return -x + a + 1; };

    // L_2^a(x) = x^2/2 - (a+2)x + (a+2)(a+1)/2
    auto l2a = delegate(real a, real x) {
        return x ^^ 2 / 2 - (a + 2) * x + (a + 2) * (a + 1) / 2;
    };

    // L_3^a(x) = -x^3/6 + (a+3)x^2/2 - (a+2)(a+3)x/2 + (a+1)(a+2)(a+3)/6
    auto l3a = delegate(real a, real x) {
        return -x ^^ 3 / 6 + (a + 3) * x ^^ 2 / 2 - (a + 2) * (a + 3) * x / 2 + (
                a + 1) * (a + 2) * (a + 3) / 6;
    };

    auto testXs = [0, 0.1, 0.5, 1, 2, 5, 10];
    auto testAs = [0, 1, 2, 3, 4, 5];

    foreach(x; testXs)
        foreach(a; testAs)
        {
            foreach(int i, f; [l0a, l1a, l2a, l3a])
                assert(laguerre(x, i, a).approxEqual(f(a, x)));
        }
}


/// orthonormal basis function
C psi(C)(uint p, C x)
    in(p%2 == 1)
    in(p > 0)
{
    int m = p/2;

    alias F = typeof(abs(x));
    immutable F normCoef = 1/sqrt(m + 1.0L);
    return (-1) ^^ m * normCoef * x * laguerre(x.sqAbs, m, 1);
}


/// Ideal predistorted amplifier
C limiter(C)(C x, real vs)
{
    auto r = x.abs;
    if(r < vs) {
        return x;
    }else{
        return C(vs * (x/r));
    }
}


/// rapp model
C rapp(C)(C x, real vs, real s)
{
    if(s == real.infinity)
        return limiter(x, vs);

    auto r = x.abs;
    auto u = r / vs;
    auto den = (u^^(2*s) + 1)^^(1.0 / s / 2);
    return C(r / den * (x/r));
}


// /// Saleh model
// C saleh(C)(C x, real vs)
// {
//     immutable r = std.complex.abs(input),
//               u = input / r;    // unit vector

//     // rが小さすぎるときに，単位ベクトルが発散するのを防ぐ
//     if(r <= 1E-6) {
//         output = input * _g;
//     } else {
//         output = _vs * normalized_saleh(r / _vs) * u;
//     }
// }


/**
飽和電圧1，小信号ゲイン1のsalehモデル
*/
Complex!F normalized_saleh(F)(F r, F phisat)
{
    immutable F aa = 1.0,
                ba = aa^^2 / 4,   // = 1.16402521
                ap = 2 * phisat,
                bp = ba;

    return Complex!F(aa * r / (1+ba*r^^2) * std.complex.expi(ap*r^^2/(1+bp*r^^2)));
}



/// convert [dBm] to [V] on 1-ohm
R dBmToVolt(R)(R dbm)
if(isFloatingPoint!R)
{
    return 10.0L^^((dbm - 30) / 20);
}


/// convert dB to gain of voltage
R dBToVGain(R)(R db)
if(isFloatingPoint!R)
{
    return 10.0L^^(db / 20);
}



auto linspace(R)(R x, R y, size_t N)
if(isFloatingPoint!R)
{
    R[] rets;
    R h = (y - x) / (N - 1);
    foreach(i; 0 .. N)
        rets ~= h * i + x;

    return rets;
}


/**
S_x(f)のn次畳込みしたスペクトルを計算します
S_x(f)は[-1/2, 1/2]で1，それ以外で0です．
*/
real nTimesConvRectSpectrum(real f, size_t n)
in{
    assert(n >= 1);
}
body {
    real powWithMax0(real x, long p)
    {
        if(x > 0)
            return x^^p;
        else
            return 0;
    }

    real dst = 0;
    foreach(k; 0 .. n+1)
        dst += (-1.0L)^^k * binom(n, k) * powWithMax0(abs(f) - k + n/2.0L, n-1);

    return dst / fact(cast(int)(n-1));
}

unittest 
{
    assert(nTimesConvRectSpectrum(0, 1).approxEqual(1));
    assert(nTimesConvRectSpectrum(0.499, 1).approxEqual(1));
    assert(nTimesConvRectSpectrum(0.501, 1).approxEqual(0));
    assert(nTimesConvRectSpectrum(-0.499, 1).approxEqual(1));
    assert(nTimesConvRectSpectrum(-0.501, 1).approxEqual(0));
    
    foreach(f; [-2, -1.501, -1.499, -1.3, -1, -0.501, -0.499, -0.3, -0.2, 0,
                +2, +1.501, +1.499, +1.3, +1, +0.501, +0.499, +0.3, +0.2])
    {
        real ans;
        if(abs(f) <= 0.5)
            ans = 3.0L/4 - f^^2;
        else if(abs(f) > 0.5 && abs(f) < 1.5)
            ans = 0.5*(-abs(f) + 1.5)^^2;
        else
            ans = 0;

        assert(ans.approxEqual(nTimesConvRectSpectrum(f, 3)));
    }
}


// Theoretical symbol error rate of M-QAM on AWGN channel
real qamSER(real snr, uint M)
{
    import std.mathspecial : erfc;

    immutable a = 1 - 1/sqrt(M * 1.0L);
    immutable b = 3.0L / (2 * (M - 1));
    immutable p = a * erfc(sqrt(b * snr));
    return 1 - (1 - p)^^2;
}





version (unittest)
{
    package bool approxEqualC(C1, C2)(C1 x, C2 y)
            // if((isComplex!C1 || isFloatingPoint!C1) && (isComplex!C2 || isFloatingPoint!C2))
    {
        static if(isFloatingPoint!C1)
            return .approxEqualC(complex(x, 0), y);
        else
        {
            static if(isFloatingPoint!C2)
                return .approxEqualC(x, complex(y, 0));
            else
            {
                return std.math.approxEqual(x.re, y.re) && std.math.approxEqual(x.im, y.im);
            }
        }
    }
}


// compute integ.integrate(func) on multi-threads
typeof(NumInt!(X, W).init.integrate(Func.init)) parallelIntegrate(X, W, Func)(ref NumInt!(X, W) integ, Func func)
{
    import std.parallelism;
    import core.atomic;

    alias Ret = typeof(return);

    shared(Ret) sum = 0;
    foreach(xw; zip(integ.xs, integ.ws).parallel)
        sum.atomicOp!"+="(func(xw[0]) * xw[1]);
    
    return sum;
}