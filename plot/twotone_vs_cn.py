import figure_env
from matplotlib import pyplot as plt

import numpy as np
import scipy.special as sp

# 数値積分の設定
deint_N = 200  # 標本点数
deint_ta = -5   # 変換後の積分開始点
deint_tb = 5    # 変換後の積分終了点

deint_ts = np.linspace(deint_ta, deint_tb, deint_N)
deint_xs = np.exp(deint_ts - np.exp(-deint_ts))
deint_ws = (1 + np.exp(-deint_ts)) * deint_xs * (deint_tb - deint_ta) / (deint_N-1)
expectX_ws = 2 * deint_xs * np.exp(-deint_xs**2) * deint_ws

def expectX(func):
    vs = np.vectorize(func)(deint_xs)
    return np.sum(vs * expectX_ws)

# P次の係数を計算
def coefP(func, p):
    m = (p-1)/2
    vs = np.vectorize(func)(deint_xs)
    ls = sp.genlaguerre(m, 1)(deint_xs**2) * deint_xs / np.sqrt(m + 1)
    return np.sum(vs * ls * expectX_ws)



def two_tone_test(func, power):
	ts = np.linspace(0, 2*np.pi, 1024, endpoint=False)
	xs = np.cos(ts) * np.sqrt(2*power)
	ys = func(xs)
	freq = np.fft.fft(ys) / len(xs)
	ps = np.abs(freq)**2
	qs = ps[0:len(ps)//2] * 2
	qs[0] = qs[0] / 2
	return qs

def complex_gaussian_test(func, power):
	amp = np.sqrt(power)
	backoffed_func = lambda x: func(x * amp)

	dst = []
	for q in range(5):
		p = q * 2 + 1
		dst.append(coefP(backoffed_func, p)**2)

	return dst

def get_linear_approximation(xs, ys):
	return (ys[1] - ys[0])/(xs[1] - xs[0]) * (xs - xs[0]) + ys[0];


# fun = np.vectorize(lambda x: x/np.sqrt(1+x*x))
fun = np.vectorize(lambda x: x/np.sqrt(np.sqrt(1+x*x*x*x)))
powers = np.logspace(-2.2, 1.2, 300)

two_tone_test_values = np.array(list(map(lambda x: 10*np.log10(two_tone_test(fun, x)), powers)))
complex_gaussian_test_values = np.array(list(map(lambda x: 10*np.log10(complex_gaussian_test(fun, x)), powers)))

powers_dB = 10*np.log10(powers)

fig = plt.figure(figsize=(4*0.80*2, 3*0.9*2))
ax = fig.add_subplot(1, 1, 1)

ax.plot(powers_dB, two_tone_test_values[:,1], label=r"1st, 2-tone", linestyle='dashed', color='C0')
ax.plot(powers_dB, two_tone_test_values[:,3], label=r"3rd, 2-tone", linestyle='dashed', color='C1')
ax.plot(powers_dB, two_tone_test_values[:,5], label=r"5th, 2-tone", linestyle='dashed', color='C4')

# # TwoToneの結果から直線近似
# ax.plot(powers_dB, get_linear_approximation(powers_dB, two_tone_test_values[:,1]), linestyle='dotted', color='C0')
# ax.plot(powers_dB, get_linear_approximation(powers_dB, two_tone_test_values[:,3]), linestyle='dotted', color='C1')
# ax.plot(powers_dB, get_linear_approximation(powers_dB, two_tone_test_values[:,5]), linestyle='dotted', color='C4')


ax.plot(powers_dB, complex_gaussian_test_values[:,0], label=r"1st, OFDM", linestyle='solid', color='C2')
ax.plot(powers_dB, complex_gaussian_test_values[:,1], label=r"3rd, OFDM", linestyle='solid', color='C3')
ax.plot(powers_dB, complex_gaussian_test_values[:,2], label=r"5th, OFDM", linestyle='solid', color='C5')



ax.legend(loc='lower right')
ax.grid(b=True, which='major', linestyle='--')
ax.set_xlabel("Input Power (dB)")
ax.set_ylabel("Output Power (dB)")
ax.set_ylim(-102, 12)
ax.set_xlim(-22, 12)

plt.savefig('twotone_vs_cn.pgf', transparent=True)
plt.savefig('twotone_vs_cn.pdf', transparent=True)

# ts = np.linspace(0, 2*np.pi, 32, endpoint=False)
# xs = np.cos(ts)

# # fun = np.vectorize(lambda x: x)
# # fun = np.vectorize(lambda x: x/np.sqrt(1+x*x))
# # fun = np.vectorize(lambda x: x/np.cbrt(np.sqrt(1+x*x*x*x*x*x)))
# fun = np.vectorize(lambda x: min(x, 0.5))

# ys = fun(xs)

# ps = np.fft.fft(ys) / len(xs)
# qs = ps[0:len(ps)//2] * 2
# qs[0] = qs[0] / 2

# # print(ps)
# # print(qs)

# rs = np.linspace(0, 1, 300)

# approx = np.polynomial.chebyshev.chebval(rs, qs)
# trues = fun(rs)

# plt.plot(rs, approx)
# plt.plot(rs, trues, linestyle="dashed")
# plt.savefig('softlimit_{0}.pdf'.format(len(ts)//2), transparent=True)
