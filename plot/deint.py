import numpy as np 
import scipy.special as sp

# 数値積分の設定
deint_N = 201   # 標本点数
deint_ta = -5   # 変換後の積分開始点
deint_tb = 5    # 変換後の積分終了点

deint_ts = np.linspace(deint_ta, deint_tb, deint_N)
deint_xs = np.exp(deint_ts - np.exp(-deint_ts))
deint_ws = (1 + np.exp(-deint_ts)) * deint_xs * (deint_tb - deint_ta) / (deint_N-1)
expectX_ws = 2 * deint_xs * np.exp(-deint_xs**2) * deint_ws

def expectX(func):
    vs = func(deint_xs)
    return np.sum(vs * expectX_ws)

def rapp(r, v, s):
    return r / (1 + (r/v)**(2*s) )**(1.0/(2*s))

def saleh(r, v):
    return r / (1 + (r/(2*v))**2 )
