import matplotlib as mpl

mpl.use("pgf")
pgf_with_pdflatex = {
    "pgf.texsystem": "pdflatex",
    # "font.family": "serif", # use serif/main font for text elements
    # "text.usetex": True,    # use inline math for ticks
    "pgf.rcfonts": False,   # don't setup fonts from rc parameters
    "pgf.preamble": [
         r"\usepackage[utf8x]{inputenc}",
         r"\usepackage[T1]{fontenc}",
         r"\usepackage{amsmath}",
         r"\usepackage[T1]{sansmath}"
         ]
}
mpl.rcParams.update(pgf_with_pdflatex)


import matplotlib.pyplot as plt
import json
import numpy as np
import math
from matplotlib import rc
import re
import os
import zipfile
import msgpack


# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)
# plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
plt.rcParams['figure.figsize'] = (8*0.87, 6*0.87)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['font.size'] = 14
# plt.rcParams['font.size'] = 16
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['lines.markersize'] = 8
plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

plt.locator_params(axis='y',nbins=6)


plt.gca().yaxis.set_tick_params(which='both', direction='in',bottom=True, top=True, left=True, right=True)



# TeXのauxを読み込み，文献の引用番号を取得します．
def load_cite(path):
    regex = re.compile(r'\{(.+?)\}')
    dst = {}
    with open(path) as auxfile:
        for line in auxfile.readlines():
            if line.startswith('\\bibcite'):
                result = regex.findall(line)
                dst[result[0]] = result[1]

    return dst

# pathにある.auxファイルを列挙します
def all_aux(path):
    return map(lambda s: os.path.join(path, s), filter(lambda s: s.endswith(".aux"), os.listdir(path)))

# pathにある.auxファイルのうち，更新時間が最も新しいものを取得します
def most_new_aux(path):
    auxs = all_aux(path)
    return max(auxs, key=lambda s: os.path.getmtime(s))

def getListFromZip(zipfname, fnamelist, func):
    dst = []
    with zipfile.ZipFile(zipfname) as results:
        for fname in fnamelist:
            with results.open(fname, 'r') as data_file:
                if fname.endswith(".json"):
                    data = json.load(data_file)
                else:
                    data = msgpack.unpackb(data_file.read(), raw=False)
                dst.append(func(data))

    return dst

def getMedianCancellationsFromZip(zipfname, fnamelist):

    def extractFunc(data):
        data = data["cancellations"]
        data = list(filter(lambda x: x != 'NaN', data))
        return np.median(data)

    return getListFromZip(zipfname, fnamelist, extractFunc)


def getMeanCancellationsFromZip(zipfname, fnamelist):

    def extractFunc(data):
        data = data["cancellations"]
        data = list(filter(lambda x: x != 'NaN', data))
        return -10*np.log10(np.mean(10**(-np.array(data)/10.0)))

    return getListFromZip(zipfname, fnamelist, extractFunc)


def getAllCancellationsFromZip(zipfname, fnamelist):

    def extractFunc(data):
        data = data["cancellations"]
        data = list(filter(lambda x: x != 'NaN', data))
        return data

    return getListFromZip(zipfname, fnamelist, extractFunc)


def sansmathFormatter(x, pos):
    return r'{\sansmath $%.0f$}' % x

def sansmathFormatterExp10(x, pos):
    s = math.log10(x)
    return r'{\sansmath $10^{%.0f}$}' % s
