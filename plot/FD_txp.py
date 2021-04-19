import figure_env
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import json
import numpy as np
import math
import zipfile
import scipy.special as sp
import scipy.optimize
import deint
from matplotlib.ticker import FuncFormatter


noise_dBm = -174 + 10*np.log10(20*10**6 * 8) + 4


def loss_func_rapp(x):
    x = 10**(x/20)
    return 10*np.log10(deint.expectX(lambda r: deint.rapp(r, x, 3)**2))

def loss_func_saleh(x):
    x = 10**(x/20)
    return 10*np.log10(deint.expectX(lambda r: deint.saleh(r, x)**2))

def txp_to_ticks_pos(txps, model):
    if model == "Rapp":
        lossfunc = loss_func_rapp
    else:
        lossfunc = loss_func_saleh

    ticks_pos = []
    for txp in txps:
        res = scipy.optimize.minimize_scalar(lambda x: (30-x+lossfunc(x) - txp)**2)
        # print("txp: ", txp, "ibo: ", res.x)
        ticks_pos.append(30 - res.x)

    return ticks_pos



def getSER_EVM_CAP(iden, ctype, txps, iso, loss):
    zipfilename = '../sim_data/{0}.zip'.format(iden)
    filenameBase = iden + '/results_txp_vs_sic/TXP{0}dBm_iso{2}dB_loss{3}dB_{1}_IRR25dB_200/all_result.bin'
    fnamelist = list(map(lambda v: filenameBase.format(v, ctype, iso, loss), txps))

    def extractSICR(data):
        data = data["cancellations"]
        data = list(filter(lambda x: x != 'NaN', data))
        data = np.sort(data)
        data = data[int(0.01*len(data)):int(len(data)*0.99)]
        return -10*np.log10(np.mean(10**(-np.array(data)/10.0)))

    def extractRemainSIPower(data):
        return np.mean(data["RemainPowers"])

    def extractSIPower(data):
        return np.mean(data["SIPowers"])

    def extractSER(data):
        data = data["sers"]
        m = np.average(data)
        if m == 0:
            return 1E-10
        else:
            return m

    def extractEVM(data):
        return np.array(data["evms"]).transpose()


    def extractCAP(data):
        data = np.array(data["evms"])
        c = np.average(2*np.log2(1 + 1.0/data.flatten()))
        return c

    sicrs = figure_env.getListFromZip(zipfilename, fnamelist, extractSICR)
    sers = figure_env.getListFromZip(zipfilename, fnamelist, extractSER)
    evms = figure_env.getListFromZip(zipfilename, fnamelist, extractEVM)
    caps = figure_env.getListFromZip(zipfilename, fnamelist, extractCAP)
    rempowers = figure_env.getListFromZip(zipfilename, fnamelist, extractRemainSIPower)
    sipowers = figure_env.getListFromZip(zipfilename, fnamelist, extractSIPower)

    return [sers, evms, caps, sicrs, 10*np.log10(np.array(sipowers)/np.array(rempowers))]



for NLFuncIden in ["Rapp", "Saleh", "Saleh_noPM"]:
    for loss in [70]:
        for iso in [50]:

            SimZipFilename = "results_Taps_{0}64_10001".format(NLFuncIden)

            sim_txps = range(10, 32, 1)
            sim_backoffs = 30 - np.array(list(sim_txps))
            # sim_results = [
            #     getSER_EVM_CAP(SimZipFilename, "L_LS", sim_txps, iso, loss),
            #     getSER_EVM_CAP(SimZipFilename, "PHPAOnly3_LS", sim_txps, iso, loss),
            #     getSER_EVM_CAP(SimZipFilename, "PHPAOnly5_LS", sim_txps, iso, loss),
            #     getSER_EVM_CAP(SimZipFilename, "PHPAOnly7_LS", sim_txps, iso, loss)
            # ]

            # BER
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)

            # ax.scatter(sim_backoffs, sim_results[0][0], label='$P$=1, sim', color='C0', marker='x')
            # ax.scatter(sim_backoffs, sim_results[1][0], label='$P$=3, sim', color='C1', marker='*')
            # ax.scatter(sim_backoffs, sim_results[2][0], label='$P$=5, sim', color='C2', marker='.')
            # ax.scatter(sim_backoffs, sim_results[3][0], label='$P$=7, sim', color='C3', marker='+')

            with open('../sweep_backoff_{0}.json'.format(NLFuncIden)) as theo_results:
                data = json.load(theo_results)
                xs = [ e["ibo_dB"] for e in data['1'] ]
                ax.plot(xs, [ (e["SER1"] + e["SER2"])/2 for e in data["1"] ], color='C0', linestyle='solid', label='$P$=1, theo')
                ax.plot(xs, [ (e["SER1"] + e["SER2"])/2 for e in data["3"] ], color='C1', linestyle='dashed', label='$P$=3, theo')
                ax.plot(xs, [ (e["SER1"] + e["SER2"])/2 for e in data["5"] ], color='C2', linestyle='dashdot', label='$P$=5, theo')
                ax.plot(xs, [ (e["SER1"] + e["SER2"])/2 for e in data["7"] ], color='C3', linestyle='dotted', label='$P$=7, theo')

            ax.legend(ncol=2, fontsize=12)
            ax.grid(b=True, which='major', linestyle='--')
            ax.set_xlabel("Input Back-Off of the TX Amplifier: $\mathrm{IBO}$ (dB)")
            ax.set_ylabel("Symbol Error Rate: $(\overline{\mathrm{SER}}_{1} + \overline{\mathrm{SER}}_{2})/2$")
            ax.set_yscale("log")
            ax.set_ylim(1E-4, 1E-0)
            ax.set_xlim(-3, 23)
            ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
            ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatterExp10))

            ax2 = ax.twiny()
            ax2.set_xlim(33, 7)

            minor_txp_ticks = range(7, 30, 1)
            ax2.xaxis.set_minor_locator(ticker.FixedLocator((txp_to_ticks_pos(minor_txp_ticks, NLFuncIden))))

            major_txp_ticks = range(10, 30, 5)
            ax2.xaxis.set_major_locator(ticker.FixedLocator((txp_to_ticks_pos(major_txp_ticks, NLFuncIden))))
            ax2.xaxis.set_major_formatter(ticker.FixedFormatter((major_txp_ticks)))
            ax2.tick_params(which='minor', length=3)

            ax2.set_xlabel("Transmission Power [dBm]")
            ax2.plot(np.linspace(10, 25, 100), np.linspace(0, 10, 100), color='#00000000')

            fig.tight_layout()
            fig.savefig('ser_FD_txp_{2}_{0}_{1}.pdf'.format(iso, loss, NLFuncIden), transparent=True)
            fig.savefig('ser_FD_txp_{2}_{0}_{1}.pgf'.format(iso, loss, NLFuncIden), transparent=True)


            # Capacity
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)

            with open('../sweep_backoff_{0}.json'.format(NLFuncIden)) as theo_results:
                data = json.load(theo_results)
                xs = [ e["ibo_dB"] for e in data["1"] ]
                ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["1"] ], color='C0', linestyle='solid', label='FD, $P$=1')
                ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["3"] ], color='C1', linestyle='dashed', label='FD, $P$=3')
                ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["5"] ], color='C2', linestyle='dashdot', label='FD, $P$=5')
                ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["7"] ], color='C3', linestyle='dotted', label='FD, $P$=7')
                ax.plot(xs, [ (e["Rate1_HD"] + e["Rate2_HD"])/2 for e in data["7"] ], color='C4', linestyle='solid', label='HD')

            ax.legend(ncol=2, fontsize=12)
            ax.grid(b=True, which='major', linestyle='--')
            ax.set_xlabel("Input Back-Off of the TX Amplifier: $\mathrm{IBO}$ (dB)")
            ax.set_ylabel("Achievable Rate [bits/s/Hz]")
            ax.set_xlim(-3, 23)
            ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
            ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

            ax2 = ax.twiny()
            ax2.set_xlim(33, 7)

            minor_txp_ticks = range(7, 30, 1)
            ax2.xaxis.set_minor_locator(ticker.FixedLocator((txp_to_ticks_pos(minor_txp_ticks, NLFuncIden))))

            major_txp_ticks = range(10, 30, 5)
            ax2.xaxis.set_major_locator(ticker.FixedLocator((txp_to_ticks_pos(major_txp_ticks, NLFuncIden))))
            ax2.xaxis.set_major_formatter(ticker.FixedFormatter((major_txp_ticks)))
            ax2.tick_params(which='minor', length=3)

            ax2.set_xlabel("Transmission Power [dBm]")
            ax2.plot(np.linspace(10, 25, 100), np.linspace(0, 10, 100), color='#00000000')

            fig.tight_layout()
            fig.savefig('cap_FD_txp_{2}_{0}_{1}.pdf'.format(iso, loss, NLFuncIden), transparent=True)
            fig.savefig('cap_FD_txp_{2}_{0}_{1}.pgf'.format(iso, loss, NLFuncIden), transparent=True)


            # SICR
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)

            # ax.scatter(sim_backoffs, sim_results[0][3], label='$P$=1, sim', color='C0', marker='x')
            # ax.scatter(sim_backoffs, sim_results[1][3], label='$P$=3, sim', color='C1', marker='*')
            # ax.scatter(sim_backoffs, sim_results[2][3], label='$P$=5, sim', color='C2', marker='.')
            # ax.scatter(sim_backoffs, sim_results[3][3], label='$P$=7, sim', color='C3', marker='+')

            with open('../sweep_backoff_{0}.json'.format(NLFuncIden)) as theo_results:
                data = json.load(theo_results)
                xs = [ e["ibo_dB"] for e in data["1"] ]
                ax.plot(xs, [ e["SICR1"] for e in data["1"] ], color='C0', linestyle='solid', label='$P$=1, theo')
                ax.plot(xs, [ e["SICR1"] for e in data["3"] ], color='C1', linestyle='dashed', label='$P$=3, theo')
                ax.plot(xs, [ e["SICR1"] for e in data["5"] ], color='C2', linestyle='dashdot', label='$P$=5, theo')
                ax.plot(xs, [ e["SICR1"] for e in data["7"] ], color='C3', linestyle='dotted', label='$P$=7, theo')

            ax.legend(ncol=2, fontsize=12, loc='lower right')
            ax.grid(b=True, which='major', linestyle='--')
            ax.set_xlabel("Input Back-Off of the TX Amplifier: $\mathrm{IBO}$ (dB)")
            ax.set_ylabel("Cancellation Performance: $\mathrm{SICR}_1$ (dB)")
            ax.set_xlim(-3, 23)
            ax.set_ylim(-3, 63)
            ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
            ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

            ax2 = ax.twiny()
            ax2.set_xlim(33, 7)

            minor_txp_ticks = range(7, 30, 1)
            ax2.xaxis.set_minor_locator(ticker.FixedLocator((txp_to_ticks_pos(minor_txp_ticks, NLFuncIden))))

            major_txp_ticks = range(10, 30, 5)
            ax2.xaxis.set_major_locator(ticker.FixedLocator((txp_to_ticks_pos(major_txp_ticks, NLFuncIden))))
            ax2.xaxis.set_major_formatter(ticker.FixedFormatter((major_txp_ticks)))
            ax2.tick_params(which='minor', length=3)

            ax2.set_xlabel("Transmission Power [dBm]")
            ax2.plot(np.linspace(10, 25, 100), np.linspace(0, 10, 100), color='#00000000')

            fig.tight_layout()
            fig.savefig('sicr_FD_txp_{2}_{0}_{1}.pdf'.format(iso, loss, NLFuncIden), transparent=True)
            fig.savefig('sicr_FD_txp_{2}_{0}_{1}.pgf'.format(iso, loss, NLFuncIden), transparent=True)
            plt.close(fig)


            # NLPowers
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)

            with open('../sweep_backoff_{0}.json'.format(NLFuncIden)) as theo_results:
                data = json.load(theo_results)
                xs = [ e["ibo_dB"] for e in data["1"] ]
                ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][1:-1])) for e in data["1"] ], color='C0', linestyle='solid', label='$P \geq 3$, theo')
                ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][2:-1])) for e in data["1"] ], color='C1', linestyle='dashed', label='$P \geq 5$, theo')
                ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][3:-1])) for e in data["1"] ], color='C2', linestyle='dashdot', label='$P \geq 7$, theo')
                ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][4:-1])) for e in data["1"] ], color='C3', linestyle='dotted', label='$P \geq 9$, theo')

            ax.legend(ncol=2, fontsize=12)
            ax.grid(b=True, which='major', linestyle='--')
            ax.set_xlabel("Input Back-Off of the TX Amplifier: $\mathrm{IBO}$ (dB)")
            ax.set_ylabel("Power (dB)")
            ax.set_xlim(-3, 23)
            ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
            ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

            ax2 = ax.twiny()
            ax2.set_xlim(33, 7)

            minor_txp_ticks = range(7, 30, 1)
            ax2.xaxis.set_minor_locator(ticker.FixedLocator((txp_to_ticks_pos(minor_txp_ticks, NLFuncIden))))

            major_txp_ticks = range(10, 30, 5)
            ax2.xaxis.set_major_locator(ticker.FixedLocator((txp_to_ticks_pos(major_txp_ticks, NLFuncIden))))
            ax2.xaxis.set_major_formatter(ticker.FixedFormatter((major_txp_ticks)))
            ax2.tick_params(which='minor', length=3)

            ax2.set_xlabel("Transmission Power [dBm]")
            ax2.plot(np.linspace(10, 25, 100), np.linspace(0, 10, 100), color='#00000000')

            fig.tight_layout()
            fig.savefig('nlpowers_FD_txp_{2}_{0}_{1}.pdf'.format(iso, loss, NLFuncIden), transparent=True)
            fig.savefig('nlpowers_FD_txp_{2}_{0}_{1}.pgf'.format(iso, loss, NLFuncIden), transparent=True)
            plt.close(fig)
