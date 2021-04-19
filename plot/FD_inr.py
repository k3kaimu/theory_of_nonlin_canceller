import figure_env
import matplotlib.pyplot as plt
import json
import numpy as np
import math
import zipfile
import scipy.special as sp
from matplotlib.ticker import FuncFormatter


noise_dBm = -174 + 10*np.log10(20*10**6 * 8) + 4
ptx_dBm = 23

@np.vectorize
def inr_to_rho11(inr_dB):
    pi_dB = inr_dB + noise_dBm
    return -(ptx_dBm - pi_dB)

def rho11_to_inr(arf_dB):
    return -inr_to_rho11(-arf_dB)


def getSER_EVM_CAP(iden, ctype, inrs, loss):
    zipfilename = '../sim_data/{0}.zip'.format(iden)
    filenameBase = iden + '/results_inr_vs_sic/TXP23dBm_inr{0}dB_loss{2}dB_{1}_IRR25dB_200/all_result.bin'
    fnamelist = list(map(lambda v: filenameBase.format(v, ctype, loss), inrs))

    def extractSICR(data):
        data = data["cancellations"]
        data = list(filter(lambda x: x != 'NaN', data))
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


iso = 50

for NLFuncIden in ["Rapp", "Saleh"]:
    for loss in [70]:

        SimZipFilename = "results_Taps_{0}64_10001".format(NLFuncIden)

        sim_inrs = range(0, 82, 3)
        sim_rho11 = inr_to_rho11(np.array(list(sim_inrs)))
        # sim_results = [
        #     getSER_EVM_CAP(SimZipFilename, "L_LS", sim_inrs, loss),
        #     getSER_EVM_CAP(SimZipFilename, "PHPAOnly3_LS", sim_inrs, loss),
        #     getSER_EVM_CAP(SimZipFilename, "PHPAOnly5_LS", sim_inrs, loss),
        #     getSER_EVM_CAP(SimZipFilename, "PHPAOnly7_LS", sim_inrs, loss)
        # ]

        # BER
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # ax.scatter(sim_rho11, sim_results[0][0], label='$P$=1, sim', color='C0', marker='x')
        # ax.scatter(sim_rho11, sim_results[1][0], label='$P$=3, sim', color='C1', marker='*')
        # ax.scatter(sim_rho11, sim_results[2][0], label='$P$=5, sim', color='C2', marker='.')
        # ax.scatter(sim_rho11, sim_results[3][0], label='$P$=7, sim', color='C3', marker='+')

        with open('../sweep_rho_11_{0}.json'.format(NLFuncIden)) as theo_results:
            data = json.load(theo_results)
            xs = [ e["rho_11_dB"] for e in data['1'] ]
            ax.plot(xs, [ (e["SER1"] + e["SER2"])/2 for e in data["1"] ], color='C0', linestyle='solid', label='$P$=1, theo')
            ax.plot(xs, [ (e["SER1"] + e["SER2"])/2 for e in data["3"] ], color='C1', linestyle='dashed', label='$P$=3, theo')
            ax.plot(xs, [ (e["SER1"] + e["SER2"])/2 for e in data["5"] ], color='C2', linestyle='dashdot', label='$P$=5, theo')
            ax.plot(xs, [ (e["SER1"] + e["SER2"])/2 for e in data["7"] ], color='C3', linestyle='dotted', label='$P$=7, theo')

        ax.legend(ncol=2, fontsize=12)

        ax.grid(b=True, which='major', linestyle='--')
        ax.set_xlabel(r"Sum of Propagation Loss and RF Cancellation: $\rho_{11}^2$ (dB)")
        ax.set_ylabel("Symbol Error Rate: $(\overline{\mathrm{SER}}_{1} + \overline{\mathrm{SER}}_{2})/2$")
        ax.set_yscale("log")
        ax.set_xlim(-105, -25)
        ax2 = ax.twiny()
        ax2.set_xlim(rho11_to_inr(-105), rho11_to_inr(-25))
        ax2.plot([rho11_to_inr(-25), rho11_to_inr(-105)], [1E-3, 1E-3], color="#00000000")
        ax2.set_xlabel("INR (dB)")

        ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatterExp10))
        ax2.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

        fig.tight_layout()
        fig.savefig('ser_FD_inr_{1}_{0}.pdf'.format(loss, NLFuncIden), transparent=True)
        fig.savefig('ser_FD_inr_{1}_{0}.pgf'.format(loss, NLFuncIden), transparent=True)


        # Capacity
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        with open('../sweep_rho_11_{0}.json'.format(NLFuncIden)) as theo_results:
            data = json.load(theo_results)
            xs = [ e["rho_11_dB"] for e in data["1"] ]
            ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["1"] ], color='C0', linestyle='solid', label='FD, $P$=1')
            ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["3"] ], color='C1', linestyle='dashed', label='FD, $P$=3')
            ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["5"] ], color='C2', linestyle='dashdot', label='FD, $P$=5')
            ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["7"] ], color='C3', linestyle='dotted', label='FD, $P$=7')
            ax.plot(xs, [ (e["Rate1_HD"] + e["Rate2_HD"])/2 for e in data["7"] ], color='C4', linestyle='solid', label='HD')


        ax.legend(ncol=2, fontsize=12)

        ax.grid(b=True, which='major', linestyle='--')
        ax.set_xlabel(r"Sum of Propagation Loss and RF Cancellation: $\rho_{11}^2$ (dB)")
        ax.set_ylabel("Achievable Rate [bits/s/Hz]")
        # ax.set_yscale("log")
        ax.set_xlim(-105, -25)
        ax2 = ax.twiny()
        ax2.set_xlim(rho11_to_inr(-105), rho11_to_inr(-25))
        ax2.plot([rho11_to_inr(-25), rho11_to_inr(-105)], [1E-3, 1E-3], color="#00000000")
        ax2.set_xlabel("INR (dB)")

        ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax2.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

        fig.tight_layout()
        fig.savefig('cap_FD_inr_{1}_{0}.pdf'.format(loss, NLFuncIden), transparent=True)
        fig.savefig('cap_FD_inr_{1}_{0}.pgf'.format(loss, NLFuncIden), transparent=True)


        # SICR
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # ax.scatter(sfvals, getBER("Nop_X", sfvals), label='0', color='C0', marker='x')
        # ax.scatter(sim_rho11, sim_results[0][4], label='$P$=1, sim', color='C0', marker='x')
        # ax.scatter(sim_rho11, sim_results[1][4], label='$P$=3, sim', color='C1', marker='*')
        # ax.scatter(sim_rho11, sim_results[2][4], label='$P$=5, sim', color='C2', marker='.')
        # ax.scatter(sim_rho11, sim_results[3][4], label='$P$=7, sim', color='C3', marker='+')

        with open('../sweep_rho_11_{0}.json'.format(NLFuncIden)) as theo_results:
            data = json.load(theo_results)
            xs = [ e["rho_11_dB"] for e in data["1"] ]
            ax.plot(xs, [ e["SICR1"] for e in data["1"] ], color='C0', linestyle='solid', label='$P$=1, theo')
            ax.plot(xs, [ e["SICR1"] for e in data["3"] ], color='C1', linestyle='dashed', label='$P$=3, theo')
            ax.plot(xs, [ e["SICR1"] for e in data["5"] ], color='C2', linestyle='dashdot', label='$P$=5, theo')
            ax.plot(xs, [ e["SICR1"] for e in data["7"] ], color='C3', linestyle='dotted', label='$P$=7, theo')

        # rho_11 = np.linspace(-110, -20, 300)
        # rho_11 = 10**(rho_11/10)
        # ax.plot()
        inrs_for_upper = np.linspace(0, 90)
        rho_11_for_upper = inr_to_rho11(inrs_for_upper)
        sicr_upper = 10*np.log10(10**(inrs_for_upper/10) + 1)
        ax.plot(rho_11_for_upper, sicr_upper, color='C4', linestyle='dotted', label=r'Upper bound: $\mathrm{INR}+1$')
        print(rho_11_for_upper)
        print(sicr_upper)
        print("")

        ax.legend(ncol=2, fontsize=12)
        ax.grid(b=True, which='major', linestyle='--')
        ax.set_xlabel(r"Sum of Propagation Loss and RF Cancellation: $\rho_{11}^2$ (dB)")
        ax.set_ylabel("Cancellation Performance: $\mathrm{SICR}_1$ (dB)")
        ax.set_xlim(-105, -25)
        ax.set_ylim(-3, 63)

        ax2 = ax.twiny()
        ax2.set_xlim(rho11_to_inr(-105), rho11_to_inr(-25))
        ax2.plot([rho11_to_inr(-25), rho11_to_inr(-105)], [1E-3, 1E-3], color="#00000000")
        ax2.set_xlabel("$\mathrm{INR}_1$ (dB)")

        ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax2.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

        fig.tight_layout()
        fig.savefig('sicr_FD_inr_{2}_{0}_{1}.pdf'.format(iso, loss, NLFuncIden), transparent=True)
        fig.savefig('sicr_FD_inr_{2}_{0}_{1}.pgf'.format(iso, loss, NLFuncIden), transparent=True)
        plt.close(fig)


        # NLPowers
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        with open('../sweep_rho_11_{0}.json'.format(NLFuncIden)) as theo_results:
            data = json.load(theo_results)
            xs = [ e["rho_11_dB"] for e in data["1"] ]
            ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][1:-1])) for e in data["1"] ], color='C0', linestyle='solid', label='$P \geq 3$, theo')
            ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][2:-1])) for e in data["1"] ], color='C1', linestyle='dashed', label='$P \geq 5$, theo')
            ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][3:-1])) for e in data["1"] ], color='C2', linestyle='dashdot', label='$P \geq 7$, theo')
            ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][4:-1])) for e in data["1"] ], color='C3', linestyle='dotted', label='$P \geq 9$, theo')

        ax.legend(ncol=2, fontsize=12)
        ax.grid(b=True, which='major', linestyle='--')
        ax.set_xlabel(r"Sum of Propagation Loss and RF Cancellation: $\rho_{11}^2$ (dB)")
        ax.set_ylabel("Power (dB)")
        ax.set_xlim(-105, -25)
        # ax.set_ylim(10**-6.2, 10**0.2)
        ax2 = ax.twiny()
        ax2.set_xlim(rho11_to_inr(-105), rho11_to_inr(-25))
        ax2.plot([rho11_to_inr(-25), rho11_to_inr(-105)], [1E-3, 1E-3], color="#00000000")
        ax2.set_xlabel("INR (dB)")

        ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax2.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

        fig.tight_layout()
        fig.savefig('nlpowers_FD_inr_{2}_{0}_{1}.pdf'.format(iso, loss, NLFuncIden), transparent=True)
        fig.savefig('nlpowers_FD_inr_{2}_{0}_{1}.pgf'.format(iso, loss, NLFuncIden), transparent=True)
        plt.close(fig)


        # Residual powers
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        with open('../sweep_rho_11_{0}.json'.format(NLFuncIden)) as theo_results:
            data = json.load(theo_results)
            xs = [ e["rho_11_dB"] for e in data["1"] ]
            ax.plot(xs, [ 10*np.log10(e["ResPower1"] + e["N_NL1"])+30 for e in data["1"] ], color='C0', linestyle='solid', label='$P=1$')
            ax.plot(xs, [ 10*np.log10(e["ResPower1"] + e["N_NL1"])+30 for e in data["3"] ], color='C1', linestyle='dashed', label='$P=3$')
            ax.plot(xs, [ 10*np.log10(e["ResPower1"] + e["N_NL1"])+30 for e in data["5"] ], color='C2', linestyle='dashdot', label='$P=5$')
            ax.plot(xs, [ 10*np.log10(e["ResPower1"] + e["N_NL1"])+30 for e in data["7"] ], color='C3', linestyle='dotted', label='$P=7$')
            ax.plot([-110, -20], [noise_dBm, noise_dBm], color='C4', linestyle='solid', label=r'$\overline{N}_\mathrm{thermal}$')

        ax.legend(ncol=2, fontsize=12)
        ax.grid(b=True, which='major', linestyle='--')
        ax.set_xlabel(r"Sum of Propagation Loss and RF Cancellation: $\rho_{11}^2$ (dB)")
        ax.set_ylabel(r"Residual Self-Interference Power: ${\smash{\overline{I}}\vphantom{I}}^R_\mathrm{11} + \overline{N}_\mathrm{NL,1}$ [dBm]")
        ax.set_xlim(-105, -25)
        ax.set_ylim(-105, -5)
        # ax.set_ylim(10**-6.2, 10**0.2)
        ax2 = ax.twiny()
        ax2.set_xlim(rho11_to_inr(-105), rho11_to_inr(-25))
        ax2.plot([rho11_to_inr(-25), rho11_to_inr(-105)], [1E-3, 1E-3], color="#00000000")
        ax2.set_xlabel("INR (dB)")

        ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax2.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

        fig.tight_layout()
        fig.savefig('residualSI_FD_inr_{2}_{0}_{1}.pdf'.format(iso, loss, NLFuncIden), transparent=True)
        fig.savefig('residualSI_FD_inr_{2}_{0}_{1}.pgf'.format(iso, loss, NLFuncIden), transparent=True)
        plt.close(fig)

