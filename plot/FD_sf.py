import figure_env
import matplotlib.pyplot as plt
import json
import numpy as np
import math
import zipfile
import scipy.special as sp
from matplotlib.ticker import FuncFormatter


noise_dBm = -174 + 10*np.log10(20*10**6 * 8) + 4


def getSER_EVM_CAP(iden, ctype, values, iso, loss):
    zipfilename = '../sim_data/{0}.zip'.format(iden)
    filenameBase = iden + '/results_sf_vs_sic/TXP23dBm_iso{2}dB_loss{3}dB_{1}_IRR25dB_200_sf{0}/all_result.bin'
    fnamelist = list(map(lambda v: filenameBase.format(v, ctype, iso, loss), values))

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


for loss in [70]:
    for iso in [50]:
        # sfvals = range(0.1, 5.1, 0.1)
        sim_sfvals_10 = np.array(list(range(2, 51, 2)))

        SimZipFilename = "results_Taps_Rapp64_10001"

        # sim_results = [
        #     getSER_EVM_CAP(SimZipFilename, "L_LS", sim_sfvals_10, iso, loss),
        #     getSER_EVM_CAP(SimZipFilename, "PHPAOnly3_LS", sim_sfvals_10, iso, loss),
        #     getSER_EVM_CAP(SimZipFilename, "PHPAOnly5_LS", sim_sfvals_10, iso, loss),
        #     getSER_EVM_CAP(SimZipFilename, "PHPAOnly7_LS", sim_sfvals_10, iso, loss)
        # ]

        # BER
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # ax.scatter(sim_sfvals_10 / 10.0, sim_results[0][0], label='$P$=1, sim', color='C0', marker='x')
        # ax.scatter(sim_sfvals_10 / 10.0, sim_results[1][0], label='$P$=3, sim', color='C1', marker='*')
        # ax.scatter(sim_sfvals_10 / 10.0, sim_results[2][0], label='$P$=5, sim', color='C2', marker='.')
        # ax.scatter(sim_sfvals_10 / 10.0, sim_results[3][0], label='$P$=7, sim', color='C3', marker='+')

        def avgSER(e):
            return (e["SER1"] + e["SER2"])/2

        with open('../sweep_sf_Rapp.json') as theo_results:
            data = json.load(theo_results)
            xs = [ e["sf"] for e in data["1"][:-1] ]
            ax.plot(xs, [ avgSER(e) for e in data["1"][:-1] ], color='C0', linestyle='solid', label='$P$=1, theo')
            ax.plot(xs, [ avgSER(e) for e in data["3"][:-1] ], color='C1', linestyle='dashed', label='$P$=3, theo')
            ax.plot(xs, [ avgSER(e) for e in data["5"][:-1] ], color='C2', linestyle='dashdot', label='$P$=5, theo')
            ax.plot(xs, [ avgSER(e) for e in data["7"][:-1] ], color='C3', linestyle='dotted', label='$P$=7, theo')
            # ax.plot(xs, [ (e["SER1_HD"] + e["SER2_HD"])/2 for e in data["1"] ], color='C4', linestyle='solid', label='HD, theo')

            pos_offset = [0.8, 0.8, 1, 0.7]
            for p in [1, 3, 5, 7]:
                infser = avgSER(data[str(p)][-1])
                print("INF SER={} at P={}".format(infser, p))
                ax.annotate(r'$P$={}'.format(p), xy=[5.2, infser], xytext=[5.4, infser * pos_offset[p//2] ], color='C{}'.format( (p-1)//2 ))
                ax.annotate('',
                        xy=[5.2, infser],
                        xytext=[5.4, infser],
                        arrowprops=dict(shrink=0, width=1, headwidth=8, 
                                    headlength=10, connectionstyle='arc3',
                                    facecolor='C{}'.format(p//2), edgecolor='C{}'.format(p//2)))

        ax.legend(ncol=2, fontsize=12)
        ax.grid(b=True, which='major', linestyle='--')
        ax.set_xlabel("Smoothness Factor of the TX Amplifier")
        ax.set_ylabel("Symbol Error Rate: $(\overline{\mathrm{SER}}_{1} + \overline{\mathrm{SER}}_{2})/2$")
        ax.set_yscale("log")
        ax.set_xlim(0, 5.2)
        ax.set_ylim(1E-5, 1E-0)
        ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatterExp10))

        fig.tight_layout()
        fig.savefig('ser_FD_sf_Rapp_{0}_{1}.pdf'.format(iso, loss), transparent=True)
        fig.savefig('ser_FD_sf_Rapp_{0}_{1}.pgf'.format(iso, loss), transparent=True)
        plt.close(fig)


        # Capacity
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        with open('../sweep_sf_Rapp.json') as theo_results:
            data = json.load(theo_results)
            xs = [ e["sf"] for e in data["1"][:-1] ]
            ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["1"][:-1] ], color='C0', linestyle='solid', label='FD, $P$=1')
            ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["3"][:-1] ], color='C1', linestyle='dashed', label='FD, $P$=3')
            ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["5"][:-1] ], color='C2', linestyle='dashdot', label='FD, $P$=5')
            ax.plot(xs, [ e["Rate1"] + e["Rate2"] for e in data["7"][:-1] ], color='C3', linestyle='dotted', label='FD, $P$=7')
            ax.plot(xs, [ (e["Rate1_HD"] + e["Rate2_HD"])/2 for e in data["1"][:-1] ], color='C4', linestyle='solid', label='HD')

        ax.legend(ncol=2, fontsize=12)
        ax.grid(b=True, which='major', linestyle='--')
        ax.set_xlabel("Smoothness Factor of the TX Amplifier")
        ax.set_ylabel("Achievable Rate [bits/s/Hz]")
        ax.set_xlim(0, 5.2)
        ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

        fig.tight_layout()
        fig.savefig('cap_FD_sf_Rapp_{0}_{1}.pdf'.format(iso, loss), transparent=True)
        fig.savefig('cap_FD_sf_Rapp_{0}_{1}.pgf'.format(iso, loss), transparent=True)
        plt.close(fig)


        # SICR
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        # ax.scatter(sim_sfvals_10 / 10.0, sim_results[0][4], label='$P$=1, sim', color='C0', marker='x')
        # ax.scatter(sim_sfvals_10 / 10.0, sim_results[1][4], label='$P$=3, sim', color='C1', marker='*')
        # ax.scatter(sim_sfvals_10 / 10.0, sim_results[2][4], label='$P$=5, sim', color='C2', marker='.')
        # ax.scatter(sim_sfvals_10 / 10.0, sim_results[3][4], label='$P$=7, sim', color='C3', marker='+')

        with open('../sweep_sf_Rapp.json') as theo_results:
            data = json.load(theo_results)
            xs = [ e["sf"] for e in data["1"][:-1] ]
            ax.plot(xs, [ e["SICR1"] for e in data["1"][:-1] ], color='C0', linestyle='solid', label='$P$=1, theo')
            ax.plot(xs, [ e["SICR1"] for e in data["3"][:-1] ], color='C1', linestyle='dashed', label='$P$=3, theo')
            ax.plot(xs, [ e["SICR1"] for e in data["5"][:-1] ], color='C2', linestyle='dashdot', label='$P$=5, theo')
            ax.plot(xs, [ e["SICR1"] for e in data["7"][:-1] ], color='C3', linestyle='dotted', label='$P$=7, theo')

            pos_offset = [-1.5, -0.5, -2, 0]
            for p in [1, 3, 5, 7]:
                infval = data[str(p)][-1]["SICR1"]
                ax.annotate(r'$P$={}'.format(p), xy=[5.2, infval], xytext=[5.4, infval + pos_offset[p//2] ], color='C{}'.format( (p-1)//2 ))
                ax.annotate('',
                        xy=[5.2, infval],
                        xytext=[5.4, infval],
                        arrowprops=dict(shrink=0, width=1, headwidth=8, 
                                    headlength=10, connectionstyle='arc3',
                                    facecolor='C{}'.format(p//2), edgecolor='C{}'.format(p//2)))

        ax.legend(ncol=2, fontsize=12)
        ax.grid(b=True, which='major', linestyle='--')
        ax.set_xlabel("Smoothness Factor of the TX Amplifier")
        ax.set_ylabel("Cancellation Performance: $\mathrm{SICR}_1$ (dB)")
        ax.set_xlim(0, 5.2)
        ax.set_ylim(-3, 63)
        ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

        fig.tight_layout()
        fig.savefig('sicr_FD_sf_Rapp_{0}_{1}.pdf'.format(iso, loss), transparent=True)
        fig.savefig('sicr_FD_sf_Rapp_{0}_{1}.pgf'.format(iso, loss), transparent=True)
        plt.close(fig)


        # NLPowers
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        with open('../sweep_sf_Rapp.json') as theo_results:
            data = json.load(theo_results)
            xs = [ e["sf"] for e in data["1"][:-1] ]
            ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][1:-1])) for e in data["1"][:-1] ], color='C0', linestyle='solid', label='$P \geq 3$, theo')
            ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][2:-1])) for e in data["1"][:-1] ], color='C1', linestyle='dashed', label='$P \geq 5$, theo')
            ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][3:-1])) for e in data["1"][:-1] ], color='C2', linestyle='dashdot', label='$P \geq 7$, theo')
            ax.plot(xs, [ 10*np.log10(np.sum(e["NLPowers1"][4:-1])) for e in data["1"][:-1] ], color='C3', linestyle='dotted', label='$P \geq 9$, theo')

        ax.legend(ncol=2, fontsize=12)
        ax.grid(b=True, which='major', linestyle='--')
        ax.set_xlabel("Smoothness Factor of the TX Amplifier")
        ax.set_ylabel(r"Sum of Nonlinear Powers: $|\alpha_{P}|^2+|\alpha_{P+2}|^2+\cdots$ [dBm]")
        ax.set_xlim(0, 5.2)

        fig.tight_layout()
        fig.savefig('nlpowers_FD_sf_Rapp_{0}_{1}.pdf'.format(iso, loss), transparent=True)
        fig.savefig('nlpowers_FD_sf_Rapp_{0}_{1}.pgf'.format(iso, loss), transparent=True)
        plt.close(fig)


        # Residual powers
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        with open('../sweep_sf_Rapp.json') as theo_results:
            data = json.load(theo_results)
            xs = [ e["sf"] for e in data["1"][10:-1] ]
            ax.plot(xs, [ 10*np.log10(e["ResPower1"])+30 for e in data["1"][10:-1] ], color='C0', linestyle='solid', label='$P=1$')
            ax.plot(xs, [ 10*np.log10(e["ResPower1"])+30 for e in data["3"][10:-1] ], color='C1', linestyle='dashed', label='$P=3$')
            ax.plot(xs, [ 10*np.log10(e["ResPower1"])+30 for e in data["5"][10:-1] ], color='C2', linestyle='dashdot', label='$P=5$')
            ax.plot(xs, [ 10*np.log10(e["ResPower1"])+30 for e in data["7"][10:-1] ], color='C3', linestyle='dotted', label='$P=7$')
            ax.plot([-1, 6], [noise_dBm, noise_dBm], color='C4', linestyle='solid', label=r'$\overline{N}_\mathrm{thermal}$')

            pos_offset = [+1, -3, +1, -3]
            for p in [1, 3, 5, 7]:
                infval = 10*np.log10(data[str(p)][-1]["ResPower1"])+30
                ax.annotate(r'$P$={}'.format(p), xy=[5.2, infval], xytext=[5.2, infval + pos_offset[p//2] ], color='C{}'.format( (p-1)//2 ))
                ax.annotate('',
                        xy=[5.2, infval],
                        xytext=[5.4, infval],
                        arrowprops=dict(shrink=0, width=1, headwidth=8, 
                                    headlength=10, connectionstyle='arc3',
                                    facecolor='C{}'.format(p//2), edgecolor='C{}'.format(p//2)))

        ax.legend(ncol=2, fontsize=12)
        ax.grid(b=True, which='major', linestyle='--')
        ax.set_xlabel("Smoothness Factor of the TX Amplifier")
        ax.set_ylabel(r"Residual Self-Interference Power: ${\smash{\overline{I}}\vphantom{I}}^{\mathrm{R}}_\mathrm{11}$ [dBm]")
        ax.set_xlim(0, 5.2)
        ax.set_ylim(-105, -35)
        ax.xaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))
        ax.yaxis.set_major_formatter(FuncFormatter(figure_env.sansmathFormatter))

        fig.tight_layout()
        fig.savefig('residualSI_FD_sf_Rapp_{0}_{1}.pdf'.format(iso, loss), transparent=True)
        fig.savefig('residualSI_FD_sf_Rapp_{0}_{1}.pgf'.format(iso, loss), transparent=True)
        plt.close(fig)
