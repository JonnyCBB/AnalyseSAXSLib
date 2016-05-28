""" Class dealing with all of the analysis of the multiple compounds used in
the SAXS analysis experiments
"""
from RadioProtectant import Compound
import pandas as pd
import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('ggplot')

class ComparisonAnalysis(object):
    """Analysis of several compounds for SAXS analysis.
    """
    # ----------------------------------------------------------------------- #
    #                         CONSTRUCTOR METHOD                              #
    # ----------------------------------------------------------------------- #
    def __init__(self, compound_list, runs_per_conc=3, buffer_sub=True,
                 average_type="mean", crop_start=1, crop_end=-1,
                 overwrite=True, use_frames=False, dose_metric="DWD",
                 dose_units="kGy", num_consec_frames=[3], frame_comp=1,
                 P_threshold=[0.01], plot_dir="Plots", dose_dir="Doses",
                 diode_dir="Diode_Readings", rp_comp_dir="RP_Comparisons",
                 plot_file_type="pdf", overwrite_doses=False,
                 rd_onset_dir="RD_Onset"):
        # Return dictionary of Compound objects
        self.compounds = process_compounds(cmpd_list=compound_list,
                                           crop_start=crop_start,
                                           crop_end=crop_end,
                                           overwrite=overwrite,
                                           buffer_subtraction=buffer_sub,
                                           average_type=average_type,
                                           use_frames=use_frames,
                                           dose_metric=dose_metric,
                                           dose_units=dose_units,
                                           plot_dir=plot_dir,
                                           dose_dir=dose_dir,
                                           diode_dir=diode_dir,
                                           overwrite_doses=overwrite_doses,
                                           rd_onset_dir=rd_onset_dir)

        # Create a dataframe with compound information
        self.cmpd_df = self.create_compound_df(num_consec_frames,
                                               P_threshold,
                                               runs_per_conc)
        # Write dataframe to csv
        self.cmpd_df.to_csv("compounds.csv")
        self.plot_RD_onset_ratios(plot_dir=plot_dir, rp_comp_dir=rp_comp_dir,
                                  plot_file_type=plot_file_type)

        # Create comparison plots
        # self.concentration_dependence_plots(plot_dir, plot_file_type)
        # self.radioprotectant_comparison_plot(plot_dir, rp_comp_dir,
        #                                      plot_file_type)

    # ----------------------------------------------------------------------- #
    #                         INSTANCE METHODS                                #
    # ----------------------------------------------------------------------- #

    def create_compound_df(self, num_consec_frames=[3], P_threshold=[0.01],
                           num_runs_per_conc=3):
        """Create a dataframe containing information about each compound's
        radioprotection ability.
        """
        num_df_rows = (len(Compound.CMPD_CONC) * num_runs_per_conc *
                       len(self.compounds))
        index = list(range(num_df_rows))
        columns = ['Dose (kGy)', 'RD Onset Frame', 'Compound',
                   'Concentration (mM)', 'Run Number', 'RD Metric',
                   'Num Consec Frames', 'P threshold']
        df = pd.DataFrame(columns=columns, index=index)
        counter = 0
        for cmpd in self.compounds.itervalues():
            warn_if_frames123_dont_correlate(cmpd, num_runs_per_conc)
            for n in num_consec_frames:
                for P in P_threshold:
                    raddam_onset_cormap = cmpd.get_raddam_onset_nums(n=n, P_thresh=P)
                    for conc in Compound.CMPD_CONC:
                        for run_num in xrange(0, num_runs_per_conc):
                            if cmpd.name == "no_protection":
                                frame_num_cormap = raddam_onset_cormap[0][run_num]
                                dose_value_cormap = cmpd.doses[0][frame_num_cormap-1, run_num]

                                frame_num_bsxcube = cmpd.raddam_onset[0][run_num]
                                dose_value_bsxcube = cmpd.doses[0][frame_num_bsxcube-1, run_num]
                            else:
                                frame_num_cormap = raddam_onset_cormap[conc][run_num]
                                dose_value_cormap = cmpd.doses[conc][frame_num_cormap-1, run_num]

                                frame_num_bsxcube = cmpd.raddam_onset[conc][run_num]
                                dose_value_bsxcube = cmpd.doses[conc][frame_num_bsxcube-1, run_num]
                            df.loc[counter] = pd.Series({columns[0]: dose_value_cormap,
                                                        columns[1]: frame_num_cormap,
                                                        columns[2]: cmpd.CMPD_INFO[cmpd.name][cmpd.LIST_INDEX["preferred_name"]],
                                                        columns[3]: conc,
                                                        columns[4]: run_num + 1,
                                                        columns[5]: "CorMap",
                                                        columns[6]: n,
                                                        columns[7]: P})
                            counter += 1
                            df.loc[counter] = pd.Series({columns[0]: dose_value_bsxcube,
                                                         columns[1]: frame_num_bsxcube,
                                                         columns[2]: cmpd.CMPD_INFO[cmpd.name][cmpd.LIST_INDEX["preferred_name"]],
                                                         columns[3]: conc,
                                                         columns[4]: run_num + 1,
                                                         columns[5]: "BsxCuBE",
                                                         columns[6]: n,
                                                         columns[7]: P})
                            counter += 1
        return df

    def plot_RD_onset_ratios(self, P_thresh=0.01, num_consec_frames=3,
                             plot_dir="Plots", rp_comp_dir="RP_Comparisons",
                             plot_file_type="pdf"):
        np_df = self.cmpd_df.loc[(self.cmpd_df['Compound'] == "No Protection") &
                                 (self.cmpd_df['RD Metric'] == "CorMap") &
                                 (self.cmpd_df['Num Consec Frames'] == num_consec_frames) &
                                 (self.cmpd_df['P threshold'] == P_thresh)]
        np_median = np_df["Dose (kGy)"].median()

        cmpd_ratios = {}
        for cmpd_data in self.compounds.itervalues():
            cmpd_name = cmpd_data.CMPD_INFO[cmpd_data.name][cmpd_data.LIST_INDEX["preferred_name"]]
            compound_df = self.cmpd_df.loc[(self.cmpd_df['Compound'] == cmpd_name) &
                                           (self.cmpd_df['RD Metric'] == "CorMap") &
                                           (self.cmpd_df['Num Consec Frames'] == num_consec_frames) &
                                           (self.cmpd_df['P threshold'] == P_thresh)]
            ratio_vals = np.zeros(len(cmpd_data.CMPD_CONC))
            for j, conc in enumerate(cmpd_data.CMPD_CONC):
                cmpd_median = compound_df[compound_df["Concentration (mM)"] == conc]["Dose (kGy)"].median()
                ratio_vals[j] = cmpd_median/np_median
            cmpd_ratios[cmpd_name] = ratio_vals
        columns = ['RD Onset Ratio', 'Compound', 'Concentration (mM)']
        index = list(range(len(cmpd_data.CMPD_CONC) * (len(cmpd_ratios) - 1)))
        df = pd.DataFrame(columns=columns, index=index)
        counter = 0
        for cmpd_name, ratio_vals in cmpd_ratios.iteritems():
            if cmpd_name not in "No Protection":
                for i, ratio in enumerate(ratio_vals):
                    df.loc[counter] = pd.Series({columns[0]: ratio,
                                                 columns[1]: cmpd_name,
                                                 columns[2]: Compound.CMPD_CONC[i]})
                    counter += 1
        fig, ax = plt.subplots(1, 1)
        ax = fig.gca()
        df.groupby("Compound").plot(x="Concentration (mM)", y="RD Onset Ratio",
                                    ax=ax, style=['-o'])
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
        plt.legend([v[0] for v in df.groupby("Compound")["Compound"]],
                   loc='center left', bbox_to_anchor=(1, 0.5))
        plt.ylabel("RD Onset Ratio")
        ax.invert_xaxis()
        plot_path = "../{}/{}/RatioPlots.{}".format(plot_dir, rp_comp_dir,
                                                    plot_file_type)
        plt.savefig(plot_path)
        plt.close('all')

    def concentration_dependence_plots(self, plot_dir, file_type):
        """Create box plots showing the dependence of compound efficacy with
        different concentrations.
        """
        for cmpd_data in self.compounds.itervalues():
            cmpd_name = cmpd_data.CMPD_INFO[cmpd_data.name][cmpd_data.LIST_INDEX["preferred_name"]]
            chart, ax = plt.subplots()
            sns_plot = sns.boxplot(x="Concentration (mM)", y="Dose (kGy)",
                                   data=self.cmpd_df[self.cmpd_df["Compound"] == cmpd_name],
                                   ax=ax, hue="RD Metric", palette="PRGn")
            sns_plot.set_title("{}".format(cmpd_name))
            make_data_dirs(plot_dir)
            file_loc = "../{}/{}.{}".format(plot_dir, cmpd_name, file_type)
            sns_plot.figure.savefig(file_loc)
            sns.plt.close(sns_plot.figure)

    def radioprotectant_comparison_plot(self, plot_dir, comp_dir, file_type):
        """Create box plots showing the efficacy of each radioprotectant
        compound
        """
        for conc in Compound.CMPD_CONC:
            chart, ax = plt.subplots()
            sns_plot = sns.boxplot(x="Compound", y="Dose (kGy)", ax=ax,
                                   data=self.cmpd_df[self.cmpd_df["Concentration (mM)"] == conc],
                                   hue="RD Metric", palette="PRGn")
            sns_plot.set_title("Concentration = {} (mM)".format(conc))
            make_data_dirs(plot_dir, comp_dir)
            file_loc = "../{}/{}/Concentration_{}.{}".format(plot_dir,
                                                             comp_dir,
                                                             conc, file_type)
            sns_plot.figure.savefig(file_loc)
            sns.plt.close(sns_plot.figure)


# ----------------------------------------------------------------------- #
#                               FUNCTIONS                                 #
# ----------------------------------------------------------------------- #


def process_compounds(cmpd_list, buffer_subtraction=True, average_type="mean",
                      crop_start=1, crop_end=-1, overwrite=True,
                      use_frames=False, dose_metric="DWD", dose_units="kGy",
                      plot_dir="Plots", dose_dir="Doses",
                      diode_dir="Diode_Readings", overwrite_doses=False,
                      rd_onset_dir="RD_Onset"):
    """Create compound objects from the list of compounds and store them in a
    dictionary
    """
    cmpd_data = {}
    for cmpd in cmpd_list:
        cmpd_data[cmpd] = Compound(compound_name=cmpd, crop_start=crop_start,
                                   crop_end=crop_end, overwrite=overwrite,
                                   buffer_subtraction=buffer_subtraction,
                                   average_type=average_type,
                                   use_frames=use_frames,
                                   dose_metric=dose_metric,
                                   dose_units=dose_units, plot_dir=plot_dir,
                                   dose_dir=dose_dir, diode_dir=diode_dir,
                                   overwrite_doses=overwrite_doses,
                                   rd_onset_dir=rd_onset_dir)
    return cmpd_data


def warn_if_frames123_dont_correlate(cmpd, num_runs_per_conc=3):
    """Function to print Warning if the first 3 frames of a run don't correlate
    with each other.
    """
    adjP_123_correlate_results = cmpd.check_frames_123_correlation()
    for j, conc in enumerate(cmpd.CMPD_CONC):
        if cmpd.name == "no_protection":
            if j == 0:
                conc = 0
            else:
                break
        for run in xrange(0, num_runs_per_conc):
            if not adjP_123_correlate_results[conc][run]:
                print '******************** WARNING *********************'
                print 'FIRST THREE FRAMES DO NOT CORRELATE'
                print """Check frames for compound: {}, concentration: {},
and run number: {}""".format(cmpd.CMPD_INFO[cmpd.name][cmpd.LIST_INDEX["preferred_name"]],
                             conc, run + 1)


def make_data_dirs(top_level_dir, second_level_dir=""):
    """Make directory structure for the comparison plots.
    """
    if not os.path.exists("../{}".format(top_level_dir)):
        os.makedirs("../{}".format(top_level_dir))
    if second_level_dir:
        if not os.path.exists("../{}/{}".format(top_level_dir,
                                                second_level_dir)):
            os.makedirs("../{}/{}".format(top_level_dir, second_level_dir))
