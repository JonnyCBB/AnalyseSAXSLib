import ipdb
""" Class dealing with all of the analysis of the multiple compounds used in
the SAXS analysis experiments
"""
from RadioProtectant import Compound
import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt


class ComparisonAnalysis(object):
    """Analysis of several compounds for SAXS analysis.
    """
    # ----------------------------------------------------------------------- #
    #                         CONSTRUCTOR METHOD                              #
    # ----------------------------------------------------------------------- #
    def __init__(self, compound_list, runs_per_conc=3, buffer_sub=True,
                 average_type="mean", crop_start=1, crop_end=-1,
                 overwrite=True, use_frames=False, dose_metric="DWD",
                 dose_units="kGy", num_consec_frames=3, frame_comp=1,
                 P_threshold=0.01, plot_dir="Plots", dose_dir="Doses",
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
                                           num_consec_frames=num_consec_frames,
                                           frame_comp=frame_comp,
                                           P_threshold=P_threshold,
                                           plot_dir=plot_dir,
                                           dose_dir=dose_dir,
                                           diode_dir=diode_dir,
                                           overwrite_doses=overwrite_doses,
                                           rd_onset_dir=rd_onset_dir)
        # Create a dataframe with compound information
        self.cmpd_df = self.create_compound_df(runs_per_conc)

        # Create comparison plots
        # self.concentration_dependence_plots(plot_dir, plot_file_type)
        # self.radioprotectant_comparison_plot(plot_dir, rp_comp_dir,
        #                                      plot_file_type)

    # ----------------------------------------------------------------------- #
    #                         INSTANCE METHODS                                #
    # ----------------------------------------------------------------------- #

    def create_compound_df(self, num_runs_per_conc=3):
        """Create a dataframe containing information about each compound's
        radioprotection ability.
        """
        num_df_rows = (len(Compound.CMPD_CONC) * num_runs_per_conc *
                       len(self.compounds))
        index = list(range(num_df_rows))
        columns = ['Dose (kGy)', 'RD Onset Frame', 'Compound',
                   'Concentration (mM)', 'Run Number', 'RD Metric']
        df = pd.DataFrame(columns=columns, index=index)
        counter = 0
        for cmpd_data in self.compounds.itervalues():
            for conc in Compound.CMPD_CONC:
                for run_num in xrange(0, num_runs_per_conc):
                    if cmpd_data.name == "no_protection":
                        frame_num_cormap = cmpd_data.raddam_onset_cormap[0][run_num]
                        dose_value_cormap = cmpd_data.doses[0][frame_num_cormap-1, run_num]

                        frame_num_bsxcube = cmpd_data.raddam_onset[0][run_num]
                        dose_value_bsxcube = cmpd_data.doses[0][frame_num_bsxcube-1, run_num]
                    else:
                        frame_num_cormap = cmpd_data.raddam_onset_cormap[conc][run_num]
                        dose_value_cormap = cmpd_data.doses[conc][frame_num_cormap-1, run_num]

                        frame_num_bsxcube = cmpd_data.raddam_onset[conc][run_num]
                        dose_value_bsxcube = cmpd_data.doses[conc][frame_num_bsxcube-1, run_num]
                    df.loc[counter] = pd.Series({columns[0]: dose_value_cormap,
                                                columns[1]: frame_num_cormap,
                                                columns[2]: cmpd_data.CMPD_INFO[cmpd_data.name][cmpd_data.LIST_INDEX["preferred_name"]],
                                                columns[3]: conc,
                                                columns[4]: run_num + 1,
                                                columns[5]: "CorMap"})
                    counter += 1
                    df.loc[counter] = pd.Series({columns[0]: dose_value_bsxcube,
                                                 columns[1]: frame_num_bsxcube,
                                                 columns[2]: cmpd_data.CMPD_INFO[cmpd_data.name][cmpd_data.LIST_INDEX["preferred_name"]],
                                                 columns[3]: conc,
                                                 columns[4]: run_num + 1,
                                                 columns[5]: "BsxCuBE"})
                    counter += 1
        return df

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
                      num_consec_frames=3, frame_comp=1, P_threshold=0.01,
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
                                   dose_units=dose_units,
                                   num_consec_frames=num_consec_frames,
                                   frame_comp=frame_comp,
                                   P_threshold=P_threshold, plot_dir=plot_dir,
                                   dose_dir=dose_dir, diode_dir=diode_dir,
                                   overwrite_doses=overwrite_doses,
                                   rd_onset_dir=rd_onset_dir)
    return cmpd_data


def make_data_dirs(top_level_dir, second_level_dir=""):
    """Make directory structure for the comparison plots.
    """
    if not os.path.exists("../{}".format(top_level_dir)):
        os.makedirs("../{}".format(top_level_dir))
    if second_level_dir:
        if not os.path.exists("../{}/{}".format(top_level_dir,
                                                second_level_dir)):
            os.makedirs("../{}/{}".format(top_level_dir, second_level_dir))
