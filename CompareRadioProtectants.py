""" Class dealing with all of the analysis of the multiple compounds used in
the SAXS analysis experiments
"""
from RadioProtectant import Compound
import pandas as pd


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
                 P_threshold=0.01):
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
                                           P_threshold=P_threshold)
        # Create a dataframe with compound information
        self.cmpd_df = self.create_compound_df(runs_per_conc)

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
        columns = ['Dose (kGy)', 'Frame Number', 'Compound',
                   'Concentration (mM)', 'Run Number']
        df = pd.DataFrame(columns=columns, index=index)
        counter = 0
        for cmpd_data in self.compounds.itervalues():
            for conc in Compound.CMPD_CONC:
                for run_num in xrange(0, num_runs_per_conc):
                    if cmpd_data.name == "no_protection":
                        frame_num = cmpd_data.merge_thresholds[0][run_num]
                    else:
                        frame_num = cmpd_data.merge_thresholds[conc][run_num]
                    dose_value = cmpd_data.doses[frame_num-1]
                    df.loc[counter] = pd.Series({columns[0]: dose_value,
                                                 columns[1]: frame_num,
                                                 columns[2]: cmpd_data.CMPD_INFO[cmpd_data.name][cmpd_data.LIST_INDEX["preferred_name"]],
                                                 columns[3]: conc,
                                                 columns[4]: run_num})
                    counter += 1
        return df

# ----------------------------------------------------------------------- #
#                               FUNCTIONS                                 #
# ----------------------------------------------------------------------- #


def process_compounds(cmpd_list, buffer_subtraction=True, average_type="mean",
                      crop_start=1, crop_end=-1, overwrite=True,
                      use_frames=False, dose_metric="DWD", dose_units="kGy",
                      num_consec_frames=3, frame_comp=1, P_threshold=0.01):
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
                                   P_threshold=P_threshold)
    return cmpd_data
