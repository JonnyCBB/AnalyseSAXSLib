"""Main script to analyse radiation damage in the SAXS data
"""
from CompareRadioProtectants import ComparisonAnalysis
# from RadioProtectant import Compound
# import pandas as pd
# import seaborn as sns
# import numpy as np
# from CorMapAnalysis import ScatterAnalysis
# from scattercurveutils import plot_1d_curve_with_filename
# FILE_NAME = "../20151214/GI_Gly/1d/GI_Gly_111_00001.dat"
# plot_1d_curve_with_filename(FILE_NAME, 1, 1300, True)

# compounds = ["np", "dtt", "suc", "tempo", "tre", "no3", "asc", "gly", "etgly"]
compounds = ["np", "asc"]
a = ComparisonAnalysis(compound_list=compounds, crop_start=25, crop_end=700,
                       overwrite=False)
# num_df_rows = len(Compound.CMPD_CONC) * num_runs_per_conc * len(compounds)
# index = list(range(num_df_rows))
# columns = ['Dose (kGy)', 'Frame Number', 'Compound', 'Concentration (mM)',
#            'Run Number']
# df = pd.DataFrame(columns=columns, index=index)
# counter = 0
# for cmpd in compounds:
#     cmpd_data = Compound(compound_name=cmpd, crop_start=25, crop_end=700,
#                          overwrite=False)
#     for conc in Compound.CMPD_CONC:
#         for run_num in xrange(0, num_runs_per_conc):
#             if cmpd_data.name == "no_protection":
#                 frame_num = cmpd_data.merge_thresholds[0][run_num]
#             else:
#                 frame_num = cmpd_data.merge_thresholds[conc][run_num]
#             dose_value = cmpd_data.doses[frame_num-1]
#             df.loc[counter] = pd.Series({columns[0]: dose_value,
#                                          columns[1]: frame_num,
#                                          columns[2]: cmpd,
#                                          columns[3]: conc,
#                                          columns[4]: run_num})
#             counter += 1

# files = "../cropped/GI_no_protection/1d/GI_np_030*.dat"
# files = "../20151214/GI_no_protection/1d/GI_np_034*.dat"
# files = "../cropped/GI_Asc/1d/GI_Asc_096*.dat"
# dose = np.linspace(1,120,120)/2.0
# metric = "dose"
# units = "kGy"
# a = ScatterAnalysis(files)
# b = ScatterAnalysis(files, dose, metric, units)
# a.plot_1d_intensity(1)
