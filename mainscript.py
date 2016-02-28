"""Main script to analyse radiation damage in the SAXS data
"""
# import numpy as np
# from CorMapAnalysis import ScatterAnalysis
from RadioProtectant import Compound
# from scattercurveutils import plot_1d_curve_with_filename
# FILE_NAME = "../20151214/GI_Gly/1d/GI_Gly_111_00001.dat"
# plot_1d_curve_with_filename(FILE_NAME, 1, 1300, True)

a = Compound(compound_name="np", crop_start=25, crop_end=700, overwrite=True)
b = Compound(compound_name="asc", crop_start=25, crop_end=700, overwrite=True)
# files = "../cropped/GI_no_protection/1d/GI_np_030*.dat"
# files = "../20151214/GI_no_protection/1d/GI_np_034*.dat"
# files = "../cropped/GI_Asc/1d/GI_Asc_096*.dat"
# dose = np.linspace(1,120,120)/2.0
# metric = "dose"
# units = "kGy"
# a = ScatterAnalysis(files)
# b = ScatterAnalysis(files, dose, metric, units)
# a.plot_1d_intensity(1)
