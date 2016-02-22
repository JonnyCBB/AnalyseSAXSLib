"""Main script to analyse radiation damage in the SAXS data
"""
from RadioProtectant import Compound
from scattercurveutils import plot_1d_curve_with_filename
FILE_NAME = "../20151214/GI_Gly/1d/GI_Gly_111_00001.dat"
plot_1d_curve_with_filename(FILE_NAME, 1, 1200, True)
