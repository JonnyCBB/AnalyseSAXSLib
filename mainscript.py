"""Main script to analyse radiation damage in the SAXS data
"""
from CompareRadioProtectants import ComparisonAnalysis

# compounds = ["np", "dtt", "suc", "tempo", "tre", "no3", "asc", "gly", "etgly"]
compounds = ["np", "asc"]
a = ComparisonAnalysis(compound_list=compounds, crop_start=25, crop_end=700,
                       overwrite=False)
