"""Main script to analyse radiation damage in the SAXS data
"""
from CompareRadioProtectants import ComparisonAnalysis

P_thresholds = [0.01, 0.05, 0.10]
num_consec_frames = [1, 3, 5, 7, 10]
# compounds = ["np", "asc", "dtt", "suc", "tempo", "tre", "no3", "gly", "etgly"]
compounds = ["np", "asc"]
a = ComparisonAnalysis(compound_list=compounds, crop_start=25, crop_end=700,
                       overwrite=False, overwrite_doses=False,
                       num_consec_frames=num_consec_frames,
                       P_threshold=P_thresholds)
