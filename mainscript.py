"""Main script to analyse radiation damage in the SAXS data
"""
from CompareRadioProtectants import ComparisonAnalysis

# compounds = ["np", "dtt", "suc", "tempo", "tre", "no3", "asc", "gly", "etgly"]
compounds = ["np", "asc"]
a = ComparisonAnalysis(compound_list=compounds, crop_start=25, crop_end=700,
                       overwrite=False)

for cmpd in a.compounds.itervalues():
    for j, conc in enumerate(cmpd.CMPD_CONC):
        if not cmpd.adjP_123_correlate:
            print '******************** WARNING *********************'
            print 'FIRST THREE FRAMES DO NOT CORRELATE'
            print """Check frames for compound: {}, concentration: {},
        and run number: {}""".format(cmpd.CMPD_INFO[cmpd.name][cmpd.LIST_INDEX["preferred_name"]],
                                     conc, j + 1)
