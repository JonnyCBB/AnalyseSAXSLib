""" Radio protectant class
"""

import sys


class Compound(object):
    """ Class containing fields and methods that can be applied to the
    radioprotectants used in the SAXS experiment.
    """
    # ----------------------------------------------------------------------- #
    #                         CLASS VARIABLES                                 #
    # ----------------------------------------------------------------------- #
    NUM_FRAMES = 120  # Number of frames taken in SAXS experiment

    RADIO_PROTECTANT_DICT = {"none": "no_protection",
                             "n/a": "no_protection",
                             "ascorbate": "Asc",
                             "sodium ascorbate": "Asc",
                             "naac": "Asc",
                             "dtt": "DTT",
                             "ethylene glycol": "EtGly",
                             "eg": "EtGly",
                             "etgly": "EtGly",
                             "ethgly": "EtGly",
                             "glycerol": "Gly",
                             "gly": "Gly",
                             "sodium nitrate": "NO3_2",
                             "nitrate": "NO3_2",
                             "no3": "NO3_2",
                             "nano3": "NO3_2",
                             "sucrose": "Sucrose",
                             "suc": "Sucrose",
                             "tempo": "TEMPO",
                             "trehalose": "trehalose_2",
                             "tre": "trehalose_2"}

    RP_INFO = {"no_protection": ["np", [29, 35], [14], "None"],
               "Asc": ["Asc", [95, 110], [10, 9], "Ascorbate"],
               "DTT": ["DTT", [191, 206], [5, 4], "DTT"],
               "EtGly": ["EtGly", [79, 94], [11, 10], "Ethylene Glycol"],
               "Gly": ["Gly", [111, 126], [9, 8], "Glycerol"],
               "NO3_2": ["NO3", [159, 174], [7, 6], "Sodium Nitrate"],
               "Sucrose": ["Sucrose", [143, 158], [8, 7], "Sucrose"],
               "TEMPO": ["TEMPO", [127, 142], [8], "TEMPO"],
               "trehalose_2": ["th", [175, 190], [6, 5], "Trehalose"]}

    LIST_INDEX = {"dat_file_prefix": 0,
                  "run_number_range": 1,
                  "bsxcube_log": 2,
                  "preferred_name": 3}


# ----------------------------------------------------------------------- #
#                         CONSTRUCTOR METHOD                              #
# ----------------------------------------------------------------------- #
    def __init__(self, radio_protectant_name):
        lowercase_rad_prot_name = radio_protectant_name.lower()
        if lowercase_rad_prot_name in self.RADIO_PROTECTANT_DICT:
            self.name = self.RADIO_PROTECTANT_DICT[lowercase_rad_prot_name]
        else:
            print '************************* ERROR **************************'
            print 'KEY: {} NOT FOUND'.format(radio_protectant_name)
            print """Please check the radioprotectant name or add to the
RADIO_PROTECTANT_DICT in the Compound class."""
            sys.exit()
