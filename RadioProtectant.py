""" Class defining the radioprotectant compounds"""
import numpy as np
import math

class Compound(object):
    """ Class containing fields and methods that can be applied to the
    radioprotectants used in the SAXS experiment.
    """
    # ----------------------------------------------------------------------- #
    #                         CLASS VARIABLES                                 #
    # ----------------------------------------------------------------------- #
    NUM_FRAMES = 120  # Number of frames taken in SAXS experiment
    NUM_POINTS_PER_FRAME = 1043  # Number of points in 1D intensity curve frame

    COMPOUND_CONC = [10, 5, 2, 1]  # In mM or v/v% if liquid.

    DATA_LOC_PREFIX = "../20151214"  # Location of SAXS data

    PROTEIN_SAMPLE = "GI"  # File prefix used for data (Glucose Isomerase)

    # Dictionary allowing user to define various names for the radioprotectant
    # compounds. This is in case they don't know the symbol used for the
    # compound in the data files.
    CMPD_DICT = {"none": "no_protection",
                 "n/a": "no_protection",
                 "ascorbate": "Asc",
                 "sodium ascorbate": "Asc",
                 "naac": "Asc",
                 "asc": "Asc",
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

    # Dictionary containing various information about each radioprotectant.
    # This includes the file prefixes used to find various files for each
    # radio protectant and the preferred name for each compound.
    CMPD_INFO = {"no_protection": ["np", [29, 35], [14], "None"],
                 "Asc": ["Asc", [95, 110], [10, 9], "Ascorbate"],
                 "DTT": ["DTT", [191, 206], [5, 4], "DTT"],
                 "EtGly": ["EtGly", [79, 94], [11, 10], "Ethylene Glycol"],
                 "Gly": ["Gly", [111, 126], [9, 8], "Glycerol"],
                 "NO3_2": ["NO3", [159, 174], [7, 6], "Sodium Nitrate"],
                 "Sucrose": ["Sucrose", [143, 158], [8, 7], "Sucrose"],
                 "TEMPO": ["TEMPO", [127, 142], [8], "TEMPO"],
                 "trehalose_2": ["th", [175, 190], [6, 5], "Trehalose"]}

    # Dictionary that gives the indices of the relevant information in the list
    # which is stored as the "value" in the CMPD_INFO dictionary variable
    # defined above.
    LIST_INDEX = {"dat_file_prefix": 0,
                  "run_number_range": 1,
                  "bsxcube_log": 2,
                  "preferred_name": 3}


# ----------------------------------------------------------------------- #
#                         CONSTRUCTOR METHOD                              #
# ----------------------------------------------------------------------- #
    def __init__(self, compound_name, buffer_subtraction=True,
                 average_type="median"):
        lowercase_cmpd_name = compound_name.lower()  # convert to lower case
        if lowercase_cmpd_name in self.CMPD_DICT:  # Check cmp is in dictionary
            self.name = self.CMPD_DICT[lowercase_cmpd_name]  # Extract name
            if buffer_subtraction:  # If user wants to subtract the buffer
                self.average_buffer = self.get_average_buffer(average_type)
            self.doses = self.get_doses()

        else:
            print '************************* ERROR **************************'
            print 'KEY: {} NOT FOUND'.format(compound_name)
            print """Please check the radioprotectant name or add to the
CMPD_DICT in the Compound class."""

# ----------------------------------------------------------------------- #
#                         INSTANCE METHODS                                #
# ----------------------------------------------------------------------- #

    def get_1d_curve_filename_prefix(self, run_number):
        """Return the filename prefix for a SAXS run
        """
        return "{}/{}_{}/{}/{}_{}_{run_num:03d}".format(self.DATA_LOC_PREFIX,
                                                        self.PROTEIN_SAMPLE,
                                                        self.name,
                                                        "1d",
                                                        self.PROTEIN_SAMPLE,
                                                        self.CMPD_INFO[self.name][self.LIST_INDEX["dat_file_prefix"]],
                                                        run_num=run_number)

    def get_run_numbers(self, buffer_runs, concentration=0.0):
        """Return the numbers of the filename prefixes that correspond to the
        runs for a given concentration
        """
        # Get number of the initial run in the relevant compound directory.
        initial_run_number = self.CMPD_INFO[self.name][self.LIST_INDEX["run_number_range"]][0]
        if buffer_runs:
            # If user wants the runs corresponding the the buffer then return
            # those.
            buffer1 = initial_run_number
            buffer2 = initial_run_number + 4
            buffer3 = initial_run_number + 8
            buffer4 = initial_run_number + 12
            return [buffer1, buffer2, buffer3, buffer4]
        else:
            # If user wants the runs corresponding to a specific concentration
            # of the compound then return these runs.
            if concentration in self.COMPOUND_CONC:
                conc_index = self.COMPOUND_CONC.index(concentration)
                run1 = (conc_index * 4) + initial_run_number + 1
                run2 = (conc_index * 4) + initial_run_number + 2
                run3 = (conc_index * 4) + initial_run_number + 3
                return [run1, run2, run3]
            else:
                print '*********************** ERROR ************************'
                print 'Concentration not given!'
                print """Rerun the get_run_numbers method with a Concentration
concentration that was used in the experiment: 1, 2, 5 or 10 mM."""

    def get_average_buffer(self, avg_type="mean"):
        """Get the average value of the buffer frames for a particular
        concentration of the radioprotectant.
        """
        # Get the run numbers corresponding to the buffer
        buffer_nums = self.get_run_numbers(buffer_runs=True)

        # Preallocate the matrix that will contain the averaged intensities
        # and variances buffer runs
        averaged_buffer_runs = np.zeros([len(buffer_nums),
                                        self.NUM_POINTS_PER_FRAME])
        averaged_var_buffer_runs = np.zeros([len(buffer_nums),
                                            self.NUM_POINTS_PER_FRAME])

        # Go through each of the runs corresponding to the buffer and average
        # all of the frames that were collected during that run (both the
        # average intensity and the corresponding averaged standard deviation).
        # Then store the result of the averaged frames in the relevant row
        # of the preallocated arrays and return those arrays at the end.
        for i, buffer_run in enumerate(buffer_nums):
            dat_file_prefix = self.get_1d_curve_filename_prefix(buffer_run)
            buffer_frames = np.zeros([self.NUM_FRAMES,
                                      self.NUM_POINTS_PER_FRAME])
            buffer_var_frames = np.zeros([self.NUM_FRAMES,
                                         self.NUM_POINTS_PER_FRAME])
            for frame in range(self.NUM_FRAMES):
                dat_file = '{prefix}_{frame_num:05d}.dat'.format(prefix=dat_file_prefix, frame_num=frame+1)
                frame_data = np.loadtxt(dat_file)
                buffer_frames[frame, :] = frame_data[:, 1]
                buffer_var_frames[frame, :] = frame_data[:, 2]
            if avg_type == "median":
                averaged_buffer_runs[i, :] = np.median(buffer_frames, axis=0)
            elif avg_type == "mean":
                averaged_buffer_runs[i, :] = np.mean(buffer_frames, axis=0)
            sum_var = np.sum(np.square(buffer_var_frames), axis=0)
            scale = 1.0 / math.pow(self.NUM_FRAMES, 2)
            averaged_var_buffer_runs[i, :] = scale * sum_var
        return averaged_buffer_runs, averaged_var_buffer_runs

    def get_doses(self):
        """Run RADDOSE-3D to get doses for the Radioprotectant compound
        """
        ###################################################################
        # NEED TO SORT THIS METHOD OUT
        ###################################################################
        return np.ones(self.NUM_FRAMES)
