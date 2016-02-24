""" Class defining the radioprotectant compounds"""
import numpy as np
import math
import os
from RunSystemCommand import run_system_command


class Compound(object):
    """ Class containing fields and methods that can be applied to the
    radioprotectants used in the SAXS experiment.
    """
    # ----------------------------------------------------------------------- #
    #                         CLASS VARIABLES                                 #
    # ----------------------------------------------------------------------- #
    NUM_FRAMES = 120  # Number of frames taken in SAXS experiment
    NUM_POINTS_PER_FRAME = 1043  # Number of points in 1D intensity curve frame

    CMPD_CONC = [10, 5, 2, 1]  # In mM or v/v% if liquid.

    DATA_LOC_PREFIX = "../20151214"  # Location of original SAXS data

    SUB_DATA_LOC = "../subtracted"  # Location of subtracted SAXS curves

    CROPPED_DATA_LOC = "../cropped"  # Location of cropped data

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
                 average_type="mean", crop_start=1, crop_end=-1,
                 overwrite=True):

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #         MANIPULATE THE DATA - SUBTRACTION, CROPPING ETC.        #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        lowercase_cmpd_name = compound_name.lower()  # convert to lower case
        if lowercase_cmpd_name in self.CMPD_DICT:  # Check cmp is in dictionary
            self.name = self.CMPD_DICT[lowercase_cmpd_name]  # Extract name

            # Subtract the buffer if the user specifies
            if buffer_subtraction:
                average_buffer = self.get_average_runs(buffer_runs=True,
                                                       avg_type=average_type)
                self.subtract_average_buffer(average_buffer, overwrite)
                data_loc = self.SUB_DATA_LOC
            else:
                data_loc = self.DATA_LOC_PREFIX

            # Crop the data
            self.crop_data(data_loc, crop_start, crop_end, overwrite)

            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
            #                 CORRELATION ANALYSIS OF FRAMES                  #
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
            self.doses = self.get_doses()

        else:
            print '************************* ERROR **************************'
            print 'KEY: {} NOT FOUND'.format(compound_name)
            print """Please check the radioprotectant name or add to the
CMPD_DICT in the Compound class."""

# ----------------------------------------------------------------------- #
#                         INSTANCE METHODS                                #
# ----------------------------------------------------------------------- #

    def get_1d_curve_filename_prefix(self, run_number, new_data_loc_prfx=""):
        """Return the filename prefix for a SAXS run
        """
        # If a file prefix location is not given then use the location where
        # original data is stored. Otherwise use the user defined location.
        if not new_data_loc_prfx:
            data_loc_prfx = self.DATA_LOC_PREFIX
        else:
            data_loc_prfx = new_data_loc_prfx
        return "{}/{}_{}/{}/{}_{}_{run_num:03d}".format(data_loc_prfx,
                                                        self.PROTEIN_SAMPLE,
                                                        self.name,
                                                        "1d",
                                                        self.PROTEIN_SAMPLE,
                                                        self.CMPD_INFO[self.name][self.LIST_INDEX["dat_file_prefix"]],
                                                        run_num=run_number)

    def get_run_numbers(self, buffers=False, concentration=0.0):
        """Return the numbers of the filename prefixes that correspond to the
        runs for a given concentration
        """
        # Get number of the initial run in the relevant compound directory.
        initial_run_number = self.CMPD_INFO[self.name][self.LIST_INDEX["run_number_range"]][0]
        if buffers:
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
            if concentration in self.CMPD_CONC:
                conc_index = self.CMPD_CONC.index(concentration)
                run1 = (conc_index * 4) + initial_run_number + 1
                run2 = (conc_index * 4) + initial_run_number + 2
                run3 = (conc_index * 4) + initial_run_number + 3
                return [run1, run2, run3]
            else:
                print '*********************** ERROR ************************'
                print 'Concentration not given!'
                print """Rerun the get_run_numbers method with a concentration
that was used in the experiment: 1, 2, 5 or 10 mM."""

    def get_average_runs(self, cmpd_concentration=0.0, buffer_runs=False,
                         avg_type="mean"):
        """Get the average value of the intensity of frames for a particular
        concentration (buffer or protein sample) of the radioprotectant.
        """
        # Get the run numbers corresponding to the required run
        run_nums = self.get_run_numbers(buffers=buffer_runs,
                                        concentration=cmpd_concentration)

        # Preallocate the matrix that will contain the averaged intensities
        # and variances of the runs
        averaged_runs = np.zeros([len(run_nums),
                                  self.NUM_POINTS_PER_FRAME])
        averaged_sig_runs = np.zeros([len(run_nums),
                                      self.NUM_POINTS_PER_FRAME])

        # Go through each of the runs corresponding to the buffer/concentration
        # and average all of the frames that were collected during that run
        # (both the average intensity and the corresponding averaged standard
        # deviation). Then store the result of the averaged frames in the
        # relevant row of the preallocated arrays and return those arrays at
        # the end.
        for i, run_num in enumerate(run_nums):
            frames = np.zeros([self.NUM_FRAMES,
                               self.NUM_POINTS_PER_FRAME])
            sig_frames = np.zeros([self.NUM_FRAMES,
                                   self.NUM_POINTS_PER_FRAME])
            dat_file_prefix = self.get_1d_curve_filename_prefix(run_num)
            for frame in range(self.NUM_FRAMES):
                frame_data = np.loadtxt(get_1d_curve_filename(dat_file_prefix,
                                                              frame+1))
                frames[frame, :] = frame_data[:, 1]
                sig_frames[frame, :] = frame_data[:, 2]
            if avg_type == "median":
                averaged_runs[i, :] = np.median(frames, axis=0)
            elif avg_type == "mean":
                averaged_runs[i, :] = np.mean(frames, axis=0)
            averaged_sig_runs[i, :] = ((1.0 / math.pow(self.NUM_FRAMES, 2)) *
                                       np.sum(np.square(sig_frames), axis=0))
        return averaged_runs, averaged_sig_runs

    def subtract_average_buffer(self, avg_buffer, overwrite_data):
        """Substract the average buffer from all of the frames for all
        concentrations. The result is a new folder containing all of the
        subtracted data.
        """
        top_level_already_existed = self.make_data_dirs(self.SUB_DATA_LOC)

        if not top_level_already_existed or overwrite_data:
            # For each frame for each of the different concentrations, subtract
            # the averaged buffer values from the intensity curves and write
            # them to file in a directory structure that mimics that of the
            # structure where the original data are stored.
            for i, conc in enumerate(self.CMPD_CONC):
                run_nums = self.get_run_numbers(concentration=conc,
                                                buffers=False)
                for run_num in run_nums:
                    dat_file_prfx = self.get_1d_curve_filename_prefix(run_num)
                    for frame in range(self.NUM_FRAMES):
                        dat_file = get_1d_curve_filename(dat_file_prfx,
                                                         frame+1)
                        avg_buff = np.column_stack((avg_buffer[0][i, :],
                                                    avg_buffer[1][i, :]))
                        subtracted_frame = subtract_frame_from_file(dat_file,
                                                                    avg_buff)

                        # Now the frames have been subtracted, write the file.
                        sub_filename = get_1d_curve_filename(self.get_1d_curve_filename_prefix(run_num, self.SUB_DATA_LOC), frame+1)
                        header = """Sample description: {}
Parent(s): {}
Compound concentration: {} mM""".format(self.PROTEIN_SAMPLE, dat_file, conc)
                        np.savetxt(fname=sub_filename, X=subtracted_frame,
                                   header=header, delimiter="\t", fmt="%.6e")
        else:
            print """'subtraction' directory already existed and it was
specified not to overwrite the data. If you want to create the data then set
the overwrite to 'True'."""

    def crop_data(self, data_location, crp_start, crp_end, overwrite_data):
        """Crop the 1d scattering curves
        """
        # If the crp_end value is below or equal to zero then set it to the max
        # value of the number of points per frame.
        if crp_end <= 0:
            crp_end = self.NUM_POINTS_PER_FRAME

        # Make directory tree for the cropped data.
        top_level_already_existed = self.make_data_dirs(self.CROPPED_DATA_LOC)

        if not top_level_already_existed or overwrite_data:
            # For each frame for each of the different concentrations, crop the
            # frame to the specified level
            for i, conc in enumerate(self.CMPD_CONC):
                run_nums = self.get_run_numbers(concentration=conc,
                                                buffers=False)
                for run_num in run_nums:
                    for frame in range(self.NUM_FRAMES):
                        # Get location of subtracted and cropped files
                        sub_filename = get_1d_curve_filename(self.get_1d_curve_filename_prefix(run_num, self.SUB_DATA_LOC), frame+1)
                        crp_filename = get_1d_curve_filename(self.get_1d_curve_filename_prefix(run_num, self.CROPPED_DATA_LOC),
                                                             frame+1)

                        # create the command string to be run
                        cmd = "datcrop {} --first {:d} --last {:d} -o {}".format(sub_filename, crp_start, crp_end, crp_filename)
                        run_system_command(command_string=cmd)

    def make_data_dirs(self, top_level_dir):
        """Make directory structure identical to the one created for the
        original data. This is so the methods written to find the correct 1d
        scattering curve data can be used without modification on the newly
        created data.
        """
        top_level_dir_already_exists = False
        if not os.path.exists(top_level_dir):
            os.makedirs(top_level_dir)
        else:
            top_level_dir_already_exists = True
        if not os.path.exists("{}/{}_{}".format(top_level_dir,
                                                self.PROTEIN_SAMPLE,
                                                self.name)):
            os.makedirs("{}/{}_{}".format(top_level_dir, self.PROTEIN_SAMPLE,
                                          self.name))
        if not os.path.exists("{}/{}_{}/{}".format(top_level_dir,
                                                   self.PROTEIN_SAMPLE,
                                                   self.name, "1d")):
            os.makedirs("{}/{}_{}/{}".format(top_level_dir,
                                             self.PROTEIN_SAMPLE, self.name,
                                             "1d"))
        return top_level_dir_already_exists

    def get_doses(self):
        """Run RADDOSE-3D to get doses for the Radioprotectant compound
        """
        ###################################################################
        # NEED TO SORT THIS METHOD OUT
        ###################################################################
        pass

# ----------------------------------------------------------------------- #
#                               FUNCTIONS                                 #
# ----------------------------------------------------------------------- #


def get_1d_curve_filename(curve_prefix, frame_num):
    """Return the filename for a 1d scattering curve
    """
    return '{prfx}_{fr:05d}.dat'.format(prfx=curve_prefix, fr=frame_num)


def subtract_frame_from_file(frame1_file, frame2):
    """Subtract frame2 to from frame1_file where frame1_file is the
    filename where the data for frame1 is stored.
    """
    frame_data = np.loadtxt(frame1_file)
    return subtract_frames(frame_data, frame2)


def subtract_frames(frame1, frame2):
    """Method used to subtract data extracted from curves of two 1D
    scattering frames. frame1 - frame2
    """
    if frame2.shape[1] == 3:
        frame2 = frame2[:, 1:]
    sub_intensities = frame1[:, 1] - frame2[:, 0]
    updated_variances = np.square(frame1[:, 2]) + np.square(frame2[:, 1])
    updated_sigmas = np.sqrt(updated_variances)
    return np.column_stack((frame1[:, 0], sub_intensities, updated_sigmas))
