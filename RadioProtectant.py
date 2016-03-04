""" Class that looks after each Radio protectant compound. This includes
sorting information about the files that store the information for each
compound, the dose accumulated for each compound and their efficacy for
protecting against radiation damage.
Since the file locations are kept here, this class also deals with manipulation
of the data in the files, i.e. buffer avergaing and subtraction and cropping.
"""
import numpy as np
import math
import os
from CorMapAnalysis import ScatterAnalysis
import matplotlib.pyplot as plt
from Raddose3D import Raddose3d


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

    # Scale factor to convert diodes readings to flux
    FLUX_SCALE_FAC = 2722950000000000

    FLUX_ADD_FAC = 5.72293

    PLOT_NUM = 0  # Keep count of number of plots

    EXP_PER_FRAME = 1  # Exposure time per frame

    MGY_TO_KGY = 1000  # Factor to convert MGy to kGy

    # Dictionary allowing user to define various names for the radioprotectant
    # compounds. This is in case they don't know the symbol used for the
    # compound in the data files.
    CMPD_DICT = {"none": "no_protection",
                 "n/a": "no_protection",
                 "np": "no_protection",
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
    CMPD_INFO = {"no_protection": ["np", [29, 35], [14], "No Protection"],
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
                 overwrite=True, use_frames=False, dose_metric="DWD",
                 dose_units="kGy", num_consec_frames=3, frame_comp=1,
                 P_threshold=0.01, plot_dir="Plots", dose_dir="Doses",
                 diode_dir="Diode_Readings", rd_onset_dir="RD_Onset",
                 overwrite_doses=False):

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #         MANIPULATE THE DATA - SUBTRACTION, CROPPING ETC.        #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        lowercase_cmpd_name = compound_name.lower()  # convert to lower case
        if lowercase_cmpd_name in self.CMPD_DICT:  # Check cmp is in dictionary
            self.name = self.CMPD_DICT[lowercase_cmpd_name]  # Extract name
            self.plot_dir = plot_dir

            # Subtract the buffer if the user specifies
            if buffer_subtraction:
                average_buffer = self.get_average_runs(buffer_runs=True,
                                                       avg_type=average_type)
                self.subtract_average_buffer(average_buffer, overwrite)
                data_loc = self.SUB_DATA_LOC
            else:
                data_loc = self.DATA_LOC_PREFIX

            # Crop the data
            self.crop_data(data_loc, crop_start, crop_end, overwrite,
                           buffer_subtraction)

            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
            #             DOSE VALUE CALCULATION WITH RADDOSE-3D              #
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
            # Only go through the long process of calculating the doses and
            # extracting the diode readings if the user actually specifies that
            # they want to overwrite the data or no data does not already exist
            if (overwrite_doses or not
                os.path.exists("../{}/{}".format(dose_dir, self.name)) or not
                    os.path.exists("../{}/{}".format(diode_dir, self.name)) or not
                    os.path.exists("../{}/{}".format(rd_onset_dir, self.name))):
                self.diode_readings, self.doses, self.raddam_onset = self.get_doses(dose_met=dose_metric)
                self.save_doses_and_diode_to_csv(dose_dir, diode_dir)
            else:
                self.diode_readings, self.doses, self.raddam_onset = self.get_saved_doses_and_diode(dose_dir, diode_dir, rd_onset_dir)

            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
            #                 CORRELATION ANALYSIS OF FRAMES                  #
            # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
            # Create ScatterAnalysis object for each run
            self.scat_analysis = self.get_data_analysis_objs(dose_metric,
                                                             dose_units)

            self.merge_thresholds, self.adjP_123_correlate = self.get_merging_thresholds(num_consec_frames, frame_comp, P_threshold)

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
            if self.name == "no_protection":
                return [29, 31, 33, 35]
            else:
                buffer1 = initial_run_number
                buffer2 = initial_run_number + 4
                buffer3 = initial_run_number + 8
                buffer4 = initial_run_number + 12
                return [buffer1, buffer2, buffer3, buffer4]
        else:
            # If user wants the runs corresponding to a specific concentration
            # of the compound then return these runs.
            if self.name == "no_protection":
                return [30, 32, 34]
            elif concentration in self.CMPD_CONC:
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
            averaged_sig_runs[i, :] = np.sqrt((1.0 / math.pow(self.NUM_FRAMES, 2)) *
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
                # Only need to loop through once for the "no_protection" sample
                # Because there aren't several concentrations of radio
                # protectant compounds in the sample.
                if i == 1 and self.name == "no_protection":
                    break
                else:
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
                            if self.name == "no_protection":
                                header = """Sample description: {}
Parent(s): {}
No protection""".format(self.PROTEIN_SAMPLE, dat_file)
                            else:
                                header = """Sample description: {}
Parent(s): {}
Compound concentration: {} mM""".format(self.PROTEIN_SAMPLE, dat_file, conc)
                            np.savetxt(fname=sub_filename, X=subtracted_frame,
                                       header=header, delimiter="\t",
                                       fmt="%.6e")
        else:
            print """'subtraction' directory already existed and it was
specified not to overwrite the data. If you want to create the data then set
the overwrite to 'True'."""

    def crop_data(self, data_location, crp_start, crp_end, overwrite_data,
                  buff_sub):
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
                # Only need to loop through once for the "no_protection" sample
                # Because there aren't several concentrations of radio
                # protectant compounds in the sample.
                if i == 1 and self.name == "no_protection":
                    break
                else:
                    for run_num in run_nums:
                        for frame in range(self.NUM_FRAMES):
                            # Get location of subtracted and cropped files
                            if buff_sub:
                                parent_filename = get_1d_curve_filename(self.get_1d_curve_filename_prefix(run_num, self.SUB_DATA_LOC), frame+1)
                            else:
                                parent_filename = get_1d_curve_filename(self.get_1d_curve_filename_prefix(run_num, self.DATA_LOC_PREFIX), frame+1)

                            crp_filename = get_1d_curve_filename(self.get_1d_curve_filename_prefix(run_num, self.CROPPED_DATA_LOC), frame+1)

                            # Crop data
                            frame_data = np.loadtxt(parent_filename)
                            cropped_frame_data = frame_data[crp_start-1:crp_end, :]
                            if self.name == "no_protection":
                                header = """Sample description: {}
Parent(s): {}
No protection""".format(self.PROTEIN_SAMPLE, parent_filename)
                            else:
                                header = """Sample description: {}
Parent(s): {}
Compound concentration: {} mM""".format(self.PROTEIN_SAMPLE,
                                        parent_filename, conc)
                            np.savetxt(fname=crp_filename,
                                       X=cropped_frame_data,
                                       header=header, delimiter="\t",
                                       fmt="%.6e")

                        # # create the command string to be run
                        # cmd = "datcrop {} --first {:d} --last {:d} -o {}".format(sub_filename, crp_start, crp_end, crp_filename)
                        # run_system_command(command_string=cmd)

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

    def get_data_analysis_objs(self, metric, units):
        """Return list of scatter analyis objects. There are three for each
        concentration corresponding to the three repeated runs for each
        concentration. The only exception is for the "no_protection" instance
        where there is only 1 lot of 3 runs.
        """
        scatter_objs = {}
        for i, conc in enumerate(self.CMPD_CONC):
            run_nums = self.get_run_numbers(concentration=conc,
                                            buffers=False)
            if i == 1 and self.name == "no_protection":
                break
            else:
                run_objs = []
                for j, run_num in enumerate(run_nums):
                    dat_file_prfx = self.get_1d_curve_filename_prefix(run_number=run_num,
                                                                      new_data_loc_prfx=self.CROPPED_DATA_LOC)
                    file_prefixes = "{}*.dat".format(dat_file_prfx)

                    if self.name == "no_protection":
                        run_objs.append(ScatterAnalysis(file_prefixes,
                                                        self.doses[0],
                                                        metric, units))
                    else:
                        run_objs.append(ScatterAnalysis(file_prefixes,
                                                        self.doses[conc][:, j],
                                                        metric, units))
                if self.name == "no_protection":
                    scatter_objs[0] = run_objs
                else:
                    scatter_objs[conc] = run_objs
        return scatter_objs

    def get_merging_thresholds(self, n=3, k=1, P_thresh=0.01):
        """Get the images corresponding to the first n consecutive frames that
        are different from the given frame. These are going to be referred to
        as the merging thresholds. Store the corresponding image for each run
        in a dictionary. The key will correspond to the concentration of the
        radioprotectant compound and the value will be a list containing the
        corresponding frames numbers of the merging threshold for each run.
        """
        merging_thresholds = {}
        adjP_123_correlate = {}
        for conc, conc_analysis in self.scat_analysis.iteritems():
            frames = list(range(3))
            for j, run in enumerate(conc_analysis):
                # Check that first three frames correlated with each other.
                adjP_1_2 = run.get_pw_data(1, 2, "adj P(>C)")
                adjP_1_3 = run.get_pw_data(1, 3, "adj P(>C)")
                adjP_2_3 = run.get_pw_data(2, 3, "adj P(>C)")
                if adjP_1_2 != 1.0 or adjP_1_3 != 1.0 or adjP_2_3 != 1.0:
                    correlation = False
                else:
                    correlation = True
                # Find frame at which we believe extensive damage has occured.
                frames[j] = run.find_first_n_diff_frames(n, k, P_thresh)
            if self.name == "no_protection":
                merging_thresholds[0] = frames
                adjP_123_correlate[0] = correlation
            else:
                merging_thresholds[conc] = frames
                adjP_123_correlate[conc] = correlation
        return merging_thresholds, adjP_123_correlate

    def get_doses(self, use_frame_nums=False, dose_met="DWD"):
        """Parse the BsxCuBE log file to get the diode readings and then
        calculate the dose for each frame from the diode readings by running
        RADDOSE-3D.
        """
        if use_frame_nums:
            return (np.zeros(self.NUM_FRAMES),
                    np.linspace(1, self.NUM_FRAMES, self.NUM_FRAMES))
        else:
            diode_readings, raddam_onset = self.parse_bsxcube()

            # Convert diode readings to flux and then calculate the dose by
            # running RADDOSE-3D
            doses = {}
            for conc, readings in diode_readings.iteritems():
                dose_values = np.zeros([self.NUM_FRAMES, readings.shape[1]])
                for run in xrange(0, readings.shape[1]):
                    raddose_obj = Raddose3d(self.diode_to_flux(readings)[:, run],
                                            self.EXP_PER_FRAME, dose_met)
                    dose_values[:, run] = raddose_obj.dose_vals
                doses[conc] = self.MGY_TO_KGY * dose_values

            return (diode_readings, doses, raddam_onset)

    def save_doses_and_diode_to_csv(self, dose_dir, diode_dir, rd_onset_dir):
        """Save doses, diode readings and radiation damage onset values
        to CSV files.
        """
        dose_file_dir = "../{}/{}".format(dose_dir, self.name)
        if not os.path.exists(dose_file_dir):
            os.makedirs(dose_file_dir)

        for conc, dose_array in self.doses.iteritems():
            dose_file_name = "{}/doses_conc_{}.csv".format(dose_file_dir, conc)
            np.savetxt(dose_file_name, dose_array, delimiter=",")

        diode_file_dir = "../{}/{}".format(diode_dir, self.name)
        if not os.path.exists(diode_file_dir):
            os.makedirs(diode_file_dir)

        for conc, diode_array in self.diode_readings.iteritems():
            diode_file_name = "{}/doses_conc_{}.csv".format(diode_file_dir,
                                                            conc)
            np.savetxt(diode_file_name, diode_array, delimiter=",")

        rd_onset_file_dir = "../{}/{}".format(rd_onset_dir, self.name)
        if not os.path.exists(rd_onset_file_dir):
            os.makedirs(rd_onset_file_dir)

        for conc, onset_array in self.raddam_onset.iteritems():
            onset_file_name = "{}/doses_conc_{}.csv".format(rd_onset_file_dir,
                                                            conc)
            np.savetxt(onset_file_name, onset_array, delimiter=",")

    def get_saved_doses_and_diode(self, dose_dir, diode_dir, onset_dir):
        """Get the saved doses and diode values
        """
        dose_file_dir = "../{}/{}".format(dose_dir, self.name)
        doses = {}
        for f in os.listdir(dose_file_dir):
            conc = float(f.split("_")[-1].split(".")[0])
            doses[conc] = np.genfromtxt("{}/{}".format(dose_file_dir, f),
                                        delimiter=',')

        diode_file_dir = "../{}/{}".format(diode_dir, self.name)
        diode_readings = {}
        for f in os.listdir(diode_file_dir):
            conc = float(f.split("_")[-1].split(".")[0])
            diode_readings[conc] = np.genfromtxt("{}/{}".format(diode_file_dir, f),
                                                 delimiter=',')

        onset_file_dir = "../{}/{}".format(onset_dir, self.name)
        rd_onset = {}
        for f in os.listdir(onset_file_dir):
            conc = float(f.split("_")[-1].split(".")[0])
            rd_onset[conc] = np.genfromtxt("{}/{}".format(onset_file_dir, f),
                                           delimiter=',')

        return (diode_readings, doses, rd_onset)

    def diode_to_flux(self, diode_readings):
        """Method to convert from diode reading to flux readings
        """
        return self.FLUX_ADD_FAC + self.FLUX_SCALE_FAC * diode_readings

    def parse_bsxcube(self, find_rad_dam=True):
        """Parse the BsxCuBE log file to get the diode readings and radiation
        damage onset information
        """
        # choose only BsxCuBE logs file for which compound present
        dat_file_prefix = self.CMPD_INFO[self.name][self.LIST_INDEX["dat_file_prefix"]]
        BsxCuBE_file_inds = self.CMPD_INFO[self.name][self.LIST_INDEX["bsxcube_log"]]
        BsxCuBE_log_files = ['{}/BsxCuBE.log.{}'.format(self.DATA_LOC_PREFIX, ind) for ind in BsxCuBE_file_inds]

        # create a dictionary 'diode_dic' with keys = concentrations
        # and subkeys = run numbers
        diode_dic = {}
        rad_dam_dic = {}
        for i, conc in enumerate(self.CMPD_CONC):
            run_nums = self.get_run_numbers(concentration=conc,
                                            buffers=False)
            if self.name == "no_protection":
                conc = 0
                if i == 1:
                    break
            diode_dic[conc] = np.zeros([self.NUM_FRAMES, len(run_nums)])
            rad_dam_dic[conc] = np.zeros(len(run_nums))
            for j, run_num in enumerate(run_nums):
                diode_vals, image_nums = [], []
                found_raddam = False
                for k, BsxCuBE_file in enumerate(BsxCuBE_log_files):
                    log_open = open(BsxCuBE_file, 'r')

                    line_num = 0
                    prev_line = ""
                    for l in log_open.readlines():
                        line_num += 1
                        # seach for diode information
                        key_words = ("CollectBrick",
                                     "{}_{}{}".format(dat_file_prefix, '0'*(3-len(str(run_num))), run_num),
                                     ".edf' was collected... (diode:")
                        if all(x in l for x in key_words):
                            diode_vals.append(float(l.split('diode:')[1].split(',')[0]))
                            image_nums.append(int(l.split('.edf')[0].split('_')[-1]))

                        # search for raddam info for current run. Performed by
                        # finding a line reporting radiation damage, checking
                        # the previous line, and keeping raddam value if
                        # previous line corresponds to the required run number
                        if find_rad_dam is True:
                            key_word = "WARNING: Radiation damage detected"
                            if key_word in l:

                                # first line of first file cannot correspond to
                                # required run number
                                if k == 1 and line_num == 1:
                                    continue

                                num_curves_merged = int(l.split('merged')[-1].split()[0])
                                rad_dam_onset = num_curves_merged + 1

                                key_words = ("CollectBrick",
                                             "-log(Fidelity) between")
                                if not all(x in prev_line for x in key_words):
                                    print 'ERROR! Unexpected line contents'
                                dat_file = prev_line.split('and')[-1].split()[0]  # get last .dat file of dataset
                                correct_dat_format = '{}_{}{}_{}{}.dat'.format(dat_file_prefix,
                                                                               '0'*(3-len(str(run_num))),
                                                                               run_num, '0'*(5-len(str(self.NUM_FRAMES))),
                                                                               self.NUM_FRAMES)
                                if correct_dat_format in dat_file:
                                    found_raddam = True
                                    break
                        prev_line = l
                    log_open.close()
                    if found_raddam:
                        break

                # check that correct number of diode values parsed
                # and ensure diode values ordered by image number
                if len(diode_vals) != self.NUM_FRAMES:
                    print 'ERROR'
                image_nums_sorted, diode_vals_sorted = (list(t) for t in zip(*sorted(zip(image_nums, diode_vals))))
                diode_dic[conc][:, j] = np.array(diode_vals_sorted)

                if found_raddam:
                    rad_dam_dic[conc][j] = rad_dam_onset
                elif find_rad_dam:
                    print 'WARNING! No radiation damage limit found for current run'

        return diode_dic, rad_dam_dic

    def plot_diode_readings(self, concentration, run_number, plot_flux=True,
                            display=True, save=False, filename="",
                            directory=""):

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                          SORT Y-VALUES                          #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        if plot_flux:
            y_vals = self.diode_to_flux(self.diode_readings[concentration])
            y_vals = y_vals[:, run_number-1]
        else:
            y_vals = self.diode_readings[concentration][:, run_number-1]

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                    PLOT DIODE/FLUX READINGS                     #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        self.PLOT_NUM += 1
        plt.figure(self.PLOT_NUM)
        plt.plot(np.linspace(1, self.NUM_FRAMES, self.NUM_FRAMES), y_vals, 'o')
        plt.xlabel("Frame number")
        if plot_flux:
            plt.ylabel("Flux (photons/sec)")
        else:
            plt.ylabel("Diode reading (cts/sec)")

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                       SAVE AND/OR DISPLAY                       #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        if save and filename:
            if directory:
                if not os.path.exists("../{}".format(directory)):
                    os.makedirs("../{}".format(directory))
                plot_path = "../{}/{}".format(directory, filename)
            else:
                if not os.path.exists("../{}".format(self.plot_dir)):
                    os.makedirs("../{}".format(self.plot_dir))
                plot_path = "../{}/{}".format(self.plot_dir, filename)
            plt.savefig(plot_path)
        elif save and not filename:
            print "********************** ERROR ***************************"
            print "COULD NOT SAVE PLOT"
            print "No filename specified. Please specify a filename if you"
            print "would like to save the plot."
        if display:
            plt.show()

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
