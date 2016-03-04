"""CorMap class implementing analyis of results from the ATSAS suite program
DATCMP.
"""
import glob
import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


class ScatterAnalysis(object):
    """CorMap class implementing analysis of data output from a run of the
    ATSAS suite program DATCMP
    """
    # ----------------------------------------------------------------------- #
    #                         CLASS VARIABLES                                 #
    # ----------------------------------------------------------------------- #
    PLOT_LABEL = {'family': 'serif',
                  'weight': 'normal',
                  'size': 16}
    mpl.rc('font', family='serif', weight='normal', size=12)

    PLOT_NUM = 0

    # ----------------------------------------------------------------------- #
    #                         CONSTRUCTOR METHOD                              #
    # ----------------------------------------------------------------------- #
    def __init__(self, scat_curve_location, x_axis_vec=[], x_metric="",
                 x_units=""):
        # Go through files and extract the frame data.
        file_list = glob.glob(scat_curve_location)
        num_frames = len(file_list)
        self.q = np.loadtxt(file_list[0])[:, 0]
        self.I = np.zeros([len(self.q), num_frames])
        for i, file in enumerate(file_list):
            frame_data = np.loadtxt(file)
            self.I[:, i] = frame_data[:, 1]

        # Run DATCMP to get pairwise comparison information.
        self.datcmp_data = self.get_datcmp_info(scat_curve_location)

        # Organise the x-axis used for the plots. Default will be the frame
        # number.
        if not isinstance(x_axis_vec, list):
            if len(x_axis_vec) == num_frames:
                self.x_axis = x_axis_vec
            else:
                print "x_axis_vec is not the same length as the number of"
                print "frames. Using frame numbers instead."
                self.x_axis = np.linspace(1, num_frames, num_frames)
        else:
            self.x_axis = np.linspace(1, num_frames, num_frames)

        if x_metric and x_units:
            self.x_metric = x_metric
            self.x_units = x_units
        else:
            self.x_metric = "Frame number"
            self.x_units = ""

    # ----------------------------------------------------------------------- #
    #                         INSTANCE METHODS                                #
    # ----------------------------------------------------------------------- #

    def get_datcmp_info(self, scattering_curve_files):
        """Method to extract the data produced by from DATCMP
        """
        cmd = "datcmp {}".format(scattering_curve_files)
        log = run_system_command(cmd)
        # define a dictionary to store the data produced from DATCMP - this
        # value will be overwritten.
        data_dict = {"1,2": 0}
        for line in iter(log.splitlines()):
            match_obj = re.match(r'\s* \d{1,} vs', line)
            if match_obj:
                data = line.split()
                if "*" in data[5]:
                    data[5] = data[5][:-1]
                data_dict["{},{}".format(data[0], data[2])] = [int(float(data[3])),
                                                               float(data[4]),
                                                               float(data[5])]
        return data_dict

    def find_diff_frames(self, frame=1, P_threshold=0.01):
        """List all statistically different frames according to the method by
        Daniel Franke, Cy M Jeffries & Dmitri I Svergun (2015)
        """
        if frame <= self.I.shape[1]:
            diff_frames = []
            for i in xrange(0, self.I.shape[1]):
                if i+1 < frame:
                    key = "{},{}".format(i+1, frame)
                elif i+1 > frame:
                    key = "{},{}".format(frame, i+1)
                else:
                    continue
                adjP = self.datcmp_data[key][2]
                if adjP < P_threshold:
                    diff_frames.append(i+1)
            return diff_frames
        else:
            print "********************** ERROR ***************************"
            print "FRAME '{}' DOES NOT EXIST".format(frame)
            print "Use different frame numbers between 1 and {}".format(self.I.shape[1])

    def find_first_n_diff_frames(self, n=1, frame=1, P_threshold=0.01):
        """Return the first frame, F, where there are n consecutive frames
        after F that are also statistically different from the chosen frame.
        """
        # Get list frames that are different
        list_of_diff_frames = self.find_diff_frames(frame, P_threshold)
        if n == 1:
            # If only looking for one frame then return the first value in
            # list of different frames.
            return list_of_diff_frames[0]
        elif n > 1:
            # If we're looking for more than one consecutive frame then we need
            # to keep track of the number of consecutive frames that we've
            # iterated through.
            consec_count = 0
            max_consec_count = 0
            fr_max_count = 0
            for i, curr_fr in enumerate(list_of_diff_frames):
                if i == 0:
                    prev_fr = curr_fr
                    consec_count = 1
                else:
                    if curr_fr == prev_fr + 1:
                        consec_count += 1
                    else:
                        consec_count = 1
                prev_fr = curr_fr
                if consec_count == n:
                    return curr_fr - n + 1
                if consec_count > max_consec_count:
                    max_consec_count = consec_count
                    fr_max_count = curr_fr - max_consec_count + 1
            print "************************ WARNING **************************"
            print "{} consecutive frames not reached!".format(n)
            print "The max number of consecutive frames was {}".format(max_consec_count)
            print "The initial frame for that run was frame {}.".format(fr_max_count)
        else:
            print "********************** ERROR ***************************"
            print "n MUST BE A POSITVE INTEGER VALUE"
            print "User chose n = {}.".format(n)
            print "Please choose a positve integer value for n."

    def similar_frames(self, frame=1, P_threshold=0.01):
        """Return list all of the frames that are similar as defined by the
        method presented in Daniel Franke, Cy M Jeffries & Dmitri I Svergun
        (2015).
        """
        list_of_diff_frames = self.find_diff_frames(frame, P_threshold)
        return [i+1 for i in xrange(0, self.I.shape[1]) if i+1 not in list_of_diff_frames]

    def get_pw_data(self, frame1, frame2, datcmp_data_type="adj P(>C)"):
        """Return C, P(>C) or Bonferroni adjusted P(>C) value from the
        DATCMP output.

        frame1: integer number of the 1st frame used for the pairwise analyis

        frame1: integer number of the 2nd frame used for the pairwise analyis

        datcmp_data_type: string specifying the pairwise result to be returned.
        The input options are:
        1) 'C' - This will return the C value i.e. the max observed patch of
        continuous runs of -1 or 1
        2) 'P(>C)' - This will return the P value of observing a patch of
        continuous runs of -1 or 1 bigger than the corresponding C value.
        3) 'adj P(>C)' - This will return the Bonferroni adjusted P value of
        observing a patch of continuous runs of -1 or 1 bigger than the
        corresponding C value.
        """
        if datcmp_data_type == "C" or datcmp_data_type == 0:
            dat_type = 0
        elif datcmp_data_type == "P(>C)":
            dat_type = 1
        elif datcmp_data_type == "adj P(>C)":
            dat_type = 2
        else:
            print "********************** ERROR ***************************"
            print "INVALID DATCMP DATA TYPE CHOSEN: '{}' DOES NOT EXIST".format(datcmp_data_type)
            print "Please choose either 'C', 'P(>C)' or 'adj P(>C)'."

        if frame1 < frame2:
            datcmp_key = "{},{}".format(frame1, frame2)
        elif frame2 < frame1:
            datcmp_key = "{},{}".format(frame2, frame1)

        if datcmp_key in self.datcmp_data:
            return self.datcmp_data[datcmp_key][dat_type]
        else:
            print "********************** ERROR ***************************"
            print "KEY '{}' DOES NOT EXIST".format(datcmp_key)
            print "Use different frame numbers between 1 and {}".format(self.I.shape[1])

    def calc_cormap(self):
        """Return CorMap matrix i.e. the (Pearson product-moment) correlation
        matrix.
        """
        return np.corrcoef(self.I)

    def calc_pwcormap(self, frame1, frame2):
        """Return the pairwise correlation matrix between frame 1 and frame 2
        """
        pw_I = np.column_stack([self.I[:, frame1-1], self.I[:, frame2-1]])
        return np.corrcoef(pw_I)

    def get_pw_data_array(self, frame=0, delete_zero_row=True):
        """Return an array of all C, P(>C) or Bonferroni adjusted P(>C) values
        from the DATCMP output for the requested frame.
        """
        if frame == 0:
            pw_data = np.zeros([len(self.datcmp_data), 3])
            for i, values in enumerate(self.datcmp_data.itervalues()):
                pw_data[i, :] = np.asarray(values)
        elif 1 <= frame <= self.I.shape[1]:
            pw_data = np.zeros([self.I.shape[1], 3])
            for i in xrange(0, self.I.shape[1]):
                if i+1 < frame:
                    key = "{},{}".format(i+1, frame)
                elif i+1 > frame:
                    key = "{},{}".format(frame, i+1)
                else:
                    continue
                pw_data[i, :] = np.asarray(self.datcmp_data[key])
        else:
            print "********************** ERROR ***************************"
            print "FRAME '{}' DOES NOT EXIST".format(frame)
            print "Use a frame number between 1 and {}".format(self.I.shape[1])

        if delete_zero_row and frame > 0:
            return np.delete(pw_data, (frame-1), axis=0)
        else:
            return pw_data


# ----------------------------------------------------------------------- #
#                        PLOT THE CORRELATION MAP                         #
# ----------------------------------------------------------------------- #
    def plot_cormap(self, colour_scheme="gray", display=True, save=False,
                    filename="", directory=""):
        """Create a CorMap (correlation map) plot
        """
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                          SET PLOT PARAMS                        #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        min_q = min(self.q)
        max_q = max(self.q)
        self.PLOT_NUM += 1

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                          PLOT CORMAP                            #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        plt.figure(self.PLOT_NUM)
        cormap = plt.imshow(self.calc_cormap(), cmap=colour_scheme,
                            extent=[min_q, max_q, min_q, max_q])
        plt.xlabel(r'Scattering Vector, q (nm$^{-1}$)',
                   fontdict=self.PLOT_LABEL)
        plt.ylabel(r'Scattering Vector, q (nm$^{-1}$)',
                   fontdict=self.PLOT_LABEL)
        plt.colorbar(cormap)

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                       SAVE AND/OR DISPLAY                       #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        if save and filename:
            if directory:
                plot_path = "{}/{}".format(directory, filename)
            else:
                plot_path = filename
            plt.savefig(plot_path)
        elif save and not filename:
            print "********************** ERROR ***************************"
            print "COULD NOT SAVE PLOT"
            print "No filename specified. Please specify a filename if you"
            print "would like to save the plot."
        if display:
            plt.show()

# ----------------------------------------------------------------------- #
#                    PLOT THE PAIRWISE CORRELATION MAP                    #
# ----------------------------------------------------------------------- #
    def plot_pwcormap(self, fr1, fr2, colour_scheme="gray", display=True,
                      save=False, filename="", directory=""):
        """Create a pairwise CorMap plot.
        """
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                          SET PLOT PARAMS                        #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        min_q = min(self.q)
        max_q = max(self.q)
        self.PLOT_NUM += 1

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                     PLOT PAIRWISE CORMAP                        #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        plt.figure(self.PLOT_NUM)
        cormap = plt.imshow(self.calc_pwcormap(frame1=fr1, frame2=fr2),
                            cmap=colour_scheme,
                            extent=[min_q, max_q, min_q, max_q])
        plt.xlabel(r'Scattering Vector, q (nm$^{-1}$)',
                   fontdict=self.PLOT_LABEL)
        plt.ylabel(r'Scattering Vector, q (nm$^{-1}$)',
                   fontdict=self.PLOT_LABEL)
        plt.colorbar(cormap)
        adjP = self.get_pw_data(fr1, fr2, "adj P(>C)")
        C = self.get_pw_data(fr1, fr2, "C")
        if self.x_units:
            change_in_x = abs(self.x_axis[fr1-1] - self.x_axis[fr2-1])
            plt.title(r'PW CorMap: frame {} vs {}. C = {}, adj P(>C) = {}, $\Delta${} = {} {}'.format(fr1, fr2, C, adjP, self.x_metric, change_in_x, self.x_units))
        else:
            plt.title("Pairwise CorMap: frame {} vs {}. C = {}, adj P(>C) = {}".format(fr1, fr2, C, adjP))

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                       SAVE AND/OR DISPLAY                       #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        if save and filename:
            if directory:
                plot_path = "{}/{}".format(directory, filename)
            else:
                plot_path = filename
            plt.savefig(plot_path)
        elif save and not filename:
            print "********************** ERROR ***************************"
            print "COULD NOT SAVE PLOT"
            print "No filename specified. Please specify a filename if you"
            print "would like to save the plot."
        if display:
            plt.show()

    # ----------------------------------------------------------------------- #
    #                  PLOT THE HISTOGRAM OF PAIRWISE DATA                    #
    # ----------------------------------------------------------------------- #
    def plot_histogram(self, frame=0, datcmp_data_type="C",
                       display=True, save=False, filename="",
                       directory="", num_bins=20):
        """Plot histogram of pairwise correlation data
        """
        if datcmp_data_type == "C" or datcmp_data_type == 0:
            dat_type = 0
        elif datcmp_data_type == "P(>C)":
            dat_type = 1
        elif datcmp_data_type == "adj P(>C)":
            dat_type = 2
        else:
            print "********************** ERROR ***************************"
            print "INVALID DATCMP DATA TYPE CHOSEN: '{}' DOES NOT EXIST".format(datcmp_data_type)
            print "Please choose either 'C', 'P(>C)' or 'adj P(>C)'."
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                          SET PLOT PARAMS                        #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        self.PLOT_NUM += 1

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                          PLOT HISTOGRAM                         #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        plt.figure(self.PLOT_NUM)
        plt.hist(self.get_pw_data_array(frame)[:, dat_type], bins=num_bins)
        plt.xlabel("{}".format(datcmp_data_type), fontdict=self.PLOT_LABEL)
        plt.ylabel(r'Frequency', fontdict=self.PLOT_LABEL)
        if frame == 0:
            plt.title("{} values for all pairwise comparisons".format(datcmp_data_type))
        else:
            plt.title("{} values for all pairwise comparisons with frame {}".format(datcmp_data_type, frame))

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                       SAVE AND/OR DISPLAY                       #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        if save and filename:
            if directory:
                plot_path = "{}/{}".format(directory, filename)
            else:
                plot_path = filename
            plt.savefig(plot_path)
        elif save and not filename:
            print "********************** ERROR ***************************"
            print "COULD NOT SAVE PLOT"
            print "No filename specified. Please specify a filename if you"
            print "would like to save the plot."
        if display:
            plt.show()

    # ----------------------------------------------------------------------- #
    #                    SCATTER PLOT OF PAIRWISE DATA                        #
    # ----------------------------------------------------------------------- #
    def plot_scatter(self, frame=1, P_threshold=0.01, markersize=60,
                     display=True, save=False, filename="", directory="",
                     legend_loc="upper left", x_change=False, use_adjP=True,
                     xaxis_frame_num=False):
        """Scatter plot of the C values for a chosen frame against all other
        frames.
        """
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                          SET PLOT PARAMS                        #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        self.PLOT_NUM += 1

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                          SCATTER PLOT                           #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        pwframe_data = self.get_pw_data_array(frame=frame,
                                              delete_zero_row=False)

        if xaxis_frame_num:
            x_axis = np.linspace(1, self.I.shape[1], self.I.shape[1])
        else:
            x_axis = self.x_axis

        if x_change:
            sub = x_axis[frame - 1]
            pwframe_data = np.column_stack([abs(x_axis - sub), pwframe_data])
        else:
            pwframe_data = np.column_stack([x_axis, pwframe_data])
        pwframe_data = np.delete(pwframe_data, (frame-1), axis=0)

        colours = ["#0072B2", "#009E73", "#D55E00"]
        if use_adjP:
            lb_dict = {0: [colours[0], "adj P(>C) == 1"],
                       1: [colours[1], "1 >= adj P(>C) >= {}".format(P_threshold)],
                       2: [colours[2], "adj P(>C) < {}".format(P_threshold)]}
            P_col = 3
        else:
            lb_dict = {0: [colours[0], "P(>C) == 1"],
                       1: [colours[1], "1 >= P(>C) >= {}".format(P_threshold)],
                       2: [colours[2], "P(>C) < {}".format(P_threshold)]}
            P_col = 2


        plt.figure(self.PLOT_NUM)
        good_points = pwframe_data[pwframe_data[:, P_col] == 1]
        plt.scatter(good_points[:, 0], good_points[:, 1], color=lb_dict[0][0],
                    s=markersize, edgecolors='#ffffff', alpha=1,
                    label=lb_dict[0][1])

        ok_points = pwframe_data[1 > pwframe_data[:, P_col]]
        ok_points = ok_points[ok_points[:, P_col] >= P_threshold]
        plt.scatter(ok_points[:, 0], ok_points[:, 1], color=lb_dict[1][0],
                    s=markersize, edgecolors='#ffffff', alpha=1,
                    label=lb_dict[1][1])

        bad_points = pwframe_data[pwframe_data[:, P_col] < P_threshold]
        plt.scatter(bad_points[:, 0], bad_points[:, 1], color=lb_dict[2][0],
                    s=markersize, edgecolors='#ffffff', alpha=1,
                    label=lb_dict[2][1])

        plt.legend(loc=legend_loc, scatterpoints=1)
        if x_change:
            if self.x_units:
                plt.xlabel(r'$\Delta${} ({})'.format(self.x_metric, self.x_units),
                           fontdict=self.PLOT_LABEL)
            else:
                plt.xlabel(r'$\Delta${}'.format(self.x_metric),
                           fontdict=self.PLOT_LABEL)
        else:
            if self.x_units:
                plt.xlabel("{} ({})".format(self.x_metric, self.x_units),
                           fontdict=self.PLOT_LABEL)
            else:
                plt.xlabel("{}".format(self.x_metric),
                           fontdict=self.PLOT_LABEL)
        plt.ylabel(r'C', fontdict=self.PLOT_LABEL)
        if xaxis_frame_num:
            plt.title("C values against frame number for frame {}".format(frame))
        else:
            plt.title("C values against {} for frame {}".format(self.x_metric,
                                                                frame))

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                       SAVE AND/OR DISPLAY                       #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        if save and filename:
            if directory:
                plot_path = "{}/{}".format(directory, filename)
            else:
                plot_path = filename
            plt.savefig(plot_path)
        elif save and not filename:
            print "********************** ERROR ***************************"
            print "COULD NOT SAVE PLOT"
            print "No filename specified. Please specify a filename if you"
            print "would like to save the plot."
        if display:
            plt.show()

    # ----------------------------------------------------------------------- #
    #                        PLOT 1D SCATTERING CURVE                         #
    # ----------------------------------------------------------------------- #
    def plot_1d_intensity(self, frames, start_point=1, end_point=-1,
                          log_intensity=False, display=True, save=False,
                          filename="", directory="", legend_loc="upper right"):
        """Plot 1d scatter curves
        """
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                          SET PLOT PARAMS                        #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        self.PLOT_NUM += 1

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                        PLOT SCATTER CURVE                       #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        if end_point == -1:
            end_point = len(self.q)
        reciprocal_resolution = self.q[start_point-1:end_point]
        if isinstance(frames, list):
            frames = [i - 1 for i in frames]
        intensity = self.I[start_point-1:end_point, frames]
        if log_intensity:
            intensity = np.log(intensity)
        plt.figure(self.PLOT_NUM)
        if len(intensity.shape) == 2:
            for i in xrange(0, intensity.shape[1]):
                if self.x_units:
                    plt.plot(reciprocal_resolution, intensity[:, i], 'o',
                             label="Frame {}, {}={} {}".format(frames[i] + 1,
                                                               self.x_metric,
                                                               self.x_axis[frames[i]],
                                                               self.x_units))
                else:
                    plt.plot(reciprocal_resolution, intensity[:, i], 'o',
                             label="Frame {}".format(frames[i] + 1))
        else:
            if self.x_units:
                plt.plot(reciprocal_resolution, intensity, 'o',
                         label="Frame {}, {}={} {}".format(frames,
                                                           self.x_metric,
                                                           self.x_axis[frames-1],
                                                           self.x_units))
            else:
                plt.plot(reciprocal_resolution, intensity, 'o',
                         label="Frame {}".format(frames))
        plt.xlabel(r'Scattering Vector, q ($nm^{-1}$)',
                   fontdict=self.PLOT_LABEL)
        if log_intensity:
            plt.ylabel('log(I) (arb. units.)', fontdict=self.PLOT_LABEL)
        else:
            plt.ylabel('Intensity (arb. units.)', fontdict=self.PLOT_LABEL)
        plt.title('1D Scattering Curve', fontdict=self.PLOT_LABEL)
        plt.legend(loc=legend_loc, scatterpoints=1)

        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        #                       SAVE AND/OR DISPLAY                       #
        # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
        if save and filename:
            if directory:
                plot_path = "{}/{}".format(directory, filename)
            else:
                plot_path = filename
            plt.savefig(plot_path)
        elif save and not filename:
            print "********************** ERROR ***************************"
            print "COULD NOT SAVE PLOT"
            print "No filename specified. Please specify a filename if you"
            print "would like to save the plot."
        if display:
            plt.show()

# --------------------------------------------------------------------------- #
#                               FUNCTIONS                                     #
# --------------------------------------------------------------------------- #


def run_system_command(command_string):
    """Function used to run the system command and return the log"""
    process = subprocess.Popen(command_string, stdout=subprocess.PIPE,
                               shell=True)  # Run system command
    output = process.communicate()  # Get the log.
    return output[0]  # return the log file
