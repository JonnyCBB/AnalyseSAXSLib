"""CorMap class implementing some analyis of results from the ATSAS suite
program DATCMP
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
    def __init__(self, scat_curve_location):
        file_list = glob.glob(scat_curve_location)
        self.q = np.loadtxt(file_list[0])[:, 0]
        self.I = np.zeros([len(self.q), len(file_list)])
        for i, file in enumerate(file_list):
            frame_data = np.loadtxt(file)
            self.I[:, i] = frame_data[:, 1]
        self.datcmp_data = self.get_datcmp_info(scat_curve_location)

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

    def calc_cormap(self):
        """Return CorMap matrix i.e. the (Pearson product-moment) correlation
        matrix.
        """
        return np.corrcoef(self.I)

    def plot_cormap(self, colour_scheme="gray", display=True, save=False,
                    filename="", directory=""):
        """Create a plot object of a CorMap
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

# --------------------------------------------------------------------------- #
#                               FUNCTIONS                                     #
# --------------------------------------------------------------------------- #


def run_system_command(command_string):
    """Function used to run the system command and return the log"""
    process = subprocess.Popen(command_string, stdout=subprocess.PIPE,
                               shell=True)  # Run system command
    output = process.communicate()  # Get the log.
    return output[0]  # return the log file
