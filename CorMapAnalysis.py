"""CorMap class implementing some analyis of results from the ATSAS suite
program DATCMP
"""
import glob
import subprocess
import re
import numpy as np
# import matplotlib.pyplot as plt


class ScatterAnalysis(object):
    """CorMap class implementing analysis of data output from a run of the
    ATSAS suite program DATCMP
    """
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


def run_system_command(command_string):
    """Function used to run the system command and return the log"""
    process = subprocess.Popen(command_string, stdout=subprocess.PIPE,
                               shell=True)  # Run system command
    output = process.communicate()  # Get the log.
    return output[0]  # return the log file
