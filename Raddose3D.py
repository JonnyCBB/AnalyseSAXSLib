"""Raddose3D class that handles methods to run RADDOSE-3D, including writing
input files, running RADDOSE-3D and parsing the output.
from
"""
import subprocess
import numpy as np
import os


class Raddose3d(object):

    RADDOSE3D_INPUT_FILENAME = "raddose3d_input.txt"
    RADDOSE3D_EXE = "raddose3d.jar"
    RADDOSE3D_OUTPUT_FILES = ["output-DoseState.csv", "output-DoseState.R",
                              "output-Summary.csv", "output-DoseState.txt",
                              "output-Summary.txt"]
    RADDOSE3D_SUMMARY_CSV = "output-Summary.csv"

    CSV_COLUMNS = {"DWD": 1,
                   "Ela": 2,
                   "Diff Eff": 3,
                   "AD-WC": 4,
                   "AD-Exp": 5,
                   "Max": 6,
                   "Dose Th": 7,
                   "Abs En": 8,
                   "TAD": 9,
                   "Dose Con": 10,
                   "Used Vol": 11,
                   "Wedge Abs": 12,
                   "Dose Ineff": 13}
    # ----------------------------------------------------------------------- #
    #                    RADDOSE-3D CONSTANT PARAMETERS                       #
    # ----------------------------------------------------------------------- #
    # CRYSTAL
    CRYST_TYPE = "Cylinder"
    DIMS = "1700 1000"
    PPM = 0.005
    COCAL = "SAXSseq"
    SEQ_FILE = "4us6_GI.txt"
    PROT_CONC = 1.0
    CONT_MAT = "Elemental"
    MAT_EL = "Si 1 O 2"
    CON_THICK = 50
    CON_DENS = 2.648

    # BEAM
    BEAM_TYPE = "ExperimentalPGM"
    BEAM_FILE = "SAXS_beam.pgm"
    ENERGY = 12.5
    PIXELSIZE = "9.836 9.836"

    # WEDGE
    WEDGE = "0 0"

    # ----------------------------------------------------------------------- #
    #                         CONSTRUCTOR METHOD                              #
    # ----------------------------------------------------------------------- #
    def __init__(self, flux_array, exposure_per_frame, concentration,
                 dose_type="DWD", heavy_atom=[]):
        self.write_raddose3d_input_file(flux_array, exposure_per_frame,
                                        concentration, heavy_atom)

        self.run_raddose3d()
        self.dose_vals = self.get_dwd_values(dose_type)
        if len(self.dose_vals) != len(flux_array):
            print '************************* ERROR **************************'
            print "Length of the dose array: {} does not equal length of flux array: {}".format(len(self.dose_vals), len(flux_array))
            print "This needs to be sorted."
        if dose_type == "DWD":
            self.dose_vals = np.divide(self.dose_vals.cumsum(),
                                       np.linspace(1, len(self.dose_vals),
                                                   len(self.dose_vals)))
        self.delete_raddose3d_files()

    # ----------------------------------------------------------------------- #
    #                         INSTANCE METHODS                                #
    # ----------------------------------------------------------------------- #

    def write_raddose3d_input_file(self, flux_array, exposure_per_frame,
                                   concentration, heavy_atom):
        """Write a RADDOSE-3D input file for a SAXS run
        """
        rd_file = open(self.RADDOSE3D_INPUT_FILENAME, "w")
        rd_file.write(self.writeCRYSTALBLOCK(concentration, heavy_atom))
        for i, flux in enumerate(flux_array):
            rd_file.write(self.writeBEAMBLOCK(flux, i+1))
            rd_file.write(self.writeWEDGEBLOCK(exposure_per_frame, i+1))
        rd_file.close()

    def run_raddose3d(self):
        """Run RADDOSE-3D
        """
        terminalCommand = "java -jar {} -i {} -r raddose.exe".format(self.RADDOSE3D_EXE, self.RADDOSE3D_INPUT_FILENAME)
        subprocess.call(terminalCommand)

    def get_dwd_values(self, dose_type="DWD"):
        """Extract DWD values from csv file
        """
        data = np.genfromtxt(self.RADDOSE3D_SUMMARY_CSV, delimiter=',')
        return data[1:, self.CSV_COLUMNS[dose_type]]

    def delete_raddose3d_files(self):
        """Delete the files used/output by RADDOSE-3D
        """
        [os.remove(f) for f in self.RADDOSE3D_OUTPUT_FILES if os.path.isfile(f)]
        if os.path.isfile(self.RADDOSE3D_INPUT_FILENAME):
            os.remove(self.RADDOSE3D_INPUT_FILENAME)

    def writeCRYSTALBLOCK(self, concentration, heavy_atom):
        """Method to write the crystal block for the RADDOSE-3D input file.
        """
    	mainCrystalBlock = """
##############################################################################
#                                 Crystal Block                              #
##############################################################################

Crystal
Type {}
Dimensions {}
PixelsPerMicron {}
AbsCoefCalc {}
SeqFile {}
ProteinConc {}
ContainerMaterialType {}
MaterialElements {}
ContainerThickness {}
ContainerDensity {}
""".format(self.CRYST_TYPE, self.DIMS, self.PPM, self.COCAL, self.SEQ_FILE,
           self.PROT_CONC, self.CONT_MAT, self.MAT_EL, self.CON_THICK,
           self.CON_DENS)
        if heavy_atom:
            heavy_solv_contribution = """SolventHeavyConc {} {}
        """.format(heavy_atom[0], concentration)
        else:
            heavy_solv_contribution = ""
        return "{}{}".format(mainCrystalBlock, heavy_solv_contribution)

    def writeBEAMBLOCK(self, flux, block_num=""):
        """Method to write the beam block for the RADDOSE-3D input file.
        """
    	raddose3dinputBEAMBLOCK = """
##############################################################################
#                             Beam Block {}                                  #
##############################################################################

Beam
Type {}
Flux {}
File {}
Energy {}
PixelSize {}

""".format(block_num, self.BEAM_TYPE, flux, self.BEAM_FILE, self.ENERGY,
           self.PIXELSIZE)
        return raddose3dinputBEAMBLOCK

    def writeWEDGEBLOCK(self, exposure_time, block_num=""):
        """Method to write the wedge block for the RADDOSE-3D input file.
        """
    	raddose3dinputWEDGEBLOCK = """
##############################################################################
#                             Wedge Block {}                                 #
##############################################################################

Wedge {}
ExposureTime {}

""".format(block_num, self.WEDGE, exposure_time)
        return raddose3dinputWEDGEBLOCK
