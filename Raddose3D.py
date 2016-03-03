"""Raddose3D class that handles methods to run RADDOSE-3D, including writing
input files, running RADDOSE-3D and parsing the output.
from
"""


class Raddose3d(object):

    # ----------------------------------------------------------------------- #
    #                    RADDOSE-3D CONSTANT PARAMETERS                       #
    # ----------------------------------------------------------------------- #
    # CRYSTAL
    CRYST_TYPE = "Cylinder"
    DIMS = "1700 1000"
    PPM = 0.05
    COCAL = "SAXSseq"
    SEQ_FILE = "4us6_GI.fasta"
    PROT_CONC = 1.0
    CONT_MAT = "Elemental"
    MAT_EL = "Si 1 O 2"
    CON_THICK = 50
    CON_DENS = 2.648

    # BEAM
    BEAM_TYPE = "Gaussian"
    FWHM = "500 500"
    ENERGY = 12.1
    COLLIMATION = "700 700"

    # WEDGE
    WEDGE = "0 0"

    # ----------------------------------------------------------------------- #
    #                         CONSTRUCTOR METHOD                              #
    # ----------------------------------------------------------------------- #
    def __init__(self, flux_array, exposure_per_frame):
        self.write_raddose3d_input_file(flux_array, exposure_per_frame)
        
    # ----------------------------------------------------------------------- #
    #                         INSTANCE METHODS                                #
    # ----------------------------------------------------------------------- #

    def write_raddose3d_input_file(self, flux_array, exposure_per_frame):
        """Write a RADDOSE-3D input file for a SAXS run
        """
        rd_file = open("raddose3d_input.txt", "w")
        rd_file.write(self.writeCRYSTALBLOCK())
        for i, flux in enumerate(flux_array):
            rd_file.write(self.writeBEAMBLOCK(flux, i+1))
            rd_file.write(self.writeWEDGEBLOCK(exposure_per_frame, i+1))
        rd_file.close()

    def writeCRYSTALBLOCK(self):
        """Method to write the crystal block for the RADDOSE-3D input file.
        """
    	raddose3dinputCRYSTALBLOCK = """
##############################################################################
#                                 Crystal Block                              #
##############################################################################

Crystal
Type {}
Dimensions {}
PixelsPerMicron {}
AbsCoefCalc {}
ProteinConc {}
ContainerMaterialType {}
MaterialElements {}
ContainerThickness {}
ContainerDensity {}

""".format(self.CRYST_TYPE, self.DIMS, self.PPM, self.COCAL, self.SEQ_FILE,
           self.PROT_CONC, self.CONT_MAT, self.MAT_EL, self.CON_THICK,
           self.CON_DENS)
        return raddose3dinputCRYSTALBLOCK

    def writeBEAMBLOCK(self, flux, block_num=""):
        """Method to write the beam block for the RADDOSE-3D input file.
        """
    	raddose3dinputBEAMBLOCK = """
##############################################################################
#                             Beam Block {}                                  #
##############################################################################

Beam
Type {}
Flux
FWHM {}
Energy {}
Collimation Rectangular {}

""".format(block_num, self.BEAM_TYPE, flux, self.FWHM, self.ENERGY,
           self.COLLIMATION)
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
