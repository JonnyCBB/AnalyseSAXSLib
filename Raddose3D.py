"""Raddose3D class that handles methods to run RADDOSE-3D, including writing
input files, running RADDOSE-3D and parsing the output.
from
"""


class Raddose3d(object):

    def writeCRYSTALBLOCK(self):
    	raddose3dinputCRYSTALBLOCK = """
##############################################################################
#                                 Crystal Block                              #
##############################################################################

Crystal

Type %s             # Cuboid or Spherical
Dimensions %s %s %s # Dimensions of the crystal in X,Y,Z in um.
                        # Z is the beam axis, Y the rotation axis and
                        # X completes the right handed set
                        #   (vertical if starting face-on).

PixelsPerMicron %s     # This needs to be at least 10x the beam
                        # FWHM for a Gaussian beam.
                        # e.g. 20um FWHM beam -> 2um voxels -> 0.5 voxels/um

AbsCoefCalc  %s       # Absorption Coefficients Calculated using
                        # RADDOSE v2 (Paithankar et al. 2009)

# Example case for insulin:
UnitCell   78.02  78.02  78.02  # unit cell size: a, b, c
                                # alpha, beta and gamma angles default to 90
NumMonomers  24                 # number of monomers in unit cell
NumResidues  51                 # number of residues per monomer
ProteinHeavyAtoms Zn 2 S 6      # heavy atoms added to protein part of the
                                # monomer, i.e. S, coordinated metals,
                                # Se in Se-Met
SolventHeavyConc P 425          # concentration of elements in the solvent
                                # in mmol/l. Oxygen and lighter elements
                                # should not be specified
SolventFraction 0.64            # fraction of the unit cell occupied by solvent
"""
        return raddose3dinputCRYSTALBLOCK

    def writeBEAMBLOCK(self):
    	raddose3dinputBEAMBLOCK = """
##############################################################################
#                                  Beam Block                                #
##############################################################################

Beam

Type %s             # can be Gaussian or TopHat
Flux %s                 # in photons per second (2e12 = 2 * 10^12)
FWHM %s %s                # in um, vertical by horizontal for a Gaussian beam
Energy %s               # in keV

Collimation Rectangular %s %s # Vertical/Horizontal collimation of the beam
                                # For 'uncollimated' Gaussians, 3xFWHM
                                # recommended
"""
        return raddose3dinputBEAMBLOCK

    def writeWEDGEBLOCK(self):
    	raddose3dinputWEDGEBLOCK = """
##############################################################################
#                                  Wedge Block                               #
##############################################################################

Wedge %s %s                # Start and End rotational angle of the crystal
                          # Start < End
ExposureTime %s           # Total time for entire angular range

# AngularResolution 2     # Only change from the defaults when using very
                          # small wedges, e.g 5.
"""
        return raddose3dinputWEDGEBLOCK
