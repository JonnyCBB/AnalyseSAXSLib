"""Script used to produce a RADDOSE-3D compatible PGM image from pinhole scans
of the X-ray beam from the SAXS experiment.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from mayavi import mlab
from BeamModule import scaleArray, writePGMFile

# Columns corresponding to the spatial coordinate and the diodebs reading from
# the data recorded after the scan.
spatial_col = 0
diodebs_col = 8

# Get the pin coordinates, and diodebs readings from the scan files
data1 = np.genfromtxt("../pinhole_scan_42.txt")
pin1 = data1[:, spatial_col]
diodebs_20 = data1[:, diodebs_col]
diodebs_21 = np.genfromtxt("../pinhole_scan_43.txt")[:, diodebs_col]
diodebs_19 = np.genfromtxt("../pinhole_scan_44.txt")[:, diodebs_col]
diodebs_17end = np.zeros(len(diodebs_20))
diodebs_23end = np.zeros(len(diodebs_20))
x1 = np.array([1.7, 1.9, 2.0, 2.1, 2.3])
z1 = np.vstack((diodebs_17end, diodebs_19, diodebs_20, diodebs_21,
                diodebs_23end))
spl1 = interpolate.RectBivariateSpline(x1, pin1, z1, kx=3, ky=3, s=0)

data2 = np.genfromtxt("../pinhole_scan_41.txt")
pin2 = data2[:, spatial_col]
diodebs_05 = data2[:, diodebs_col]
diodebs_06 = np.genfromtxt("../pinhole_scan_45.txt")[:, diodebs_col]
diodebs_04 = np.genfromtxt("../pinhole_scan_46.txt")[:, diodebs_col]
diodebs_08end = np.zeros(len(diodebs_06))
diodebs_02end = np.zeros(len(diodebs_06))
x2 = np.array([-0.8, -0.6, -0.5, -0.4, -0.2])
z2 = np.vstack((diodebs_08end, diodebs_06, diodebs_05, diodebs_04,
                diodebs_02end))
spl2 = interpolate.RectBivariateSpline(x2, pin2, z2, kx=2, ky=2, s=0)

# Create spline grid for horizontal scans
horz_grid = spl1(pin2, pin1)  # Knows between 1.7 and 2.3

# Create spline grid for vertical scans
vert_grid = spl2(pin1, pin2)  # Knows between -0.8 and 0.2

# Create weight arrays used for averaging each interpolated beam array
X, Y = np.meshgrid(pin1, pin2)
strong_wt = 0.5  # This value should be 0.5 <= x <= 1.0
horz_grid_wts = 0.5 * np.ones(X.shape)
vert_grid_wts = 0.5 * np.ones(X.shape)
for i in xrange(0, X.shape[0]):
    for j in xrange(0, X.shape[1]):
        if (X[i, j] < -0.6 or X[i, j] > -0.4) and 1.9 <= Y[i, j] <= 2.1:
            vert_grid_wts[i, j] = strong_wt
            horz_grid_wts[i, j] = 1.0 - vert_grid_wts[i, j]
        elif (X[i, j] < -0.6 or X[i, j] > -0.4) and (Y[i, j] < 1.9 or Y[i, j] > 2.1):
            vert_grid_wts[i, j] = 0.5
            horz_grid_wts[i, j] = 0.5
        elif -0.6 <= X[i, j] <= -0.4 and (Y[i, j] < 1.9 or Y[i, j] > 2.1):
            horz_grid_wts[i, j] = strong_wt
            vert_grid_wts[i, j] = 1.0 - horz_grid_wts[i, j]

# Perform weighted averages of the beam arrays
beam_array = np.multiply(horz_grid.transpose(), horz_grid_wts) + \
             np.multiply(vert_grid, vert_grid_wts)

# Scale array for conversion to pgm format
scaled_beam = scaleArray(beam_array)

# Write PGM file
writePGMFile(scaled_beam, "SAXS_beam.pgm")

pixelSizeX = (pin1[-1] - pin1[0])/len(pin1)
pixelSizeY = (pin2[-1] - pin2[0])/len(pin2)

print "Pixel Size X = {:.3f} microns".format(pixelSizeX * 1000)
print "Pixel Size Y = {:.3f} microns".format(pixelSizeY * 1000)

# fig = plt.figure()
# ax = fig.gca(projection='3d')
# surf = ax.plot_surface(X, Y, beam_array, rstride=1, cstride=1,
#                        cmap=cm.coolwarm, linewidth=0, antialiased=True)
# plt.show()

# plt.figure()
# plt.plot(pin2, diodebs_20, '-o')
# plt.show()

mlab.surf(X.transpose(), Y.transpose(), 5e4*beam_array.transpose())
mlab.scalarbar()
mlab.axes()
