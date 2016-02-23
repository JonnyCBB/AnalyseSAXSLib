"""Module to perform manipulations and analysis on 1D SAXS scattering curves
"""
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def plot_1d_curve_with_filename(filename, start_point=1, end_point=-1,
                                log_intensity=False):
    """Plot 1D SAXS scatter curve
    """
    font = {'family': 'serif',
            'weight': 'normal',
            'size': 16}
    scatter_curve_data = np.loadtxt(filename)
    reciprocal_resolution = scatter_curve_data[start_point-1:end_point, 0]
    intensity = np.log(scatter_curve_data[start_point-1:end_point, 1])
    if log_intensity:
        intensity = np.log(intensity)
    plt.plot(reciprocal_resolution, intensity, 'o')
    plt.xlabel(r'Scattering Vector, q ($\AA^{-1}$)', fontdict=font)
    if log_intensity:
        plt.ylabel('log(I) (arb. units.)', fontdict=font)
    else:
        plt.ylabel('Intensity (arb. units.)', fontdict=font)
    plt.title('1D Scattering Curve', fontdict=font)
    plt.show()
