
#
# Date: 11 August 2023
#
# Author: Andrew Rigby
#
# Purpose: Code to produce LTE column density cubes from 12CO/13CO cubes
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy import constants as ac
from astropy import units as u
from scipy import interpolate
from scipy.signal import medfilt2d
from scipy.ndimage import binary_dilation, generic_filter
from tqdm import tqdm
from matplotlib.backend_bases import MouseButton
import json
import re
import pkg_resources


# ===== Convenience functions =====

def vaxis(header):
    """
    Purpose: return the velocity axis from a given header
    """
    NAX = header['NAXIS3']
    CDELT = header['CDELT3']
    CRPIX = header['CRPIX3']
    CRVAL = header['CRVAL3']
    vaxis = CDELT * (np.arange(NAX) - CRPIX + 1) + CRVAL

    if 'CUNIT3' in header:
        vunit = u.Unit(header['CUNIT3'])
    else:
        vunit = u.m / u.s

    return vaxis * vunit


def get_einsteinA(filename, J):
    """
    Read the Einstein A coefficient from a LAMDA .dat file
    """
    file = open(filename)
    # read the content of the file opened
    content = file.readlines()
    line = content[J + 50]
    EA = re.split(' +', line)[4]
    return np.float32(EA) / u.s


def get_freq(filename, J):
    """
    Read the rest frequency from a LAMDA .dat file
    """
    file = open(filename)
    # read the content of the file opened
    content = file.readlines()
    line = content[J + 50]
    freq = re.split(' +', line)[5]
    return np.float32(freq) * u.GHz


def get_rms_map(cube):
    """
    Returns the RMS map of a cube by inverting the negative pixel values,
    and assuming that the noise is normally distributed around 0.
    """
    values = cube.copy()
    values[cube.value > 0] = np.nan
    noise_cube = np.concatenate([values, -1 * values])
    rms_map = np.sqrt(np.nanmean(noise_cube**2, axis=0))
    return rms_map


# ==== Functions for LTE calculations, in order of appearance =====

def get_T0(nu):
    """
    Returns the T_0 quantity (h * nu / k_B) for a particular frequency
    """
    return (ac.h * nu / ac.k_B).to('K')


def Jfunc(nu, Tin, Tbg=(2.7 * u.K)):
    """
    The J function used in many LTE equations
    """
    J_Tin = 1 / (np.exp(get_T0(nu) / Tin) - 1)
    J_Tbg = 1 / (np.exp(get_T0(nu) / Tbg) - 1)

    return J_Tin - J_Tbg


def get_Tex(cube, nu, Tbg=2.7 * u.K):
    """
    Returns the excitation temperature cube for an emission cube that
    is assumed to be optically thick
    """
    T0 = get_T0(nu)
    Tex = T0 * (np.log(1 + (T0 / 
                (cube + (T0 / (np.exp(T0 / Tbg) - 1))))))**-1
    return Tex


def get_Tau(cube, nu, Tex_cube):
    """
    Returns a cube of optical depth given an emission cube and accompanying 
    excitation temperature cube
    """
    T0 = get_T0(nu)
    Tau = - np.log(1 - (cube / T0) * Jfunc(nu, Tex_cube)**-1)
    return Tau


def binary_dilate_cube(cube, **kwargs):
    """
    Perform a plane-wise binary dilation on a cube
    """
    dilated = np.zeros_like(cube)
    for zz, plane in enumerate(cube):
        dilated[zz] = binary_dilation(plane, **kwargs)
    return dilated


def fill_cube_interp(cube, emissionmask, **kwargs):
    """
    Interpolate over missing values in each plane of a cube, following:
    https://stackoverflow.com/questions/37662180/interpolate-missing-values-2d-python
    """
    interpolated = cube.copy()

    for pp, plane in enumerate(tqdm(interpolated)):
        if len(np.where(~np.isnan(plane))[0]) >= 4:
            xindices = np.arange(0, plane.shape[1])
            yindices = np.arange(0, plane.shape[0])
            
            #  Mask invalid values
            plane = np.ma.masked_invalid(plane)
            xx, yy = np.meshgrid(xindices, yindices)
            
            # Get only the valid values
            x1 = xx[~plane.mask]
            y1 = yy[~plane.mask]
            newplane = plane[~plane.mask]

            # Perform the interpolation
            newplane = interpolate.griddata((x1, y1), newplane.ravel(),
                                       (xx, yy), **kwargs)  
            
            # Mask to only include values with not-nan value in nanmask cube
            newplane[np.isnan(emissionmask[pp].value)] = np.nan   
            
            # Assign plane to output cube
            interpolated[pp] = newplane   
    
    return interpolated


def fill_cube_median(cube, mask, **kwargs):
    """
    Fills in nan values within a cube using a median filter, and specifying
    an initial mask.
    """
    masked_cube = cube.copy()

    masked_cube[np.isnan(mask)] = np.nan

    # Now clean it up with a median filter
    for i in tqdm(range(np.shape(masked_cube)[0]), desc='Median-filling Tau values'):
        masked_cube[i] = generic_filter(masked_cube[i], 
                                        np.nanmedian, **kwargs)
    masked_cube[np.isnan(mask)] = np.nan

    filled_cube = cube.copy()
    fill_pixels = np.isnan(cube) * np.isfinite(masked_cube)
    filled_cube[fill_pixels] = masked_cube[fill_pixels]

    return filled_cube


def fill_cube_rescale(cube1, cube2, **kwargs):
    """
    Fills missing values in cube1 based on values in cube2, and assuming
    that there is a constant scaling relationship between the two
    """
    # cube1_fill = cube1.copy()
    # ratio = np.nanmedian(cube2.value / cube1.value)
    # cube1_fill[np.isnan(cube1)] = cube2[np.isnan(cube1)].value / ratio

    # # Now clean it up with a median filter
    # for i in tqdm(range(np.shape(cube1_fill)[0]), desc='Extrapolating Tau values'):
    #     cube1_fill[i] = generic_filter(cube1_fill[i], 
    #                                         np.nanmedian, **kwargs)
    
    # cube1_fill[np.isnan(cube2)] = np.nan

    ratio = cube2 / cube1
    wmean = np.nansum(ratio * cube2) / np.nansum(cube2)

    rescaled_cube1 = cube2 / wmean

    filled_cube1 = cube1.copy()
    filled_cube1[np.isnan(filled_cube1)] = rescaled_cube1[np.isnan(filled_cube1)]

    return filled_cube1


def get_columndensity(nu, Tex_cube, Tau_cube, chanwidth, mol, transition,
                      rad_atomic=0.112 * u.nm):
    """
    Arguments:
        chanwidth - channel width
        mol - either '12CO', '13CO', or 'C18O'
        J_lo - J value of the lower energy level
    Returns:
        Column density cube
    """
    if mol == '12CO':
        amu1 = 12
        amu2 = 16
    elif mol == '13CO':
        amu1 = 13
        amu2 = 16
    elif mol == 'C18O':
        amu1 = 12
        amu2 = 18
    else:
        raise Exception('Molecule not recognised')
    
    Jup, Jlo = np.int32(transition.split('-'))

    reduced_mass = (amu1 * amu2) * ac.m_p / (amu1 + amu2)

    # Moment of inertia = reduced mass x mean atomic separation
    inertia = reduced_mass * rad_atomic.to('cm')**2

    # Rotational constant
    B_rot = ac.h / (8 * np.pi**2 * inertia)

    # Partion function (dimensionless)
    Zfunc = ((ac.k_B / (B_rot * ac.h)) * 
             (Tex_cube + (ac.h * B_rot) / (3 * ac.k_B)))
        
    # Einstein A coefficient
    einsteinA = get_einsteinA(f'{DATA_PATH}/{mol.lower()}.dat', Jup)
        
    # Statistical weights
    g_lo = 2 * Jlo + 1
    g_up = 2 * Jup + 1
        
    # Column density in the lower (i.e. [J_up - 1]th) state
    N_lo = ((8 * np.pi / ac.c**3) *
            (g_lo / g_up) * 
            (nu**3 / einsteinA) * 
            (1 - np.exp(-ac.h * nu / (ac.k_B * Tex_cube)))**-1 * 
            Tau_cube)
        
    # Total column density
    N_total = (N_lo *
                (Zfunc / (2 * Jlo + 1)) * 
                np.exp(ac.h * B_rot * Jlo * (Jlo + 1) / 
                       (ac.k_B * Tex_cube))) * chanwidth
    
    return N_total.to('cm-2')


def wmean_map(values, weights):
    return np.nansum(values * weights, axis=0) / np.nansum(weights, axis=0)


def plot_results(velocityplane, velocityaxis,
                 CO12, CO13, CO18, Tex, Tau, coldens):
    """
    Plots some results for sanity checking.
    Arguments:
        velocityplane - if given an integer, just plots the plane with that
                        index. If a float or an astropy quantity, interprets
                        this as a velocity, and locates the appropriate plane.

    """
    if (type(velocityplane) == np.int32) | (type(velocityplane) == np.int64):
        plane = velocityplane
    else:
        plane = np.argmin(np.abs(velocityplane - velocityaxis))

    plt.rcParams['image.origin'] = 'lower'

    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_subplot(321)
    ax1 = ax
    ax.set_title('Tmb_12')
    im = ax.imshow(CO12[plane].value, vmin=0)
    fig.colorbar(im)

    ax = fig.add_subplot(322, sharex=ax1, sharey=ax1)
    ax.set_title('Tex')
    im = ax.imshow(Tex[plane].value, vmin=0, vmax=45)
    fig.colorbar(im)

    ax = fig.add_subplot(323, sharex=ax1, sharey=ax1)
    ax.set_title('Tmb_13')
    im = ax.imshow(CO13[plane].value, vmin=0)
    fig.colorbar(im)

    ax = fig.add_subplot(324,sharex=ax1, sharey=ax1)
    ax.set_title('Tau_13')
    im = ax.imshow(Tau[plane].value, vmin=0, vmax=5)
    fig.colorbar(im)

    ax = fig.add_subplot(325, sharex=ax1, sharey=ax1)
    ax.set_title('Tmb_18')
    im = ax.imshow(CO18[plane].value, vmin=0)
    fig.colorbar(im)

    ax = fig.add_subplot(326, sharex=ax1, sharey=ax1)
    ax.set_title('N_total')
    im = ax.imshow(coldens[plane].value, vmin=1E15, vmax=1E16)
    fig.colorbar(im)

    plt.tight_layout()
    plt.connect('button_press_event', click_spec)

    return fig



def click_spec(event):
    """
    https://matplotlib.org/stable/gallery/event_handling/coords_demo.html
    """
    # print(f'Clicked {event.button} button')
    if event.button is MouseButton.LEFT:
        ypix = int(event.ydata)
        xpix = int(event.xdata)

    print(f'\nValues for pixel {xpix}, {ypix}')
    print(f'12CO: {CO12[testplane, ypix, xpix]}')
    print(f'13CO: {CO13[testplane, ypix, xpix]}')
    print(f'C18O: {CO18[testplane, ypix, xpix]}')
    print(f'Tex_12_13: {Tex_12_13[testplane, ypix, xpix]}')
    print(f'Tau_13_18_fill: {Tau_13_18_fill[testplane, ypix, xpix]}')
    print(f'N_total: {N_total_13[testplane, ypix, xpix]}')


def plot_plane(cube, **kwargs):
    fig = plt.figure()
    plt.imshow(cube[testplane], **kwargs)
    plt.colorbar()

# ===== Main script ======

# if __name__ == "__main__":
#     paramfile = '../testparams.json'

def make_cubes(paramfile):
    plt.ion()
    debug_figures = False

    # ==== Load parameters =====
    # Opening JSON file
    params = json.loads(open(paramfile).read())
    
    data_in_12CO = params['data_in_12CO']
    data_in_13CO = params['data_in_13CO']
    data_in_C18O = params['data_in_C18O']

    data_out = params['data_out']

    snrlim12 = params['snrlim12']
    snrlim13 = params['snrlim13']
    snrlim18 = params['snrlim18']

    Tex_mode = params['Tex_mode']
    single_Tex = params['single_Tex'] * u.K

    Tau_mode = params['Tau_mode']
    R1318 = params['R1318']
    single_Tau = params['single_Tau']
    fill_method = params['fill_method']

    transition = params['transition']
    J_up, J_lo = np.int32(transition.split('-'))

    plot_testplane = params['plot_results']
    Tbg = 2.7 * u.K


    # ===== Load rest frequencies
    global DATA_PATH
    DATA_PATH = pkg_resources.resource_filename('colte', 'data/')
    nu_12 = get_freq(f'{DATA_PATH}/co.dat', J_up)
    nu_13 = get_freq(f'{DATA_PATH}/13co.dat', J_up)
    nu_18 = get_freq(f'{DATA_PATH}/c18o.dat', J_up)

    # ===== Load in data cubes =====
    hdr = fits.getheader(data_in_13CO)
    wcs = WCS(hdr).celestial
    CO12 = fits.getdata(data_in_12CO) * u.K
    CO13 = fits.getdata(data_in_13CO) * u.K
    CO18 = fits.getdata(data_in_C18O) * u.K

    # Identify the brightest 13CO plane for figures
    testplane = np.where(CO13 == np.nanmax(CO13))[0][0]
    vax = vaxis(hdr)

    if 'CUNIT3' in hdr:
        vunit = u.Unit(hdr['CUNIT3'])
    else:
        vunit = u.m / u.s

    dv = (np.abs(hdr['CDELT3']) * vunit).to('km/s')
    
    # ======= Mask cubes based on RMS maps =========
    RMS12 = get_rms_map(CO12)
    RMS13 = get_rms_map(CO13)
    RMS18 = get_rms_map(CO18)

    SNR12 = CO12 / RMS12[np.newaxis, :, :]
    SNR13 = CO13 / RMS13[np.newaxis, :, :]
    SNR18 = CO18 / RMS18[np.newaxis, :, :]

    CO12[SNR12 < snrlim12] = np.nan
    CO13[SNR13 < snrlim13] = np.nan
    CO18[SNR18 < snrlim18] = np.nan

    # ===== First scenario: 13CO coldens from 12CO Tex =====
    if Tex_mode == 1:
        # Excitation temperature from optically thick 12CO (K)
        Tex_12 = get_Tex(CO12, nu_12)

        # Excitation temperature from optically thick 13CO (K)
        Tex_13 = get_Tex(CO13, nu_13)

    elif Tex_mode == 2:
        # Adopt a single excitation temperature
        Tex_12 = np.ones_like(CO12.value) * single_Tex
        Tex_13 = np.ones_like(CO13.value) * single_Tex
    else:
        raise Exception('Tex_mode mot valid')
        
    # Create a Tex map that takes the brighter of the 12CO/13CO values (K)
    # (Assumes 13CO is optically thick where it is brighter than 12CO)
    Tex_max = np.maximum(Tex_12, Tex_13)

    if Tau_mode == 1:
        # Create a mask of pixels that are likely to be problematic
        mask = (1 * CO13 > CO12)

        # Expand by one since edge effects can cause issues
        dilated_mask = binary_dilate_cube(mask, iterations=1)

        # Calculating 13CO Tau directly from excitation temperature.
        # 13CO Optical depth (dimensionless)
        Tau_13 = get_Tau(CO13, nu_13, Tex_max)

        # Calculating 13CO Tau directly from excitation temperature.
        # C18O optical depth (dimensionless)
        Tau_18 = get_Tau(CO18, nu_18, Tex_max)

        # Combine into a hybrid map
        # Use the 13CO optical map, but fill in scaled-up values from C18O where available
        Tau_13_18 = Tau_13.copy()
        Tau_13_18[dilated_mask == 1] = np.nan

        # Converting C18O tau to 13CO assumnig constant abundance because
        # Tau_13 is not always reliable due to non-LTE effects.
        Tau_13_18[~np.isnan(Tau_18)] = R1318 * Tau_18[~np.isnan(Tau_18)]


        # Now fill in missing values
        if fill_method == "interpolate":
            Tau_13_18_fill = fill_cube_interp(Tau_13_18, CO13, method="linear")
        elif fill_method == "median":
            Tau_13_18_fill = fill_cube_median(Tau_13_18, CO13, size=5)
        elif fill_method == "rescale":
            Tau_13_18_fill = fill_cube_rescale(Tau_13_18, CO13, size=3)
        elif fill_method == "none":
            Tau_13_18_fill = Tau_13_18.copy()
    elif Tau_mode == 2:
        print(f'Tau mode {Tau_mode} in development! May give bogus results')
        print('Using low Tau assumption: ')
        Tau0 = np.ones_like(CO13).value * single_Tau
        Tau_13_18_fill = (Tau0 / (1 - np.exp(-Tau0))) * CO13.value
        if single_Tau > 2:
            raise Exception('Cant use low Tau assumption for single_Tau > 2')  
    elif Tau_mode == 3:
        print(f'Tau mode {Tau_mode} in development! May give bogus results')
        print('Using optically thin assumption')
        # In the optically thin limit
        # the integral of tau dv = integral of Tmb dv
        # (Tools of Radio Astronomy, above eq 15.39)
        Tau_13_18_fill = CO13.value
    else:
        raise Exception('Value of optically_thin parameter not recognised. Must be 0 or 1')

    # Now calculate the column density
    N_total_13 = get_columndensity(nu_13, Tex_max, Tau_13_18_fill, 
                                    dv, '13CO', transition)
    
    if plot_testplane == 1:
        plot_results(testplane, vax, CO12, CO13, CO18, Tex_max, 
                    u.Quantity(Tau_13_18_fill), N_total_13)
    
    hdr_out = WCS(hdr).to_header()
    hdr_Tex = hdr_out.copy()
    hdr_Tex['BUNIT'] = 'K'
    hdr_Tau = hdr_out.copy()
    hdr_N = hdr_out.copy()
    hdr_N['BUNIT'] = 'cm-2'

    fits.writeto(f'{data_out}Tex.fits', Tex_max.value, hdr_Tex, overwrite=True)
    fits.writeto(f'{data_out}Tau13.fits', Tau_13_18_fill.value, hdr_Tau, overwrite=True)
    fits.writeto(f'{data_out}N13CO.fits', N_total_13.value, hdr_N, overwrite=True)
