#!/usr/local/epd/bin/python
# -*- coding: utf-8 -*-
#
# $Id$
#
"""
Python module to support kcorrect operations.  Designed to replace IDL
functions.  Currently just supports the solar_magnitudes function.
"""

__author__ = 'Benjamin Weaver <benjamin.weaver@nyu.edu>'

__version__ = '$Revision$'.split(': ')[1].split()[0]

__all__ = [ 'read_basel', 'wavelength_to_edges', 'load_filters',
    'projection_table', 'project_filters', 'solar_magnitudes']

#
# Modules
#
import pylab
import os
import yanny

#
#
#
def read_basel(**kwargs):
    """
    Read a spectrum from a Basel spectrum file.
    """
    if 'solarname' not in kwargs:
        raise TypeError('Invalid keyword arguments passed to read_basel')
    if 'silent' in kwargs:
        silent = kwargs['silent']
    else:
        silent = False
    filename = os.getenv('KCORRECT_DIR') + '/data/basel/' + kwargs['solarname']
    try:
        f = open(filename,'r')
    except IOError:
        print "Error opening %s!" % filename
        raise
    else:
        returndata = {'flux':[], 'model':[], 'teff':[], 'logg':[], 'mh':[],
            'vturb':[], 'xh':[]}
        rawdata = f.read()
        f.close()
        alldata = rawdata.replace("\n", '').split()
        returndata['wavelength'] = pylab.array(map(float,alldata[0:1221]))
        del alldata[0:1221]
        nunits = len(alldata)/(1227)
        if not silent:
            print "%d block(s) of spectra" % nunits
        for u in range(nunits):
            returndata['model'].append(int(alldata[0]))
            returndata['teff'].append(int(alldata[1]))
            returndata['logg'].append(float(alldata[2]))
            returndata['mh'].append(float(alldata[3]))
            returndata['vturb'].append(float(alldata[4]))
            returndata['xh'].append(float(alldata[5]))
            returndata['flux'].append(pylab.array(map(float,alldata[6:1227])))
            del alldata[0:1227]
        return returndata

#
#
#
def wavelength_to_edges(centers):
    """
    Convert set of pixel centers to equivalent edges.
    """
    n = len(centers)
    edges = pylab.zeros((n+1,),centers.dtype)
    edges[1:n] = 0.5*(centers[0:(n-1)] + centers[1:n])
    edges[0] = centers[0] - (edges[1] - centers[0])
    edges[n] = centers[n-1] + (centers[n-1] - edges[n-1])
    return edges

#
#
#
def load_filters(**kwargs):
    """
    Load filter information from a list of files.

    Filter files should contain one structure with columns 'lambda' &
    'pass'.  The name of the table can be arbitrary though.
    """
    if 'filterlist' not in kwargs:
        raise TypeError('Invalid keyword arguments passed to load_filters')
    if 'filterpath' not in kwargs:
        kwargs['filterpath'] = os.getenv('KCORRECT_DIR')+'/data/filters'
    filterdata = {'lambda':[], 'pass':[]}
    for fil in kwargs['filterlist']:
        f = yanny.yanny(kwargs['filterpath']+'/'+fil)
        if str(f) == '':
            raise IOError("Could not read %s" % (kwargs['filterpath']+'/'+fil))
        t = f.tables()
        filterdata['lambda'].append(pylab.array(f[t[0]]['lambda']))
        filterdata['pass'].append(pylab.array(f[t[0]]['pass']))
    return filterdata

#
#
#
def locate(xx, x):
    """
    Find the index of an array whose value at that index is closest to x.
    """
    jl = -1L
    ju = long(len(xx))
    ascnd = xx[-1] > xx[0]
    while ju-jl > 1:
        jm = (ju+jl) >> 1
        if x > xx[jm] and ascnd:
            jl = jm
        else:
            ju = jm
    return jl

#
#
#
def cr_filter(wavelength, filter_wavelength, filter_pass, z, matrix):
    i = locate(filter_wavelength,wavelength)
    if i >= (len(filter_wavelength)-1) or i < 0:
        return 0.0
    ip1 = i + 1
    sl = (wavelength - filter_wavelength[i])/(filter_wavelength[ip1]-filter_wavelength[i])
    filt = filter_pass[i] + sl*(filter_pass[ip1]-filter_pass[i])
    rfw = wavelength/(1.0+z)
    j = locate(wavelength, rfw)
    if i >= (len(wavelength)-1) or i < 0:
        return 0.0
    jp1 = j + 1
    sl = (rfw-wavelength[j])/(wavelength[jp1]-wavelength[j])
    spectrum = matrix[i] + sl*(matrix[jp1]-matrix[j])
    filt = filt*wavelength*spectrum/(1.0+z)
    return filt

#
#
#
def cr_integrate(wavelength, filter_wavelength, filter_pass, z, matrix):
    total = 0.0
    for ip in range(len(wavelength)-1):
        mw = (wavelength[ip+1]+wavelength[ip])*0.5*(1.0+z)
        il = locate(filter_wavelength,mw)
        if il < (len(filter_wavelength)-1) and il >= 0:
            ilp1 = il + 1
            sl = (mw - filter_wavelength[il])/(filter_wavelength[ilp1]-filter_wavelength[il])
            filt = filter_pass[il]+sl*(filter_pass[ilp1] - filter_pass[il])
            dl = pylab.fabs(wavelength[ip+1] - wavelength[ip])
            total += mw*filt*dl*matrix[ip]
    return total

#
#
#
def projection_table(wavelength,flux,**kwargs):
    """
    Create lookup table for calculating maggies corresponding to spectra.
    """
    if 'band_shift' not in kwargs:
        kwargs['band_shift'] = 0.0
    if 'filterpath' not in kwargs:
        kwargs['filterpath'] = os.getenv('KCORRECT_DIR')+'/data/filters'
    if 'silent' in kwargs:
        silent = kwargs['silent']
    else:
        silent = False
    if 'zmin' not in kwargs:
        kwargs['zmin'] = 0.0
    if 'zmax' not in kwargs:
        kwargs['zmax'] = 2.0
    if 'nz' not in kwargs:
        kwargs['nz'] = 1000
    if 'zvals' in kwargs:
        kwargs['nz'] = len(kwargs['zvals'])
    else:
        kwargs['zvals'] = (kwargs['zmin'] +
            (kwargs['zmax'] - kwargs['zmin']) *
            (pylab.arange(kwargs['nz'],dtype=float)+0.5)/float(kwargs['nz']))
    try:
        filterdata = load_filters(**kwargs)
    except IOError:
        print 'Some filters appear to be unreadable!'
        raise
    else:
        if not silent:
            print 'Creating rmatrix ...'
        abfnu = 3.631e-20  # ergs/s/cm^2/Hz
        c = 2.99792458e+18 # angstrom sec^-1
        cl = 0.5*(wavelength[0:len(wavelength)-1] + wavelength[1:len(wavelength)])
        gmatrix = abfnu*c/(cl**2)
        rmatrix = pylab.zeros((kwargs['nz'],len(flux),len(filterdata['lambda']),),dtype=float)
        for k in range(len(filterdata['lambda'])):
            cr_filter_wavelength = filterdata['lambda'][k]/(1.0+kwargs['band_shift'])
            cr_filter_pass = filterdata['pass'][k]
            cr_z = 0.0
            cr_project_matrix = gmatrix
            scale = 1.0/cr_integrate(wavelength,cr_filter_wavelength,cr_filter_pass,cr_z,cr_project_matrix)
            for l in range(kwargs['nz']):
                cr_z = kwargs['zvals'][l]
                for m in range(len(flux)):
                    cr_project_matrix = flux[m]
                    rmatrix[l,m,k] = scale*cr_integrate(wavelength,cr_filter_wavelength,cr_filter_pass,cr_z,cr_project_matrix)
        if not silent:
            print 'Done'
        maggies = rmatrix
        return maggies

#
#
#
def project_filters(wavelength,flux,**kwargs):
    """
    Project a flux onto a set of bandpasses.
    """
    if 'band_shift' not in kwargs:
        kwargs['band_shift'] = 0.0
    if 'filterlist' not in kwargs:
        kwargs['filterlist'] = [("sdss_%s0.par" % f) for f in ('u','g','r','i','z')]
    kwargs['zmin'] = 0.0
    kwargs['zmax'] = 0.0
    kwargs['nz'] = 1
    maggies = projection_table(wavelength,flux,**kwargs)
    return maggies

#
#
#
def solar_magnitudes(**kwargs):
    """
    Calculate the Solar absolute magnitudes in various bandpasses.
    """
    if 'band_shift' not in kwargs:
        kwargs['band_shift'] = 0.0
    if 'solarname' not in kwargs:
        kwargs['solarname'] = 'lcbsun.ori'
    if 'filterlist' not in kwargs:
        kwargs['filterlist'] = [("sdss_%s0.par" % f) for f in ('u','g','r','i','z')]
    #
    # Read in the Sun & put it at 10 pc.
    #
    basel = read_basel(**kwargs)
    nspectra = len(basel['model'])
    angstroms = basel['wavelength']*10.0
    c = 2.99792458e+18 # angstrom sec^-1
    pctocm = 3.085677581e+18 # cm pc^-1
    Rsun = 6.96e+10 # cm
    flux = [ 4.0*pylab.pi*f*c*((0.1*Rsun/pctocm)**2)/(angstroms**2) for f in basel['flux'] ]
    edges = wavelength_to_edges(angstroms)
    maggies = project_filters(edges,flux,**kwargs)
    mag = -2.5*pylab.log10(maggies.reshape((maggies.size,)))
    return mag

#
#
#
def version():
    """
    Returns the output of kcorrect_version.
    """
    return subprocess.Popen("kcorrect_version", stdout=subprocess.PIPE).communicate()[0].strip()

#
#
#
def main():
    """
    Put any tests in this function.
    """
    mag = solar_magnitudes()
    print mag
    return

#
# Run tests here
#
if __name__ == '__main__':
    main()
