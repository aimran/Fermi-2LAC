#!/usr/bin/env python
# -*- coding: utf8 -*-

"""
 This module defines functions for parsing the 2nd-year Fermi LAT AGN Catalog.
 The web page for the catalog is

   http://heasarc.gsfc.nasa.gov/W3Browse/fermi/fermilac.html

 and a full description of the point source analysis and the contents of the
 catalog is provided in http://adsabs.harvard.edu/abs/2011ApJ...743..171A

 The accompanying fits file: fermi_2lac_v01.fits generated using Table 3,
 Table 4, Table 7, Table 8 of the above paper.

 Please contact Asif Imran (asif.imran at icecube.wisc.edu) with questions
 or issues.

 A class is provided to read in the catalog in FITS format and give access to
 the data columns as NumPy arrays.

 ---------------------------------------------------------------------------
 NOTE ON MAKING CUTS:

 We use PyFITS to extract column values from the main catalog data table.  The
 values are extracted into NumPy arrays, making it very easy to apply cuts on
 variables and then apply them across the whole data table using array
 slicing.

 For example, if you want to print the right ascension of all AGN also present
 in the 2FGL catalog, you would do this:

 >>cat = FLACatalog('fermi_2lac_v01.fits')
 >>in2FGL = cat.f2gl_glag
 >>dec = cat.dec[in2FGL]

 Combining several cuts requires use of the logical_and, logical_or, etc.
 functions from the NumPy library.  E.g., to cut on sources present in the
 2FGL catalog as well as declination, do this:

 import numpy
 decCut = cat.dec()
 cut = numpy.logical_and(in2FGL == 1,decCut > 35.0]
 ra = cat.ra[cut]
"""

try:
    from astropy.io import fits as pyfits
    import numpy as np
except ImportError, e:
    print e
    raise SystemExit

class Units:
    """List of base units for the FLAC catalog, in case one wants to convert
    energies or fluxes.
    """

    # Angular
    degree = 3.14159265359 / 180.

    # Energy
    MeV = 1.
    GeV = 1e3*MeV
    TeV = 1e3*GeV
    erg = 0.624150974*TeV

    # Length/area
    cm  = 1.
    cm2 = cm*cm
    m   = 1e2*cm
    m2  = m*m

    # Time
    second = 1.
    minute = 60.*second
    hour   = 60.*minute
    day    = 24.*hour

    # Flux
    inv_cm2sMeV = 1./(cm2*second*MeV)
    inv_cm2sGeV = 1./(cm2*second*GeV)
    inv_cm2sTeV = 1./(cm2*second*TeV)
    inv_m2sTeV  = 1./(m2*second*TeV)
    inv_cm2s    = 1./(cm2*second)
    inv_cm2sErg = 1./(cm2*second*erg)


class FLACatalog:
    """
    Parse the FLAC fits file and extract data from the main table.
    """

    def __init__(self, filename):
        hdulist = pyfits.open(filename)
        tbl = hdulist[1].data  # main table
        try:
            self._sn      = tbl.field("Fermi_Name")
            self._assoc   = tbl.field("CName")
            self._ra      = tbl.field("RAJ2000") * Units.degree
            self._de      = tbl.field("DEJ2000") * Units.degree
            self._l2      = tbl.field("GLON") * Units.degree
            self._b2      = tbl.field("GLAT") * Units.degree
            self._an      = tbl.field("AngSep") * Units.degree
            self._an95    = tbl.field("PosErr95deg") * Units.degree
            self._class1  = tbl.field("Class")
            self._z       = tbl.field("Redshift")
            self._sed     = tbl.field("SED_Type")
            self._bayprob = tbl.field("BayProb")
            self._rgprob  = tbl.field("LR_RGProb")
            self._xgprob  = tbl.field("LR_XGProb")
            self._nsprob  = tbl.field("LogNLogSProb")
            self._rg      = tbl.field("LR_RG")
            self._xg      = tbl.field("LR_XG")
            self._F1e3    = tbl.field("Flux1000_100000") * Units.inv_cm2s
            self._dF1e3   = tbl.field("Unc_Flux1000_100000") * Units.inv_cm2s
            self._PIdx    = tbl.field("Spectral_Index")
            self._dPIdx   = tbl.field("Unc_Spectral_Index")
            self._ts      = tbl.field("TS")
            self._radFl   = tbl.field("Radio_Flux")
            self._xrayFl  = tbl.field("Xray_Flux")
            self._usnov   = tbl.field("USNOVmag")
            self._sdss    = tbl.field("SDSSVmag")
            self._aro     = tbl.field("ARO")
            self._aox     = tbl.field("AOX")
            self._clean   = tbl.field("CLEAN")
            self._in2fgl  = tbl.field("IN_2FGL")

            self._lgnupeak = self.logNuPeak(self.aox, self.aro)

        except KeyError as detail:
            print detail
            pass

    def logNuPeak(self, aox, aro):
        """
        See 3.3.2 of http://iopscience.iop.org/0004-637X/715/1/429/fulltext/
        for a full description of the analytic relatioship
        """

        # Masked array necessary to be able to do arithmatic with missing
        # values without raising runtime errors
        aox = np.ma.masked_array(aox, np.isinf(aox))
        aro = np.ma.masked_array(aro, np.isinf(aro))

        X = 0.565 - 1.433*aro + 0.155*aox
        Y = 1.00 - 0.661*aro - 0.339*aox

        # Define the two cases
        case1 = np.logical_and(X < 0, Y < 0.3)
        case2 = ~case1

        peak = np.zeros_like(X)
        peak[case1] = 13.85 + 2.30*X[case1]
        peak[case2] = 13.15 + 6.58*Y[case2]

        return peak.filled(np.inf)

    @property
    def catalog_name(self):
        """
        Name of form 2FGL JHHMM.m+DDMM.
        """
        return self._sn

    @property
    def counter_part(self):
        """
        Return the name of the most likely association
        """
        return self._assoc

    @property
    def ra(self):
        """
        List of J2000 right ascensions for the counterpart.
        """
        return self._ra

    @property
    def dec(self):
        """
        List of J2000 declinations for the counterpart.
        """
        return self._de

    @property
    def l2(self):
        """
        List of galactic longitudes for the counterpart
        """
        return self._l2

    @property
    def b2(self):
        """
        List of galactic latitudes for the counterpart
        """
        return self._b2

    @property
    def angular_sep(self):
        """
        Angular separation between the 2FGL soure and its counterpart
        in degrees
        """
        return self._an

    @property
    def error_radius(self):
        """
        The 95% Error radius
        """
        return self._an95

    @property
    def source_class(self):
        """
        Primary source class designation.

        1. FSRQ, BL Lac object, radio galazy, steep-spectrum radio quasar
        (SSRQ), Seyfert, NLS1, starburst galazy-for sources with well-
        established classes in literature and/or through an optical
        spectrum with a good evaluation of emission lines.

        2. AGU-for sources without a good optical spectrum or without an
        optical spectrum at all.
         a) BZU objects in the BZCAT list.
         b) Sources in AT20G, VCS, CRATES, FRBA, PMN-CA, CRATES-Gaps, or
            CLASS lists, selected by the log N - log S method (see Section
            3.3) and the LR method (see Section 3.2)
         c) Coincident radio and X-ray sources selected by the LR method (see
            Section 3.2)

        3.AGN-this class is more generic than AGU. These sources are not
        confirmed blazars nor blazar candidates (such as AGU). Although they
        may have had evidence for their flatness in radio emission or
        broadband emission, our intensive optical follow-up program did not
        provide a clear evidence for optical blazar characteristics.

        The LAT 2FLAC source classes are described in detail here:
        http://iopscience.iop.org/0004-637X/743/2/171/fulltext/

        Available types are:
        'FSRQ'
        'Radio Galaxy'
        'SSRQ'
        'Starburst galaxy'
        'BL Lac'
        'AGU'
        'AGN'
        'Unidentified'
        """
        return self._class1

    @property
    def redshift(self):
        """
        Redshift of the source
        """
        return self._z

    @property
    def sed_type(self):
        """
        The spectral energy distribution class (based on the synchrotron
        peak frequency)
        """
        return self._sed

    @property
    def bayesian_prob(self):
        """
        Posterior probability that a source from a catalog of candidate
        counterparts is truly an emitter of gamma-rays detected by the LAT
        """
        return self._bayprob

    @property
    def rg_prob(self):
        """
        Reliability value that radio counterpart and gamma ray values match
        """
        return self._rgprob

    @property
    def xg_prob(self):
        """
        Reliability value that x-ray counterpart and gamma ray values match
        """
        return self._xgprob

    @property
    def logNS_prob(self):
        """
        Probability of association between gamma-ray source and a candidate
        counterpart using the logN-LogS method (see text for detail)
        """
        return self._nsprob

    @property
    def rg_ll_ratio (self):
        """
        Likelihood ratio that radio counterpart and gamma ray values match
        """
        return self._rg

    @property
    def xg_ll_ratio (self):
        """
        Likelihood ratio that x-ray counterpart and gamma ray values match
        """
        return self._xg

    @property
    def flux1000(self):
        """
        Photon flux for 1 GeV - 100 GeV obtained by summing the photon flux
        values form the likelihood analysis in the three bands from 1 GeV
        to 100 GeV.  If dFlux/Flux exceeds 0.5, then Flux + 2*dFlux is
        given instead as an approximate 2-sigma upper limit.
        """
        return self._F1e3

    @property
    def flux1000_err (self):
        """
        One-sigma error on integral flux from 1 to 100 GeV obtained from
        summing in quadrature the errors from the three energy bands from 1
        GeV to 100 GeV.  If dFlux/Flux > 0.5, this is set to zero.
        """
        return self._dF1e3

    @property
    def spectral_index(self):
        """
        Best fit for the photon number power-law index.  For LogParabola
        spectra, index at Pivot Energy; for PLExpCutoff spectra, low-energy
        index.  Derived from the likelihood analysis for 100 MeV - 100 GeV.
        """
        return self._PIdx

    @property
    def spectral_index_error(self):
        """
        Get the 1-sigma error on the spectral index.
        """
        return self._dPIdx

    @property
    def ts(self):
        """
        The test statistic is defined as TS = 2(log L(source)-logL(nosource)),
        where L represents the likelihood of the data given the model with or
        without a source present at a given position on the sky.
        """
        return self._ts

    @property
    def radio_flux(self):
        """
        Flux from the source in radio band
        (in units of mili-Jansky)
        """
        return self._radFl

    @property
    def xray_flux(self):
        """
        Flux from the gamma-ray source region in X-ray band
        (in units of erg/cm2/sec)
        """
        return self._xrayFl

    @property
    def usno_mag(self):
        """
        Magnitude of the optical flux (taken from USNO-B1.0)
        """
        return self._usnov

    @property
    def sdss_mag(self):
        """
        Magnitude of the optical flux (taken from Sloan Digital Sky Survey)
        """
        return self._sdss

    @property
    def aox(self):
        """
        Broadband spectral index between optical and X-ray band (between
        5 GHz and 5000 Angstrom), used for determining the synchrotron-peak
        frequency.
        """
        return self._aox

    @property
    def aro(self):
        """
        Broadband spectral index between radio and optical band (between
        5000 Angstrom and 1 KeV), used for determining the synchrotron-peak
        frequency.
        """
        return self._aro

    @property
    def log_nu_peak(self):
        """
        The log of the synchrotron peak (in Hz)
        """
        return self._lgnupeak

    @property
    def clean_flag(self):
        """
        Flag for a source in the Clean sample. See text for detail but based on
        the 2FGL analysis flag
        1: Clean 0:No
        """
        return self._clean

    @property
    def f2gl_flag(self):
        """
        Flag for a source in the 2FGL List
        1: Yes 0:No
        """
        return self._in2fgl
