#!/usr/bin/env python
"""
Estimate the ecosystem conductance (Gs) from inverting the penman-monteith
against eddy covariance flux data. Finally, make a 1:1 plot of VPD_leaf vs
VPD_atmospheric

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (20.06.2023)"
__email__ = "mdekauwe@gmail.com"

#import matplotlib
#matplotlib.use('agg') # stop windows popping up

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import calendar
import datetime as dt
from scipy.stats import pearsonr
from rmse import rmse
import re
from datetime import datetime


sys.path.append('src')
import constants as c
from penman_monteith import PenmanMonteith
from estimate_pressure import estimate_pressure


def main(fname, hour=False):

    site_name = "Alice_Holt"


    date_parse = lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M:%S')

    df = pd.read_csv(fname, index_col='DateTime',
                     parse_dates=['DateTime'],
                     date_parser=date_parse)
    df.index.names = ['date']


    # Convert units ...

    # hPa -> Pa
    df.loc[:, 'VPD'] *= c.KPA_TO_PA


    # W m-2 to kg m-2 s-1
    lhv = latent_heat_vapourisation(df['Tair'])
    df.loc[:, 'ET'] = df['LE_uStar_f'] / lhv

    # kg m-2 s-1 to mol m-2 s-1
    conv = c.KG_TO_G * c.G_WATER_TO_MOL_WATER
    df.loc[:, 'ET'] *= conv

    # screen for bad data
    df = df[df['Rg'] > -900.0]

    (df, no_G) = filter_dataframe(df, hour)

    if no_G:
        G = None

    df = df.replace('NaN', np.nan)
    df = df.replace('nan', np.nan)
    df = df.replace('#na', np.nan)
    df = df.dropna()


    PM = PenmanMonteith(use_ustar=False)

    # some issue with the wind array that needs checking, one element is a str
    wind = df['wind_speed'].values.astype(float)

    # Height from Wilkinson, M., Eaton, E. L., Broadmeadow, M. S. J., and
    # Morison, J. I. L.: Inter-annual variation of carbon uptake by a
    # plantation oak woodland in south-eastern England, Biogeosciences, 9,
    # 5373–5389, https://doi.org/10.5194/bg-9-5373-2012, 2012.
    (df['Gs'],
     df['VPDl'])  = PM.invert_penman(df['VPD'].values, wind, df['Rnet'].values,
                                     df['Tair'].values, df['Psurf'].values,
                                     df['ET'].values, canht=28., G=G)

    # screen for bad data
    df = df[(df['Gs'] > 0.0) & (np.isnan(df['Gs']) == False)]

    VPDa = df['VPD'] * c.PA_TO_KPA
    VPDl = df['VPDl'] * c.PA_TO_KPA

    plot_vpd(VPDa, VPDl, site_name)

    print(df)

def plot_vpd(VPDa, VPDl, site_name):

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.1)
    fig.subplots_adjust(wspace=0.05)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 14
    plt.rcParams['font.size'] = 14
    plt.rcParams['legend.fontsize'] = 14
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    ax = fig.add_subplot(111)
    ax.set_aspect('equal', adjustable='box')

    ax.plot(VPDa, VPDl, "k.", alpha=0.5)

    one2one = np.array([0, 350])
    ax.plot(one2one, one2one, ls='--', color="grey", label="1:1 line")

    r, pval = pearsonr(VPDa, VPDl)
    #print(r, pval)
    if pval <= 0.05:
        m,c = np.polyfit(VPDa, VPDl, 1)
        ax.plot(VPDa, VPDa*m+c, ls="-", c="red")


    ax.set_ylim(0, 6)
    ax.set_xlim(0, 6)
    ax.set_xlabel("VPD$_a$ (kPa)")
    ax.set_ylabel("VPD$_l$ (kPa)")

    ax.text(3.5, 0.5, 'R$^{2}$ = %0.2f' % r**2)
    ax.text(3.5, 0.25, 'm = %0.2f; c = %0.2f' % (m, c))

    odir = "plots"
    ofname = "%s.pdf" % (site_name)
    fig.savefig(os.path.join(odir, ofname),
                bbox_inches='tight', pad_inches=0.1)


def get_site_info(df_site, fname):

    d = {}
    s = os.path.basename(fname).split(".")[0].split("_")[1].strip()

    d['site'] = s
    d['yrs'] = os.path.basename(fname).split(".")[0].split("_")[5]
    d['lat'] = df_site.loc[df_site.SiteCode == s,'SiteLatitude'].values[0]
    d['lon'] = df_site.loc[df_site.SiteCode == s,'SiteLongitude'].values[0]
    d['pft'] = df_site.loc[df_site.SiteCode == s,\
                           'IGBP_vegetation_short'].values[0]
    d['pft_long'] = df_site.loc[df_site.SiteCode == s,\
                                'IGBP_vegetation_long'].values[0]

    # remove commas from country tag as it messes out csv output
    name = df_site.loc[df_site.SiteCode == s,'Fullname'].values[0]
    d['name'] = name.replace("," ,"")
    d['country'] = df_site.loc[df_site.SiteCode == s,'Country'].values[0]
    d['elev'] = df_site.loc[df_site.SiteCode == s,'SiteElevation'].values[0]
    d['Vegetation_description'] = df_site.loc[df_site.SiteCode == s,\
                                        'VegetationDescription'].values[0]
    d['soil_type'] = df_site.loc[df_site.SiteCode == s,\
                                        'SoilType'].values[0]
    d['disturbance'] = df_site.loc[df_site.SiteCode == s,\
                                        'Disturbance'].values[0]
    d['crop_description'] = df_site.loc[df_site.SiteCode == s,\
                                        'CropDescription'].values[0]
    d['irrigation'] = df_site.loc[df_site.SiteCode == s,\
                                        'Irrigation'].values[0]
    d['measurement_ht'] = -999.9
    try:
        ht = float(df_site.loc[df_site.SiteCode == s, \
                   'MeasurementHeight'].values[0])
        if ~np.isnan(ht):
            d['measurement_ht'] = ht
    except IndexError:
        pass

    d['tower_ht'] = -999.9
    try:
        ht = float(df_site.loc[df_site.SiteCode == s, \
                   'TowerHeight'].values[0])
        if ~np.isnan(ht):
            d['tower_ht'] = ht
    except IndexError:
        pass

    d['canopy_ht'] = -999.9
    try:
        ht = float(df_site.loc[df_site.SiteCode == s, \
                   'CanopyHeight'].values[0])
        if ~np.isnan(ht):
            d['canopy_ht'] = ht
    except IndexError:
        pass

    return (d)


def read_file(fname):

    date_parse = lambda x: datetime.strptime(x, '%Y%m%d%H%M%S')

    df = pd.read_csv(fname, index_col='TIMESTAMP_START',
                     parse_dates=['TIMESTAMP_START'],
                     date_parser=date_parse)
    df.index.names = ['date']

    # Using ERA interim filled met vars ... _F
    df = df.rename(columns={'LE_F_MDS': 'Qle', 'H_F_MDS': 'Qh',
                            'VPD_F_MDS': 'VPD', 'TA_F': 'Tair',
                            'NETRAD': 'Rnet',
                            'G_F_MDS': 'Qg',
                            'WS_F': 'Wind', 'P_F': 'Precip',
                            'USTAR': 'ustar', 'LE_CORR': 'Qle_cor',
                            'H_CORR': 'Qh_cor', 'CO2_F_MDS': 'CO2air',
                            'CO2_F_MDS_QC': 'CO2air_qc', 'PA_F': 'Psurf',
                            'G_F_MDS_QC': 'Qg_qc',
                            'LE_F_MDS_QC': 'Qle_qc', 'H_F_MDS_QC': 'Qh_qc',
                            'LE_CORR_JOINTUNC': 'Qle_cor_uc',
                            'H_CORR_JOINTUNC': 'Qh_cor_uc',
                            'GPP_NT_VUT_REF': 'GPP'})


    df = df[['Qle', 'Qh', 'VPD', 'Tair', 'Rnet', 'Qg', 'Wind', \
             'Precip', 'ustar', 'Qle_cor', 'Qh_cor', 'Psurf',\
             'CO2air', 'CO2air_qc', 'Qg_qc', 'Qle_qc', 'Qh_qc', \
             'Qle_cor_uc', 'Qh_cor_uc', 'GPP']]

    # Convert units ...

    # hPa -> Pa
    df.loc[:, 'VPD'] *= c.HPA_TO_KPA * c.KPA_TO_PA

    # kPa -> Pa
    df.loc[:, 'Psurf'] *= c.KPA_TO_PA

    # W m-2 to kg m-2 s-1
    lhv = latent_heat_vapourisation(df['Tair'])
    df.loc[:, 'ET'] = df['Qle'] / lhv

    # Use EBR value instead - uncomment to output this correction
    #df.loc[:, 'ET'] = df['Qle_cor'] / lhv


    # kg m-2 s-1 to mol m-2 s-1
    conv = c.KG_TO_G * c.G_WATER_TO_MOL_WATER
    df.loc[:, 'ET'] *= conv

    # screen by low u*, i.e. conditions which are often indicative of
    # poorly developed turbulence, after Sanchez et al. 2010, HESS, 14,
    # 1487-1497. Some authors use 0.3 m s-1 (Oliphant et al. 2004) or
    # 0.35 m s-1 (Barr et al. 2006) as a threshold for u*
    df = df[df.ustar >= 0.25]

    # screen for bad data
    df = df[df['Rnet'] > -900.0]

    return (df)

def latent_heat_vapourisation(tair):
    """
    Latent heat of vapourisation is approximated by a linear func of air
    temp (J kg-1)

    Reference:
    ----------
    * Stull, B., 1988: An Introduction to Boundary Layer Meteorology
      Boundary Conditions, pg 279.
    """
    return (2.501 - 0.00237 * tair) * 1E06

def filter_dataframe(df, hour):
    """
    Filter data only using QA=0 (obs) and QA=1 (good)
    """
    no_G = False

    # filter daylight hours
    #
    # If we have no ground heat flux, just use Rn
    no_G = True
    df = df[(df.index.hour >= 7) &
            (df.index.hour <= 18) &
            (df['ET'] > 0.01 / 1000.) & # check in mmol, but units are mol
            (df['VPD'] > 0.05)]

    return (df, no_G)

if __name__ == "__main__":

    fname = "/Users/xj21307/research/Alice_holt/data/alice_holt_met_data.csv"
    main(fname, hour=False)
