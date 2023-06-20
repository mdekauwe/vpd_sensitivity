# VPD sensitivity

## Overview ##

# Estimate leaf VPD from the eddy covariance measurements:

Invert ecosystem conductance (Gs) from the eddy covariance data and use this to estimate the (big) leaf VPD_leaf (eqn 10a, Monteith, 1965; eqn 3, Lin et al. 2018)

$ python src/invert_canopy_gs_vpdl.py

## Datasets

* Eddy covariance dataset: [FLUXNET](http://www.fluxdata.org/DataInfo/default.aspx)


## References

* Monteith, J.L., 1965. Evaporation and environment. Symp. Soc. Exp. Biol. 4.
* Lin, C. et al. 2018. Diel ecosystem conductance response to vapor pressure deficit is suboptimal and independent of soil moisture. Agricultural and Forest Meteorology, 250–251, 24-34.
