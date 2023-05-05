
# Terrestrial System Modeling Platform - TSMP

[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/HPSCTerrSys/TSMP/RenderMasterSphinxDocumentation.yml?label=documentation)](https://hpscterrsys.github.io/TSMP/index.html)
[![Latest release](https://img.shields.io/github/v/tag/HPSCTerrSys/TSMP.svg?color=brightgreen&label=latest%20release&sort=semver)](https://github.com/HPSCTerrSys/TSMP/tags) 
[![GitHub last commit](https://img.shields.io/github/last-commit/HPSCTerrSys/TSMP)](https://github.com/HPSCTerrSys/TSMP/commits/master)
[![Twitter Follow](https://img.shields.io/twitter/follow/HPSCTerrSys?style=social)](https://twitter.com/HPSCTerrSys)

## Introduction 

The Terrestrial System Modeling Platform (TSMP or TerrSysMP, https://www.terrsysmp.org) is an open source scale-consistent, highly modular, massively parallel regional Earth system model. TSMP essentially consists of an interface which couples dedicated versions of the Consortium for Small-scale Modeling (COSMO, http://www.cosmo-model.org) or ICOsahedral Nonhydrostatic (ICON, https://code.mpimet.mpg.de/projects/iconpublic) atmospheric model in NWP or climate mode, the Community Land Model (CLM, http://www.cesm.ucar.edu/models/clm/), and the hydrologic model ParFlow (https://www.parflow.org) through the OASIS3-MCT coupler (https://oasis.cerfacs.fr/en/, https://www.mcs.anl.gov/research/projects/mct/).

TSMP allows for a physically-based representation of transport processes of mass, energy and momentum and interactions between the different compartments of the geo-ecosystem across scales, explicitly reproducing feedbacks in the hydrological cycle from the groundwater into the atmosphere.

TSMP is extensively used for idealized and real data process and sensitivity studies in water cycle research, for climate change simulations, data assimilation studies including reanalyses, as well as experimental real time forecasting and monitoring simulations, ranging from individual catchments to continental model domains. TSMP runs on notebooks as well on latest supercomputers using a range of compilers.

TSMP development has been driven by groups within the [Center for High-Performance Scientific Computing in Terrestrial Systems](http://www.hpsc-terrsys.de) (HPSC-TerrSys), as part of the [Geoverbund ABC/J](http://www.geoverbund-abcj.de/geoverbund/EN/Home/home_node.html), the geoscientific network of the University of Cologne, Bonn University, RWTH Aachen University, and the Research Centre Jülich. The current team is anchored in Jülich and Bonn in Germany.

**Visit**

**https://www.terrsysmp.org**

**for information on the features of TSMP, ongoing developments, citation, usage examples, links to documentation, the team, contact information and publications.**

## Quick Start on Linux

Please see [getting started section](https://hpscterrsys.github.io/TSMP/content/gettingstarted.html) for guided steps on how the model can be setup and configured for *one* specific experiment, which we use as one of the default test cases. To get an overview on possible TSMP applications refer to the [TSMP website](https://www.terrsysmp.org) and the [TSMP documention](https://hpscterrsys.github.io/TSMP/index.html).

## TSMP version history
The model components used in TSMP are OASIS3-MCT v2, COSMO v5.01, CLM v3.5, ParFlow 3.2 for TSMP versions v1.2.1, v1.2.2 and v1.2.3, and ParFlow 3.9 for version v1.3.3. TSMP supports ParFlow 3.7 onwards from version v1.3.3 onward.

Those who need to work with ParFlow 3.2, should use the branch `TSMP_pdaf`.

## Citing TSMP

If you use TSMP in a publication, please cite the these papers that describe the model's basic functionalities:

* Shrestha, P., Sulis, M., Masbou, M., Kollet, S., and Simmer, C. (2014). A Scale-Consistent Terrestrial Systems Modeling Platform Based on COSMO, CLM, and ParFlow. Monthly Weather Review, 142(9), 3466–3483. doi:[10.1175/MWR-D-14-00029.1](https://dx.doi.org/10.1175/MWR-D-14-00029.1).
* Gasper, F., Goergen, K., Kollet, S., Shrestha, P., Sulis, M., Rihani, J., and Geimer, M. (2014). Implementation and scaling of the fully coupled Terrestrial Systems Modeling Platform (TerrSysMP) in a massively parallel supercomputing environment &ndash; a case study on JUQUEEN (IBM Blue Gene/Q). Geoscientific Model Development, 7(5), 2531-2543. doi:[10.5194/gmd-7-2531-2014](https://dx.doi.org/10.5194/gmd-7-2531-2014).

## License
TSMP is open source software and is licensed under the [MIT-License](https://github.com/HPSCTerrSys/TSMP/blob/master/LICENSE.txt).

## To come

HPSC-TerrSys will also stage run control scripts for different usage scenarios such as long climate runs, and a variety of pre- and postprocessing tools we developed specifically for TSMP. Documented test case experiments also encompass convection permitting simulations at 1km for COSMO and 0.5km for CLM and ParFlow, idealized experiments, a data assimilation experiment based on the TSMP-PDAF version. We will also update the generic machinefiles to make an adjustment of the built system more straightforward for external users.
