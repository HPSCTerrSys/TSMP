# Introduction to TSMP & TSMP-PDAF

## TSMP 

The Terrestrial System Modeling Platform (TSMP or TerrSysMP, https://www.terrsysmp.org) is an open source scale-consistent, highly modular, massively parallel regional Earth system model. TSMP essentially consists of an interface which couples dedicated versions of the Consortium for Small-scale Modeling (COSMO, http://www.cosmo-model.org) atmospheric model in NWP or climate mode, the Community Land Model (CLM, http://www.cesm.ucar.edu/models/clm/), and the hydrologic model ParFlow (https://www.parflow.org) through the OASIS3-MCT coupler (https://portal.enes.org/oasis, https://www.mcs.anl.gov/research/projects/mct/).

TSMP allows as a fully integrated soil-vegetation-atmosphere modeling system for a physically-based representation of transport processes of mass, energy and momentum and interactions between the different compartments of the geo-ecosystem across scales, explicitly reproducing feedbacks in the hydrological cycle from the groundwater into the atmosphere.

TSMP is extensively used for idealized and real data process and sensitivity studies in water cycle research, for climate change simulations, data assimilation studies including reanalyses, as well as experimental real time forecasting and monitoring simulations, ranging from individual catchments to continental model domains. TSMP runs on notebooks as well on latest supercomputers using a range of compilers.

TSMP development has been driven by groups within the [Center for High-Performance Scientific Computing in Terrestrial Systems](http://www.hpsc-terrsys.de) (HPSC-TerrSys), as part of the [Geoverbund ABC/J](http://www.geoverbund-abcj.de/geoverbund/EN/Home/home_node.html), the geoscientific network of the University of Cologne, Bonn University, RWTH Aachen University, and the Research Centre Jülich. The current team is anchored in Jülich and Bonn in Germany.

**Visit**

**https://www.terrsysmp.org**

**for information on the features of TSMP, ongoing developments, citation, usage examples, links to documentation, the team, contact information and publications.**

## TSMP-PDAF 

TSMP-PDAF describes the branches of the TSMP-repositories usable for TSMP simulations coupled with the Parallel Data Assimilation Framework ([PDAF](http://pdaf.awi.de/trac/wiki)).

## Web Resources

The main remote repository of TSMP is located on Github:

- Github
  ([https://github.com/HPSCTerrSys/TSMP](https://github.com/HPSCTerrSys/TSMP))

Useful information for users, including a wiki and exercises, can be
found among the resources for the Fall School HPSC Terrestrial
Systems:
[https://gitlab.jsc.fz-juelich.de/sdlts/FallSchool_HPSC_TerrSys](https://gitlab.jsc.fz-juelich.de/sdlts/FallSchool_HPSC_TerrSys).

Additional information on TSMP can be found on the TSMP-website:
[https://www.terrsysmp.org/](https://www.terrsysmp.org/).

### Issues

If you encounter any kind of problem with TSMP or TSMP-PDAF, do not hesitate to open in issue on Github under [https://github.com/HPSCTerrSys/TSMP/issues](https://github.com/HPSCTerrSys/TSMP/issues).

## TSMP: Introduction from Website

TSMP implements a physically-based representation of transport processes of water, energy and momentum across scales down to sub-km
resolution including explicit feedbacks between groundwater, the land surface and atmosphere with a focus on the terrestrial hydrological and energy cycles. TSMP currently comprises three component models:
The Consortium for Small-scale Modeling (COSMO) atmospheric model, the Comunity Land Model (CLM), and the hydrologic model ParFlow ParFlow coupled with OASIS3-MCT coupler and data assimilation capabilities through the Parallel Data Assimilation Framework (PDAF).

TSMP is extensively used for idealized and real data process and sensitivity studies in water cycle research, for climate change
simulations, data assimilation studies including reanalyses, as well as experimental real time forecasting and monitoring simulations, ranging from individual catchments to continental model domains.

## Citing TSMP

If you use TSMP in a publication, please cite the these papers that describe the model's basic functionalities:

* Shrestha, P., Sulis, M., Masbou, M., Kollet, S., and Simmer, C. (2014). A Scale-Consistent Terrestrial Systems Modeling Platform Based on COSMO, CLM, and ParFlow. Monthly Weather Review, 142(9), 3466–3483. doi: [10.1175/MWR-D-14-00029.1](https://dx.doi.org/10.1175/MWR-D-14-00029.1).
* Gasper, F., Goergen, K., Kollet, S., Shrestha, P., Sulis, M., Rihani, J., and Geimer, M. (2014). Implementation and scaling of the fully coupled Terrestrial Systems Modeling Platform (TerrSysMP) in a massively parallel supercomputing environment &ndash; a case study on JUQUEEN (IBM Blue Gene/Q). Geoscientific Model Development, 7(5), 2531-2543. doi: [10.5194/gmd-7-2531-2014](https://dx.doi.org/10.5194/gmd-7-2531-2014).

For TSMP-PDAF additionally 
* Kurtz, W., He, G., Kollet, S. J., Maxwell, R. M., Vereecken, H., & Hendricks Franssen, H. J. (2016). TerrSysMP–PDAF (version 1.0): a modular high-performance data assimilation framework for an integrated land surface–subsurface model. Geoscientific Model Development, 9(4), 1341-1360. doi: [10.5194/gmd-9-1341-2016](http://dx.doi.org/10.5194/gmd-9-1341-2016) 
   

## Relevant publications

Papers describing the model physics or the technical implementation of the TSMP-PDAF components are summarized in the following table:

| Model     | References    | Link (doi)                                                                                                                                 |
|-----------|---------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| ParFlow   | @Ashby1996    | [http://dx.doi.org/10.13182/nse96-a24230](http://dx.doi.org/10.13182/nse96-a24230)                                                         |
|           | @Jones2001    | [http://dx.doi.org/10.1016/s0309-1708(00)00075-0](http://dx.doi.org/10.1016/s0309-1708(00)00075-0)                                         |
|           | @Kollet2006   | [http://dx.doi.org/10.1016/j.advwatres.2005.08.006](http://dx.doi.org/10.1016/j.advwatres.2005.08.006)                                     |
|           | @Maxwell2013  | [http://dx.doi.org/10.1016/j.advwatres.2012.10.001](http://dx.doi.org/10.1016/j.advwatres.2012.10.001)                                     |
| CLM3.5    | @Olefson2004  | [http://dx.doi.org/10.5065/D6N877R0](http://dx.doi.org/10.5065/D6N877R0)                                                                   |
|           | @Oleson2008   | [http://dx.doi.org/10.1029/2007jg000563](http://dx.doi.org/10.1029/2007jg000563)                                                           |
| COSMO     | @Baldauf2011  | [http://dx.doi.org/10.1175/mwr-d-10-05013.1](http://dx.doi.org/10.1175/mwr-d-10-05013.1)                                                   |
| OASIS-MCT | @Valcke2013   | [http://dx.doi.org/10.5194/gmd-6-373-2013](http://dx.doi.org/10.5194/gmd-6-373-2013)                                                       |
|           | @Valcke2013a  | [http://dx.doi.org/https://cerfacs.fr/code_coupling_with_oasis3-mct/](http://dx.doi.org/https://cerfacs.fr/code_coupling_with_oasis3-mct/) |
| PDAF      | @Nerger2013   | [http://dx.doi.org/10.1016/j.cageo.2012.03.026](http://dx.doi.org/10.1016/j.cageo.2012.03.026)                                             |


## Contributors

The following persons contributed signficantly to this documentation before
contributions were made transparent in the git history (alphabetic order):

Fabian Gasper, Abouzar Ghasemi, Johannes Keller, Wolfgang Kurtz, Stefan Poll, Mukund Pondkule, Prabhakar Shrestha

## Additional remark 

This documentnation will be further populated. It is planned to contain an extensive user guide, which describes the TSMP features, the technical implementation of the coupling, the different configurations, adjustments to different HPC systems as well as various test cases incl. all needed forcing datasets is currently under revision and will be available soon with the model system. See the TSMP website for a list of the individual websites and manuals of the component models. To operate TSMP, the functioning and configuration of the component models has to be understood.
