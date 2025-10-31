# Obtaining the code 

## TSMP

The source code of TSMP (Terrestrial System Modeling Platform) can be obtained from a git-repository hosted on Github:
<https://github.com/HPSCTerrSys/TSMP>. On this website, a README file explains how to run TSMP for the first time on Linux.

This git-repository provides wiki pages on how to download, configure and install TSMP. The latter two points are also explained in this manual in Sections [Build Examples](./build_examples.md) and [Setup](./../setup_tsmp/setup_examples.md) for convenience. Additionally, a few modifications need to be done for TSMP in order to work within TSMP-PDAF (*this patching is possibly not necessary anymore in the current version*).

Commands for cloning from these repositories (https or ssh):

```bash
	git clone https://github.com/HPSCTerrSys/TSMP TSMP
	git clone git@github.com:HPSCTerrSys/TSMP.git TSMP
```

More information about TSMP is available from the webpage:
<https://www.terrsysmp.org/>
and from the home page of the Centre for High-Performance Scientific Computing in Terrestrial Systems (HPSC-TerrSys):
[http://www.hpsc-terrsys.de](http://www.hpsc-terrsys.de)

There is also a public TSMP Virtual Machine available for download. This is a Ubuntu virtual machine containing TSMP and a tutorial on data assimilation. The webpage for downloading (TSMP Virtual Machine): <https://datapub.fz-juelich.de/slts/>

## TSMP-PDAF

Commands for cloning TSMP-PDAF (https or ssh):
```bash
	git clone -b TSMP-pdaf https://github.com/HPSCTerrSys/TSMP TSMP
	git clone -b TSMP-pdaf git@github.com:HPSCTerrSys/TSMP.git TSMP
```

A second git-repository contains internal branches with current developments and is hosted by the Institute of Bio- and Geosciences
Agrosphere (IBG-3) at the Forschungszentrum Jülich (registration required, contact jo.keller@fz-juelich.de):
<https://icg4geo.icg.kfa-juelich.de>. 
```bash
	git clone https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/tsmp.git TSMP
	git clone git@icg4geo.icg.kfa-juelich.de:ExternalRepos/tsmp-pdaf/tsmp.git TSMP
```

## Component Models

The Github-page of TSMP does not include the source code of component models. This makes TSMP independent of third party code and licenses.

Thus, with:

``` bash
	git clone https://github.com/HPSCTerrSys/TSMP TSMP
```

or (deprecated)

``` bash
	git clone https://git2.meteo.uni-bonn.de/git/terrsysmp TSMP
```

The standard root folder in this document is called `TSMP`.

You are creating your root-folder `TSMP` including the subfolder `TSMP/bldsva` that holds the entire TSMP interface.

All the source codes of component models should ideally be placed into one installation directory (for example the root folder `TSMP`).

The directory structure of `TSMP` is discussed in more detail in Section [Structuring TSMP-PDAF](./structure.md).

### HPSC-TerrSys users

Authenticate with your GitLab web GUI user name and password and clone
the repositories (instead of "fresh", also "legacy" repositories with
specific code modifications may be retrieved):

```shell
git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/cosmo5.01_fresh.git  cosmo5_1
git clone -b v3.9.0 https://github.com/parflow/parflow.git                              parflow
git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/clm3.5_fresh.git     clm3_5
git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/oasis3-mct.git       oasis3-mct
```

### External users

Make sure directories are renamed as shown above. The component models can be downloaded from the respective websites as indicated below. Due to licensing reasons HPSC-TerrSys cannot provide these source codes to external users.

#### COSMO v5.01

Available from http://www.cosmo-model.org. A license agreement is needed.

#### CLM v3.5

Available from https://github.com/HPSCTerrSys/CLM3.5/tree/clm3.5_rel3.code.c070524.

#### ParFlow v3.2

Available from https://github.com/parflow/.

```shell
git clone --branch v3.2.0 https://github.com/parflow/parflow.git parflow3_2
```
#### ParFlow v3.9
ParFlow 3.9  is available from
```shell
git clone -b v3.9.0 https://github.com/parflow/parflow.git
```

#### OASIS3-MCT v2.0

Available from https://oasis.cerfacs.fr/en/downloads/.

#### depracted?

Include obtaining from original sources and from Gitlab page
<https://icg4geo.icg.kfa-juelich.de/ModelSystems>.

For information on component models of TSMP-PDAF, see
[https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/TSMP-PDAF-component-models](https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/TSMP-PDAF-component-models)


#### PDAF

The source code of PDAF is available from the following web page under Content-link `Software download` (registration required):
<http://pdaf.awi.de/trac/wiki>. If possible, ask a senior institute member for existing versions on JSC systems.

Install instructions are provided in the README files included in the downloadable archives. The installation procedure is currently not described in this user guide. The current version of TSMP-PDAF is build with PDAF version 1.13.2.

The source code of TSMP-PDAF is available form the HPSC-TerrSys web site: [http://www.hpsc-terrsys.de](http://www.hpsc-terrsys.de)\ The
new GIT-project "TSMP\" does not include the source code of the component models anymore.

