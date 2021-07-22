#! /bin/ksh


# Please edit this lists if you add new versions/platforms to the script
# The keys must not exceed 20 characters !!! 
typeset -A platforms
typeset -A availability
typeset -A versions
typeset -A combinations
typeset -A modelVersion

# list of platforms with descriptions
platforms+=(
        ["CLUMA2"]="MIUB (Uni Bonn) - general purpose Linux Cluster"
        ["AGROCLUSTER"]="IBG3 (FZ-Juelich) - general purpose Linux Cluster"
        ["CCA2"]="ECMWF (Reading, UK) - general purpose Linux Cluster"
        ["JURECA"]="JSC (FZ-Juelich) - general purpose Linux Cluster"
        ["JUQUEEN"]="JSC (FZ-Juelich) - high scale machine"
        ["JUWELS"]="JSC (FZ-Juelich) - general purpose Linux Cluster (big jobs)"
)

# list of available versions for a platform
# IMPORTANT: add a leading and trailing " "(space)
availability+=(
        ["JUWELS"]=" 1.1.0 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT 3.0.0 3.0.0MCT 3.0.0MCTPDAF 3.1.0 3.1.0MCT 1.1.0MCTPDAF "
        ["JURECA"]=" 1.1.0 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT 3.0.0 3.0.0MCT 3.0.0MCTPDAF 3.1.0 3.1.0MCT 1.1.0MCTPDAF "
        ["JUQUEEN"]=" 1.1.0 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT 3.0.0 3.0.0MCT 3.0.0MCTPDAF 3.1.0 3.1.0MCT 1.1.0MCTPDAF "
        ["AGROCLUSTER"]=" 1.1.0 1.1.0MCT 1.1.0MCTPDAF 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT 3.0.0 3.0.0MCT 3.0.0MCTPDAF 3.1.0 3.1.0MCT"
        ["CCA2"]=" 1.1.0 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT 3.0.0 3.0.0MCT 3.1.0 3.1.0MCT "
        ["CLUMA2"]=" 1.1.0 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT 3.0.0 3.0.0MCT 3.1.0 3.1.0MCT "
)

# list of versions with descriptions
versions+=(
        ["1.2.0"]="1.2.0 (Cosmo5.1) without modifications"
        ["1.2.0MCT"]="1.2.0 (Cosmo5.1) with Oasis3-MCT"
        ["1.1.0"]="1.1.0 without modifications"
        ["1.1.0MCT"]="1.1.0 with Oasis3-MCT"
        ["1.1.0MCTPDAF"]="1.1.0 with Oasis3-MCT and PDAF Data Assimilation"
        ["2.1.0"]="2.1.0 (Cosmo5.1 & CESM) without modifications"
        ["2.1.0MCT"]="2.1.0 (Cosmo5.1 & CESM) with Oasis3-MCT"
        ["2.0.5"]="2.0.5 (Cosmo4.21 & CESM) without modifications"
        ["2.0.5MCT"]="2.0.5 (Cosmo4.21 & CESM) with Oasis3-MCT"
        ["3.0.0MCT"]="3.0.0 old models (clm3_5 and cosmo4_32) but new Parflow3_2 and with Oasis3-MCT"
        ["3.0.0MCTPDAF"]="3.0.0 old models (clm3_5 and cosmo4_32) but new Parflow3_2 and with Oasis3-MCT and PDAF Data Assimilation"
        ["3.1.0MCT"]="3.1.0 old clm3_5 but new cosmo5_1 and Parflow3_2 and with Oasis3-MCT"
        ["3.0.0"]="3.0.0 old models (clm3_4 and cosmo4_32) but new Parflow3_2"
 	["3.1.0"]="3.1.0 old clm3_5 but new cosmo5_1 and Parflow3_2"
        ["4.0.0MCT"]="4.0.0 clm3_5-icon icon-lem Parflow3_2 Oasis3-MCT"
        ["4.1.0MCT"]="4.1.0 clm4_0 icon-lem Parflow3_2 Oasis3-MCT"
)


# The model versions that correspond to a release version
# order: Oasis, CLM , COSMO, Parflow, ICON !!!
# Important: this order must be fulfilled. If one of it is not supported, leave a "" at its place.
modelVersion+=(
        ["1.2.0"]="oasis3 clm3_5 cosmo5_1 parflow"
        ["1.2.0MCT"]="oasis3-mct clm3_5 cosmo5_1 parflow"
        ["1.1.0"]="oasis3 clm3_5 cosmo4_21 parflow"
        ["1.1.0MCT"]="oasis3-mct clm3_5 cosmo4_21 parflow"
        ["1.1.0MCTPDAF"]="oasis3-mct clm3_5 cosmo4_21 parflow pdaf1_1"
        ["2.1.0"]="oasis3 clm4_0 cosmo5_1 parflow"
        ["2.1.0MCT"]="oasis3-mct clm4_0 cosmo5_1 parflow"
        ["2.0.5"]="oasis3 clm4_0 cosmo4_21 parflow"
        ["2.0.5MCT"]="oasis3-mct clm4_0 cosmo4_21 parflow"
        ["3.0.0MCT"]="oasis3-mct clm3_5 cosmo4_21 parflow3_2"
        ["3.0.0MCTPDAF"]="oasis3-mct clm3_5 cosmo4_21 parflow3_2 pdaf1_1"
        ["3.0.0"]="oasis3 clm3_5 cosmo4_21 parflow3_2"
        ["3.1.0MCT"]="oasis3-mct clm3_5 cosmo5_1 parflow3_2"
        ["3.1.0"]="oasis3 clm3_5 cosmo5_1 parflow3_2"
        ["4.0.0MCT"]="oasis3-mct clm3_5-icon icon-lem parflow"
        ["4.1.0MCT"]="oasis3-mct clm4_0 icon-lem parflow"
)

# list of model combinations that are available for a version. (first is default) 
# usecase: if a platform does not support a model component, create a new version
# with that model limitation/combination and make it available only for this platform 
# IMPORTANT: add a leading and trailing " "(space)
combinations+=(
        ["1.2.0"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["1.2.0MCT"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["1.1.0"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["1.1.0MCT"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["1.1.0MCTPDAF"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["2.1.0"]=" clm cos pfl clm-cos "
        ["2.1.0MCT"]=" clm cos pfl clm-cos "
        ["2.0.5"]=" clm cos pfl clm-cos "
        ["2.0.5MCT"]=" clm cos pfl clm-cos "
        ["3.0.0MCT"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["3.0.0MCTPDAF"]=" clm-cos-pfl clm pfl clm-cos clm-pfl "
        ["3.0.0"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["3.1.0MCT"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["3.1.0"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["4.0.0MCT"]=" clm-icon-pfl clm icon pfl clm-icon clm-pfl "
        ["4.1.0MCT"]=" clm-icon-pfl clm icon pfl clm-icon clm-pfl "
)

#list of supported testcases for a certain machine.
setups+=(
        ["cordex"]="444x432 (12km res) atmosphere 436x424 (12km res) land domain of Europe and northern Afrika"
	["nrw"]="150x150 (1km res) atmosphere 300x300 (0.5km res) land domain of North-Rhine-Westphalia"
	["ideal300150"]="idealized domain with gridsize scaled to 150x150 (atmosphere) 300x300 (land)"
	["ideal600300"]="idealized domain with gridsize scaled to 300x300 (atmosphere) 600x600 (land)"
	["ideal1200600"]="idealized domain with gridsize scaled to 600x600 (atmosphere) 1200x1200 (land)"
	["ideal24001200"]="idealized domain with gridsize scaled to 1200x1200 (atmosphere) 2400x2400 (land)"
        ["idealRTD"]="idealized domain 20x20 (atmosphere) 16x16 (land) for land-atmosphere-interaction and DA test"
        ["multi-scale"]="real data simulation over multiple scale Rur subcatchment simulation"
        ["idealLES"]="idealized LES runs"
        ["rur"]="Reanalysis over Rur"
        ["bonnRadar"]="two moment microphysics"
        ["seabreeze"]="sea breeze simulation"
        ["smresponse"]="soil moisture response"
        ["bonn"]="flood area of interest"
        ["icon-ccs"]="icon non-hydrostatic convective boundary layer (Anurag et al. 2015)"
        ["nrw-icon"]="150x150 (1.2km res) atmosphere 300x300 (0.5km res) land domain of North-Rhine-Westphalia"
        ["ENS1"]="Neckar Ensemble part1"
        ["ENS11"]="Neckar Ensemble part2"
        ["VR1"]="Neckar VR version 1"
        ["VR2"]="Neckar VR version 2"
	["idealscale"]="Vincents ideal scaling tests"
        ["S3Mnew"]="Intermediate SM Spinup 2018"
        ["S3Mnewagg"]="Intermediate SM Spinup 2018 coarser grid"
        ["S3MnewaggENS"]="Intermediate SM ENSEMBLE part0"
        ["S3MnewaggENS1"]="Intermediate SM ENSEMBLE part1"
        ["S3MnewaggENS2"]="Intermediate SM ENSEMBLE part2"
        ["S3MnewaggENS3"]="Intermediate SM ENSEMBLE part3"
        ["S3MnewaggENSindi"]="Intermediate SM ENSEMBLE individual runs"
        ["S3MENShomoDA"]="Intermediate SM ENSEMBLE for DA"
        ["S3MENShomoDA1"]="Intermediate SM ENSEMBLE for DA part1"
        ["S3MENShomoDA2"]="Intermediate SM ENSEMBLE for DA part2"
        ["S3MENShomoDA3"]="Intermediate SM ENSEMBLE for DA part3"
        ["S3MENShomoDA4"]="Intermediate SM ENSEMBLE for DA part4"
        ["S3MENShomoDAfull_test"]="Intermediate SM ENSEMBLE for DA just for testing"
        ["S3MENShomoDAsel"]="Intermediate SM ENSEMBLE for DA complete"
        ["testcase_day4"]="Fall school test case"
        ["S3MENShomoDAcorr"]="corrected Intermediate SM ENSEMBLE for DA"
)

# list of setups that are available on a machine. (first is default)
# IMPORTANT: add a leading and trailing " "(space)
setupsAvail+=(
	["JUWELS"]=" nrw ideal300150 ideal600300 ideal1200600 ideal24001200 cordex idealRTD multi-scale rur bonnRadar bonn seabreeze smresponse ENS1 ENS11 VR1 VR2 idealscale S3Mnew S3Mnewagg S3MnewaggENS S3MnewaggENSindi S3MnewaggENS1 S3MnewaggENS2 S3MnewaggENS3 S3MENShomoDA S3MENShomoDA1 S3MENShomoDA2 S3MENShomoDA3 S3MENShomoDA4 S3MENShomoDAfull_test S3MENShomoDAsel testcase_day4 S3MENShomoDAcorr "
	["JURECA"]=" nrw ideal300150 ideal600300 ideal1200600 ideal24001200 cordex idealRTD multi-scale rur bonnRadar bonn seabreeze smresponse "
        ["JUQUEEN"]=" nrw ideal300150 ideal600300 ideal1200600 ideal24001200 cordex idealRTD "
	["CLUMA2"]=" nrw idealRTD multi-scale idealLES ENS1 "
	["AGROCLUSTER"]=" nrw "
        ["CCA2"]=" nrw cordex "
)

