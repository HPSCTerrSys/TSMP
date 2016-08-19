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
)

# list of available versions for a platform
# IMPORTANT: add a leading and trailing " "(space)
availability+=(
        ["JURECA"]=" 1.1.0 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT "
        ["JUQUEEN"]=" 1.1.0 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT "
        ["AGROCLUSTER"]=" 1.1.0 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT "
        ["CCA2"]=" 1.1.0 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT "
        ["CLUMA2"]=" 1.1.0 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0 2.1.0MCT 2.0.5 2.0.5MCT "
)

# list of versions with descriptions
versions+=(
        ["1.2.0"]="1.2.0 (Cosmo5.1) without modifications"
        ["1.2.0MCT"]="1.2.0 (Cosmo5.1) with Oasis3-MCT"
        ["1.1.0"]="1.1.0 without modifications"
        ["1.1.0MCT"]="1.1.0 with Oasis3-MCT"
        ["2.1.0"]="2.1.0 (Cosmo5.1 & CESM) without modifications"
        ["2.1.0MCT"]="2.1.0 (Cosmo5.1 & CESM) with Oasis3-MCT"
        ["2.0.5"]="2.0.5 (Cosmo4.21 & CESM) without modifications"
        ["2.0.5MCT"]="2.0.5 (Cosmo4.21 & CESM) with Oasis3-MCT"
)


# The model versions that correspond to a release version
# order: Oasis, CLM , COSMO, Parflow !!!
# Important: this order must be fulfilled. If one of it is not supported, leave a "" at its place.
modelVersion+=(
        ["1.2.0"]="oasis3 clm3_5 cosmo5_1 parflow"
        ["1.2.0MCT"]="oasis3-mct clm3_5 cosmo5_1 parflow"
        ["1.1.0"]="oasis3 clm3_5 cosmo4_21 parflow"
        ["1.1.0MCT"]="oasis3-mct clm3_5 cosmo4_21 parflow"
        ["2.1.0"]="oasis3 clm4_0 cosmo5_1 parflow"
        ["2.1.0MCT"]="oasis3-mct clm4_0 cosmo5_1 parflow"
        ["2.0.5"]="oasis3 clm4_0 cosmo4_21 parflow"
        ["2.0.5MCT"]="oasis3-mct clm4_0 cosmo4_21 parflow"
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
        ["2.1.0"]=" clm cos pfl clm-cos "
        ["2.1.0MCT"]=" clm cos pfl clm-cos "
        ["2.0.5"]=" clm cos pfl clm-cos "
        ["2.0.5MCT"]=" clm cos pfl clm-cos "

)

#list of supported testcases for a certain machine. (first is default)
setups+=(
        ["cordex"]="444x432 (12km res) atmosphere 436x424 (12km res) land domain of Europe and northern Afrika"
	["nrw"]="150x150 (1km res) atmosphere 300x300 (0.5km res) land domain of North-Rhine-Westphalia"
	["ideal300150"]="idealized domain with gridsize scaled to 150x150 (atmosphere) 300x300 (land)"
	["ideal600300"]="idealized domain with gridsize scaled to 300x300 (atmosphere) 600x600 (land)"
	["ideal1200600"]="idealized domain with gridsize scaled to 600x600 (atmosphere) 1200x1200 (land)"
	["ideal24001200"]="idealized domain with gridsize scaled to 1200x1200 (atmosphere) 2400x2400 (land)"
        ["idealRTD"]="idealized domain 20x20 (atmosphere) 16x16 (land) for land-atmosphere-interaction and DA test"
        ["multi-scale"]="real data simulation over multiple scale Rur subcatchment simulation"
)

# list of setups that are available on a machine. (first is default)
# IMPORTANT: add a leading and trailing " "(space)
setupsAvail+=(
	["JURECA"]=" nrw ideal300150 ideal600300 ideal1200600 ideal24001200 cordex "
        ["JUQUEEN"]=" nrw ideal300150 ideal600300 ideal1200600 ideal24001200 cordex "
	["CLUMA2"]=" nrw idealRTD multi-scale "
	["AGROCLUSTER"]=" nrw "
        ["CCA2"]=" nrw "
)

