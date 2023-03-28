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
        ["GENERIC_X86"]="Generic Linux x86 machine"
        ["JURECA"]="JSC (FZ-Juelich) - general purpose Linux Cluster"
        ["JUWELS"]="JSC (FZ-Juelich) - general purpose Linux Cluster"
        ["DEEP"]="JSC (FZ-Juelich) - general purpose Linux Cluster"
        ["JUSUF"]="JSC (FZ-Juelich) - general purpose Linux Cluster"
)

# list of available versions for a platform
# IMPORTANT: add a leading and trailing " "(space)
availability+=(
        ["JURECA"]=" 1.1.0MCT 1.2.0MCT 2.1.0MCT 2.0.5MCT  3.1.0MCT 3.1.0MCTPDAF 1.1.0MCTPDAF \
                     1.4.0MCT 1.5.0MCT 4.1.0MCT  5.0.0 5.0.0MCT "
        ["JUWELS"]=" 1.1.0MCT 1.2.0MCT 2.1.0MCT 2.0.5MCT  3.1.0MCT  3.1.0MCTPDAF 1.1.0MCTPDAF 
                     1.4.0MCT  1.5.0MCT 1.6.0MCT 4.1.0MCT 5.0.0 5.0.0MCT "
        ["DEEP"]="   1.1.0MCT 1.2.0 1.2.0MCT 2.1.0MCT 2.0.5MCT  3.1.0MCT 1.1.0MCTPDAF 
                     1.4.0MCT  "
        ["JUSUF"]="  1.1.0MCT 1.2.0 1.2.0MCT 2.1.0MCT 2.0.5MCT  3.1.0MCT 1.1.0MCTPDAF 
                     1.4.0MCT  5.0.0 5.0.0MCT "
        ["GENERIC_X86"]=" 1.1.0MCT 1.2.0 1.2.0MCT 2.1.0MCT 2.0.5MCT  3.1.0MCT 1.1.0MCTPDAF "
)

# list of versions with descriptions
versions+=(
        ["1.1.0MCT"]="1.1.0 with Oasis3-MCT and Parflow >=3.10"
        ["1.1.0MCTPDAF"]="1.1.0 with Oasis3-MCT and PDAF Data Assimilation"
        ["1.2.0MCT"]="1.2.0 (Cosmo5.1) and Parflow >=3.10 with Oasis3-MCT"
        ["2.1.0MCT"]="2.1.0 (Cosmo5.1 & CESM) and Parflow >=3.10 with Oasis3-MCT"
        ["2.0.5MCT"]="2.0.5 (Cosmo4.21 & CESM) and Parflow >=3.10 with Oasis3-MCT"
        ["3.1.0MCT"]="3.1.0 clm3_5 cosmo5_1 and Parflow >=3.10 and with Oasis3-MCT"
        ["3.1.0MCTPDAF"]="3.1.0 with Oasis3-MCT and PDAF Data Assimilation"
		["1.4.0MCT"]="1.4.0 clm3_5 and Parflow >=3.10 icon2-1 with Oasis3-MCT"
		["1.5.0MCT"]="1.5.0 clm3_5 new icon version 2.622 and Parflow >=3.10 with Oasis3-MCT"
		["1.6.0MCT"]="1.6.0 clm3_5 new icon version 2.622 and Parflow >=3.10 with Oasis3-MCT v4.0"
        ["4.1.0MCT"]="4.1.0 clm3_5-icon icon2-1 and Parflow >=3.10 with Oasis3-MCT"
        ["5.0.0"]="Standalone eCLM"
        ["5.0.0MCT"]="eCLM with Oasis3-MCT"
)


# The model versions that correspond to a release version
# order: Oasis, CLM , COSMO, Parflow !!!
# Important: this order must be fulfilled. If one of it is not supported, leave a "" at its place.
modelVersion+=(
        ["1.1.0MCT"]="oasis3-mct clm3_5 cosmo4_21 parflow"
        ["1.1.0MCTPDAF"]="oasis3-mct clm3_5 cosmo4_21 parflow pdaf1_1"
        ["1.2.0MCT"]="oasis3-mct clm3_5 cosmo5_1 parflow"
        ["2.0.5MCT"]="oasis3-mct clm4_0 cosmo4_21 parflow"
        ["2.1.0MCT"]="oasis3-mct clm4_0 cosmo5_1 parflow"
        ["3.1.0MCT"]="oasis3-mct clm3_5 cosmo5_1 parflow"
        ["3.1.0MCTPDAF"]="oasis3-mct clm3_5 cosmo5_1 parflow pdaf1_1"
        ["1.4.0MCT"]="oasis3-mct clm3_5-icon icon2-1 parflow"
		["1.5.0MCT"]="oasis3-mct clm3_5-icon icon2-622 parflow"
		["1.6.0MCT"]="oasis3-mct4 clm3_5-icon icon2-622 parflow"
        ["4.1.0MCT"]="oasis3-mct clm3_5-icon icon2-1 parflow"
        ["5.0.0"]="mct eclm cosmo5_1 parflow"
        ["5.0.0MCT"]="oasis3-mct eclm cosmo5_1 parflow"
)

# list of model combinations that are available for a version. (first is default) 
# usecase: if a platform does not support a model component, create a new version
# with that model limitation/combination and make it available only for this platform 
# IMPORTANT: add a leading and trailing " "(space)
combinations+=(
        ["1.2.0MCT"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["1.1.0MCT"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["1.1.0MCTPDAF"]=" clm-cos-pfl clm pfl clm-cos clm-pfl "
        ["2.1.0MCT"]=" clm cos pfl clm-cos "
        ["2.0.5MCT"]=" clm cos pfl clm-cos "
        ["3.1.0MCT"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["3.1.0MCTPDAF"]=" clm-cos-pfl clm pfl clm-cos clm-pfl "
		["1.4.0MCT"]=" clm-icon-pfl clm icon pfl clm-icon clm-pfl "
		["1.5.0MCT"]=" clm-icon-pfl clm icon pfl clm-icon clm-pfl "
		["1.6.0MCT"]=" clm-icon-pfl clm icon pfl clm-icon clm-pfl "
        ["4.1.0MCT"]=" clm-icon-pfl clm icon pfl clm-icon clm-pfl "
        ["5.0.0"]=" clm "
        ["5.0.0MCT"]=" clm "
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
    ["icon-ccs"]="icon non-hydrostatic convective boundary layer (Anurag et al. 2015)"
)

# list of setups that are available on a machine. (first is default)
# IMPORTANT: add a leading and trailing " "(space)
setupsAvail+=(
	["JUWELS"]=" nrw ideal300150 ideal600300 ideal1200600 ideal24001200 cordex idealRTD icon-ccs "
    ["JUSUF"]=" nrw ideal300150 ideal600300 ideal1200600 ideal24001200 cordex idealRTD icon-ccs "
	["JURECA"]=" nrw ideal300150 ideal600300 ideal1200600 ideal24001200 cordex idealRTD icon-ccs "
    ["DEEP"]=" nrw ideal300150 ideal600300 ideal1200600 ideal24001200 cordex idealRTD "
	["GENERIC_X86"]=" nrw ideal300150 ideal600300 ideal1200600 ideal24001200 cordex idealRTD "
)

