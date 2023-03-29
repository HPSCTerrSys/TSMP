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
        ["JURECA"]=" clm3-cosmo4-parflow clm3-cosmo4-parflow-pdaf clm3-cosmo5-parflow clm4-cosmo4-parflow clm4-cosmo5-parflow clm3-cosmo5-parflow-pdaf clm3-icon21-parflow \
                     clm3-icon26-parflow eclm eclm-mct "
        ["JUWELS"]=" clm3-cosmo4-parflow clm3-cosmo4-parflow-pdaf clm3-cosmo5-parflow clm4-cosmo4-parflow clm4-cosmo5-parflow clm3-cosmo5-parflow-pdaf clm3-icon21-parflow \
                     clm3-icon26-parflow eclm eclm-mct "
        ["DEEP"]=" clm3-cosmo4-parflow clm3-cosmo4-parflow-pdaf clm3-cosmo5-parflow clm4-cosmo4-parflow clm4-cosmo5-parflow clm3-cosmo5-parflow-pdaf clm3-icon21-parflow \
                     clm3-icon26-parflow eclm eclm-mct "
        ["JUSUF"]=" clm3-cosmo4-parflow clm3-cosmo4-parflow-pdaf clm3-cosmo5-parflow clm4-cosmo4-parflow clm4-cosmo5-parflow clm3-cosmo5-parflow-pdaf clm3-icon21-parflow \
                     clm3-icon26-parflow eclm eclm-mct "
        ["GENERIC_X86"]=" clm3-cosmo4-parflow clm3-cosmo4-parflow-pdaf clm3-cosmo5-parflow clm4-cosmo4-parflow clm4-cosmo5-parflow clm3-cosmo5-parflow-pdaf clm3-icon21-parflow \
                     clm3-icon26-parflow eclm eclm-mct "
)

# list of versions with descriptions
versions+=(
        ["clm3-cos4-pfl"]=" clm3_5, Cosmo4.21, Oasis3-MCT and Parflow >=3.10"
        ["clm3-cos4-pfl-pdaf"]="with PDAF Data Assimilation"
        ["clm3-cos5-pfl"]="Cosmo5.1 and Parflow >=3.10 with Oasis3-MCT"
        ["clm4-cos4-pfl"]="Cosmo4.21 & CESM and Parflow >=3.10 with Oasis3-MCT"
        ["clm4-cos5-pfl"]="Cosmo5.1 & CESM and Parflow >=3.10 with Oasis3-MCT"
        ["clm3-cos5-pfl-pdaf"]="with PDAF Data Assimilation"
		["clm3-icon21-pfl"]="clm3_5, Parflow >=3.10 and icon2-1 with Oasis3-MCT"
		["clm3-icon26-pfl"]="clm3_5, Parflow >=3.10 and icon2-622 with Oasis3-MCT"
        ["eclm"]="Standalone eCLM"
        ["eclm-mct"]="eCLM with Oasis3-MCT"
)


# The model versions that correspond to a release version
# order: Oasis, CLM , COSMO, Parflow !!!
# Important: this order must be fulfilled. If one of it is not supported, leave a "" at its place.
modelVersion+=(
        ["clm3-cos4-pfl"]="oasis3-mct clm3_5 cosmo4_21 parflow"
        ["clm3-cos4-pfl-pdaf"]="oasis3-mct clm3_5 cosmo4_21 parflow pdaf1_1"
        ["clm3-cos5-pfl"]="oasis3-mct clm3_5 cosmo5_1 parflow"
        ["clm4-cos4-pfl"]="oasis3-mct clm4_0 cosmo4_21 parflow"
        ["clm4-cos5-pfl"]="oasis3-mct clm4_0 cosmo5_1 parflow"
        ["clm3-cos5-pfl-pdaf"]="oasis3-mct clm3_5 cosmo5_1 parflow pdaf1_1"
        ["clm3-icon21-pfl"]="oasis3-mct clm3_5-icon icon2-1 parflow"
		["clm3-icon26-pfl"]="oasis3-mct clm3_5-icon icon2-622 parflow"
        ["eclm"]="eclm cosmo5_1 parflow"
        ["eclm-mct"]="oasis3-mct eclm cosmo5_1 parflow"
)

# list of model combinations that are available for a version. (first is default) 
# usecase: if a platform does not support a model component, create a new version
# with that model limitation/combination and make it available only for this platform 
# IMPORTANT: add a leading and trailing " "(space)
combinations+=(
        ["clm3-cos4-pfl"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["clm3-cos4-pfl-pdaf"]=" clm-cos-pfl clm cos pfl clm-cos clm-pfl "
        ["clm3-cos5-pfl"]=" clm-cos-pfl clm pfl clm-cos clm-pfl "
        ["clm4-cos4-pfl"]=" clm cos pfl clm-cos "
        ["clm4-cos5-pfl"]=" clm cos pfl clm-cos "
        ["clm3-cos5-pfl-pdaf"]=" clm-cos-pfl clm pfl clm-cos clm-pfl "
		["clm3-icon21-pfl"]=" clm-icon-pfl clm icon pfl clm-icon clm-pfl "
		["clm3-icon26-pfl"]=" clm-icon-pfl clm icon pfl clm-icon clm-pfl "
        ["eclm"]=" clm "
        ["eclm-mct"]=" clm "
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

