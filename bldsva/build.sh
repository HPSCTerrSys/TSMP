#!/usr/bin/env bash
set -eo pipefail

git submodule update --init eTSMP
source eTSMP/env/jsc.2023_Intel.sh
# Models to build
MODEL_ID="CLM3.5-ParFlow-COSMO"

CLM35_SRC=$(realpath ../CLM3.5)
PARFLOW_SRC=$(realpath ../parflow)
COSMO_SRC=$(realpath ../cosmo5.01_fresh)
ICON_SRC=$(realpath ../icon2.6.4_oascoup)
eCLM_SRC=$(realpath ../eCLM)
OASIS_SRC=$(realpath ../oasis3-mct)
string_cmake=""

if echo "$MODEL_ID" | grep -qE '-'; then

        git submodule update --init --remote $OASIS_SRC
	string_cmake="-DOASIS_SRC=${OASIS_SRC}"

fi

if echo "$MODEL_ID" | grep -qE 'COSMO'; then
	
	git submodule update --init --remote $COSMO_SRC
	cd $COSMO_SRC
	git switch tsmp-oasis
	cd ../bldsva/
        string_cmake="${string_cmake} -DCOSMO_SRC=${COSMO_SRC}"

fi

if echo "$MODEL_ID" | grep -qE 'ParFlow'; then
       
	git submodule update --init --remote $PARFLOW_SRC
	cd $PARFLOW_SRC
	git switch --detach v3.12.0
	cd ../bldsva/
	string_cmake="${string_cmake} -DPARFLOW_SRC=${PARFLOW_SRC}"

fi

if echo "$MODEL_ID" | grep -qE 'CLM3.5'; then

        git submodule update --init --remote $CLM35_SRC
	string_cmake="${string_cmake} -DCLM35_SRC=${CLM35_SRC}"

fi

if echo "$MODEL_ID" | grep -qE 'eCLM'; then

        git submodule update --init --remote $eCLM_SRC
	string_cmake="${string_cmake} -DeCLM_SRC=${eCLM_SRC}"

fi

if echo "$MODEL_ID" | grep -qE 'icon'; then

        git submodule update --init --remote $ICON_SRC
	string_cmake="${string_cmake} -DICON_SRC=${ICON_SRC}"

fi

echo $string_cmake
# Where build artifacts and binaries will be saved
BUILD_DIR="../../bld/${SYSTEMNAME^^}_${MODEL_ID}"
INSTALL_DIR="../../bin/${SYSTEMNAME^^}_${MODEL_ID}"
BUILD_LOG="$(dirname ${BUILD_DIR})/${MODEL_ID}_$(date +%Y.%m.%d_%H.%M).log"
CLM35_SRC="/CLM3.5"
PARFLOW_SRC="/parflow"
COSMO_SRC="/cosmo5.01_fresh"
ICON_SRC="/icon2.6.4_oascoup"
eCLM_SRC="/eCLM"
OASIS_SRC="/oasis3-mct"



cd eTSMP
rm -rf ${BUILD_DIR}
mkdir -p ${BUILD_DIR} ${INSTALL_DIR}
cmake -S . -B ${BUILD_DIR}  $string_cmake   \
      -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      |& tee ${BUILD_LOG}

# Build and install
cmake --build ${BUILD_DIR} |& tee -a ${BUILD_LOG}
cmake --install ${BUILD_DIR} |& tee -a ${BUILD_LOG}

echo ""
echo "Successfully installed ${MODEL_ID} to ${INSTALL_DIR}"
