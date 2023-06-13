#!/usr/bin/env bash
set -eo pipefail

# Models to build
MODEL_ID="CLM3.5-ParFlow-COSMO"
CLM35_SRC="../CLM3.5"
PARFLOW_SRC="../parflow"
COSMO_SRC="../cosmo5.01_fresh"

if if echo "$MODEL_ID" | grep -qE 'COSMO'; then
	cd $COSMO_SRC
	git checkout tsmp-oasis
        cd ../bldsva
if if echo "$MODEL_ID" | grep -qE 'ParFlow'; then
        cd $PARFLOW_SRC
        git checkout v3.12.0 
        cd ../bldsva


# Where build artifacts and binaries will be saved
BUILD_DIR="../bld/${SYSTEMNAME^^}_${MODEL_ID}"
INSTALL_DIR="..//bin/${SYSTEMNAME^^}_${MODEL_ID}"
BUILD_LOG="$(dirname ${BUILD_DIR})/${MODEL_ID}_$(date +%Y.%m.%d_%H.%M).log"

rm -rf ${BUILD_DIR}
mkdir -p ${BUILD_DIR} ${INSTALL_DIR}
#cd ..
cmake -S . -B ${BUILD_DIR}                  \
      -DCLM35_SRC=${CLM35_SRC}              \
      -DCOSMO_SRC=${COSMO_SRC}              \
      -DPARFLOW_SRC=${PARFLOW_SRC}          \
      -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      |& tee ${BUILD_LOG}

# Build and install
cmake --build ${BUILD_DIR} |& tee -a ${BUILD_LOG}
cmake --install ${BUILD_DIR} |& tee -a ${BUILD_LOG}

echo ""
echo "Successfully installed ${MODEL_ID} to ${INSTALL_DIR}"
