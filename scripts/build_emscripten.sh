#!/bin/bash

set -e
rm -rf cm* CM* lib* cryptomini* Testing* tests* pycryptosat include tests utils Make*
emcmake cmake -DCMAKE_INSTALL_PREFIX=$EMINSTALL -DENABLE_TESTING=OFF -DPOLYS=OFF ..
emmake make -j26
emmake make install
cp cmsgen.wasm ../html
cp $EMINSTALL/bin/cmsgen.js ../html
