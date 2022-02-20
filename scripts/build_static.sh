#!/bin/bash

set -e

rm -rf cm* CM* lib* cryptomini* Testing* tests* pycryptosat include tests
cmake -DSTATICCOMPILE=ON ..
make -j26
strip cmsgen
