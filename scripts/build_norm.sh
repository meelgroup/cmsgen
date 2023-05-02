#!/bin/bash

set -e

rm -rf src-cmsgen cm* CM* lib* cryptomini* Testing* tests* pycryptosat include tests cusp* scalmc* utils drat-trim *.cmake Make*
cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make -j26
