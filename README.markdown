[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![build](https://github.com/meelgroup/cmsgen/actions/workflows/build.yml/badge.svg)](https://github.com/meelgroup/cmsgen/actions/workflows/build.yml)

CMSGen uniform-like sampler
===========================================

This system provides CMSGen, a fast uniform-like sampler.

When citing, always reference our [FMCAD'21 paper](https://TODO), bibtex record is [here](https://TODO).

Command-line usage
-----
Let's take a DIMACS CNF file `input.cnf`. To get 50 uniform-like samples, run:

```
./cmsgen input.cnf --samplesfile mysamples.out --samples 50
 Writing samples to file: mysamples.out
```

Compiling in Linux
-----

To build and install, issue:

```
sudo apt-get install build-essential cmake
# not required but very useful
sudo apt-get install zlib1g-dev libboost-program-options-dev help2man
tar xzvf cryptominisat-version.tar.gz
cd cryptominisat-version
mkdir build && cd build
cmake ..
make
sudo make install
sudo ldconfig
```

Compiling in Mac OSX
-----

First, you must get Homebew from https://brew.sh/ then:

```
brew install cmake boost zlib
tar xzvf cryptominisat-version.tar.gz
cd cryptominisat-version
mkdir build && cd build
cmake ..
make
sudo make install
```

Compiling in Windows
-----

You will need python installed, then for Visual Studio 2015:

```
C:\> [ download cryptominisat-version.zip ]
C:\> unzip cryptominisat-version.zip
C:\> rename cryptominisat-version cms
C:\> cd cms
C:\cms> mkdir build
C:\cms> cd build

C:\cms\build> [ download http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.zip ]
C:\cms\build> unzip boost_1_59_0.zip
C:\cms\build> mkdir boost_1_59_0_install
C:\cms\build> cd boost_1_59_0
C:\cms\build\boost_1_59_0> bootstrap.bat --with-libraries=program_options
C:\cms\build\boost_1_59_0> b2 --with-program_options address-model=64 toolset=msvc-14.0 variant=release link=static threading=multi runtime-link=static install --prefix="C:\cms\build\boost_1_59_0_install" > boost_install.out
C:\cms\build\boost_1_59_0> cd ..

C:\cms\build> git clone https://github.com/madler/zlib
C:\cms\build> cd zlib
C:\cms\build\zlib> git checkout v1.2.8
C:\cms\build\zlib> mkdir build
C:\cms\build\zlib> mkdir myinstall
C:\cms\build\zlib> cd build
C:\cms\build\zlib\build> cmake -G "Visual Studio 14 2015 Win64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=C:\cms\build\zlib\myinstall ..
C:\cms\build\zlib\build> msbuild /t:Build /p:Configuration=Release /p:Platform="x64" zlib.sln
C:\cms\build\zlib\build> msbuild INSTALL.vcxproj
C:\cms\build> cd ..\..

C:\cms\build> cmake -G "Visual Studio 14 2015 Win64" -DCMAKE_BUILD_TYPE=Release -DSTATICCOMPILE=ON -DZLIB_ROOT=C:\cms\build\zlib\myinstall -DBOOST_ROOT=C:\cms\build\boost_1_59_0_install ..
C:\cms\build> cmake --build --config Release .
```

You now have the static binary under `C:\cms\build\Release\cryptominisat5.exe`


CryptoMiniSat vs CMSGen
-----
CMSGen is a version of CryptoMiniSat that's made simpler to work with for researchers. But it is possible to get a version of CryptoMiniSat, build it, and run it with a specific command line set to achieve the _same exact behaviour_ as per the paper published. To do this:

* Clone this repository
* Execute: `git checkout 028357ee4b29a1da36e4d0929cc44138e5ae36e9`
* Execute: `./cryptominisat5 --maxsol $1  --nobansol --restart fixed --maple 0 --verb 0 --scc 1 -n 1  --presimp 0 --polar rnd --freq 0.9999 --fixedconfl  $2 --random $3 --dumpresult $4 [CNFFILE]` where `--random` is the seed and `--maxsol` is the number of samples.

