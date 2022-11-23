[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![build](https://github.com/meelgroup/cmsgen/actions/workflows/build.yml/badge.svg)](https://github.com/meelgroup/cmsgen/actions/workflows/build.yml)

CMSGen a fast uniform-like sampler
===========================================

This system provides CMSGen, a fast uniform-like sampler. While we give no guarntees that the sampling is uniform, it is currently the best non-guaranteed uniform sampler as per our testing with [Barbarik](https://github.com/meelgroup/barbarik). In case you need guaranteed uniform sampling, please check out [UniGen](https://github.com/meelgroup/unigen). When citing CMSGen, always reference our [FMCAD'21 paper](https://meelgroup.github.io/files/publications/fmcad21_shakuni.pdf) (bibtex [here](https://meelgroup.github.io/publication/fmcad21/cite.bib)).

Command-line usage
-----
Let's take a DIMACS CNF file `input.cnf`. To get 50 uniform-like samples, run:

```
./cmsgen input.cnf --samplefile mysamples.out --samples 50
 Writing samples to file: mysamples.out
```

You can add weights for polarities like this:
```
p cnf 2 1
w 1 0.1
2 0
```

This will give solutions where variable 2 is TRUE and where variable 1 is TRUE with a probability of 0.1. This is indeed the case:

```
$ echo "p cnf 2 1
w 1 0.1
2 0" | ./cmsgen
c Writing samples to file: samples.out
c Finished generating all 100 samples
c Total time: 0.0016 s
$ sort -n samples.out | uniq -c
     92 -1 2 0
      8 1 2 0
```

In other words, we got 8% samples where we had variable 1 as TRUE.


Compiling in Linux
-----

To build and install, issue:

```
sudo apt-get install build-essential cmake
# not required but very useful
sudo apt-get install zlib1g-dev libboost-program-options-dev help2man
git clone https://github.com/msoos/cryptominisat
cd cryptominisat
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
tar xzvf cmsgen-[version].tar.gz
cd cmsgen-[version]
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
C:\> unzip cmsgen-[version].zip
C:\> rename cmsgen-[version] cmsgen
C:\> cd cmsgen
C:\cms> mkdir build
C:\cms> cd build

C:\cms\build> [ download http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.zip ]
C:\cms\build> unzip boost_1_59_0.zip
C:\cms\build> mkdir boost_1_59_0_install
C:\cms\build> cd boost_1_59_0
C:\cms\build\boost_1_59_0> bootstrap.bat --with-libraries=program_options
C:\cms\build\boost_1_59_0> b2 --with-program_options address-model=64 toolset=msvc-14.0 variant=release link=static threading=multi runtime-link=static install --prefix="C:\cms\build\boost_1_59_0_install" > boost_install.out
C:\cms\build\boost_1_59_0> cd ..

C:\cmsgen\build> git clone https://github.com/madler/zlib
C:\cmsgen\build> cd zlib
C:\cmsgen\build\zlib> git checkout v1.2.8
C:\cmsgen\build\zlib> mkdir build
C:\cmsgen\build\zlib> mkdir myinstall
C:\cmsgen\build\zlib> cd build
C:\cmsgen\build\zlib\build> cmake -G "Visual Studio 14 2015 Win64" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=C:\cms\build\zlib\myinstall ..
C:\cmsgen\build\zlib\build> msbuild /t:Build /p:Configuration=Release /p:Platform="x64" zlib.sln
C:\cmsgen\build\zlib\build> msbuild INSTALL.vcxproj
C:\cmsgen\build> cd ..\..

C:\cmsgen\build> cmake -G "Visual Studio 14 2015 Win64" -DCMAKE_BUILD_TYPE=Release -DSTATICCOMPILE=ON -DZLIB_ROOT=C:\cms\build\zlib\myinstall -DBOOST_ROOT=C:\cms\build\boost_1_59_0_install ..
C:\cmsgen\build> cmake --build --config Release .
```

You now have the static binary under `C:\cmsgen\build\Release\cmsgen.exe`


CryptoMiniSat vs CMSGen
-----
CMSGen is a version of CryptoMiniSat that's made simpler to work with for researchers. But it is possible to get a version of CryptoMiniSat, build it, and run it with a specific command line set to achieve the _same exact behaviour_ as per the paper published. To do this:

* Clone this repository
* Execute: `git checkout 028357ee4b29a1da36e4d0929cc44138e5ae36e9`
* Build as per above
* Execute: `./cryptominisat5 --maxsol $1  --nobansol --restart fixed --maple 0 --verb 0 --scc 1 -n 1  --presimp 0 --polar rnd --freq 0.9999 --fixedconfl  $2 --random $3 --dumpresult $4 [CNFFILE]`, where `--random` is the seed and `--maxsol` is the number of samples.

