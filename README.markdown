[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![build](https://github.com/meelgroup/cmsgen/actions/workflows/build.yml/badge.svg)](https://github.com/meelgroup/cmsgen/actions/workflows/build.yml)

# CMSGen a fast weighted, uniform-like sampler

This system provides CMSGen, a fast weighted uniform-like sampler. While we give no
guarantees that the sampling follows the desired distribution, it is currently
the best non-guaranteed (uniform) sampler as per our testing with
[Barbarik](https://github.com/meelgroup/barbarik). In case you need guaranteed
uniform sampling, please check out
[UniGen](https://github.com/meelgroup/unigen). When citing CMSGen, always
reference our [FMCAD'21
paper](https://meelgroup.github.io/files/publications/fmcad21_shakuni.pdf)
(bibtex [here](https://meelgroup.github.io/publication/fmcad21/cite.bib)).

## Command-line usage
Let's take a DIMACS CNF file `input.cnf`. To get 50 uniform-like samples, run:

```
./cmsgen input.cnf --samplefile mysamples.out --samples 50 --seed 0
 Writing samples to file: mysamples.out
```

You can add weights for polarities like this:
```
p cnf 2 1
w 1 0.1
2 0
```

The above gives solutions where variable 2 is TRUE always, and where variable 1 is TRUE
with a probability of 0.1. This is indeed the case:

```bash
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

## Python usage

Install via pip:
```bash
pip install pycmsgen
```

Then:
```python
import pycmsgen
solver = pycmsgen.Solver(seed=0)
solver.add_clause([1,2,4])
solver.add_clause([1,2,-5])
sat, sol = solver.solve()
```

Where the return value `sat` will be `True`, indicating there is a solution found (i.e. it's not unsatisfiable), and `sol[1]`, `sol[2]`, etc. will indicate the solution to variables 1, 2, etc.

## Compiling in Linux

To build and install, issue:

```bash
sudo apt-get install build-essential cmake
# not required but very useful
sudo apt-get install zlib1g-dev help2man
git clone https://github.com/meelgroup/cmsgen
cd cmsgen
mkdir build && cd build
cmake ..
make
sudo make install
sudo ldconfig
```

## Compiling in Mac OSX

First, you must get Homebew from https://brew.sh/ then:

```bash
brew install cmake zlib
tar xzvf cmsgen-[version].tar.gz
cd cmsgen-[version]
mkdir build && cd build
cmake ..
make
sudo make install
```
