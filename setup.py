#
# CryptoMiniSat
#
# Copyright (c) 2009-2017, Mate Soos. All rights reserved.
# Copyright (c) 2017, Pierre Vignet
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


from setuptools import Extension, setup
import toml
import pathlib

def _parse_toml(pyproject_path):
    pyproject_text = pyproject_path.read_text()
    pyproject_data = toml.loads(pyproject_text)
    return pyproject_data['project']['version']

def gen_modules(version):
    modules = Extension(
        name = "pycmsgen",
        include_dirs = ["src/"],
        sources = ["python/src/pycmsgen.cpp",
                   "python/src/GitSHA1.cpp",
                   "src/clauseallocator.cpp",
                   "src/clausecleaner.cpp",
                   "src/clauseusagestats.cpp",
                   "src/cnf.cpp",
                   "src/completedetachreattacher.cpp",
                   "src/cmsgen.cpp",
                   "src/distillerlong.cpp",
                   "src/distillerlongwithimpl.cpp",
                   "src/gaussian.cpp",
                   "src/hyperengine.cpp",
                   "src/implcache.cpp",
                   "src/intree.cpp",
                   "src/matrixfinder.cpp",
                   "src/occsimplifier.cpp",
                   "src/packedrow.cpp",
                   "src/prober.cpp",
                   "src/propengine.cpp",
                   "src/reducedb.cpp",
                   "src/sccfinder.cpp",
                   "src/searcher.cpp",
                   "src/searchstats.cpp",
                   "src/solutionextender.cpp",
                   "src/solverconf.cpp",
                   "src/solver.cpp",
                   "src/stamp.cpp",
                   "src/str_impl_w_impl_stamp.cpp",
                   "src/subsumeimplicit.cpp",
                   "src/subsumestrengthen.cpp",
                   "src/varreplacer.cpp",
                   "src/xorfinder.cpp"
               ],
        extra_compile_args = ['-I../', '-Isrc/', '-std=c++17'],
        define_macros=[("USE_GAUSS", "1"), ("TRACE", ""), ("CMSGEN_FULL_VERSION", "\""+version+"\"")],
        language = "c++",
    )
    return modules

if __name__ == '__main__':
    pyproject_path = pathlib.Path('pyproject.toml')
    version = _parse_toml(pyproject_path)
    modules = gen_modules(version)
    setup(
        ext_modules =  [modules],
    )
