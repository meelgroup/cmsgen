# pycmsgen uniform-like sampler

This directory provides Python bindings to the CMSGen unform-like sampler,

## Installing

```
pip install pycmsgen
```

## Compiling
If you don't want to use the pip package, you can compile it as:

```
apt-get install python-dev
python -m build
```

## Usage

The `pycmsgen` module has one object, `Solver` that has three main functions
`solve`, `add_clause`, and `get_model`.

The function `add_clause()` takes an iterable list of literals such as
`[1, 2]` which represents the truth `1 or 2 = True`. For example,
`add_clause([1])` sets variable `1` to `True`.

The function `solve()` solves the system of equations that have been added
with `add_clause()`:

```
>>> from pycmsgen import Solver
>>> s = Solver()
>>> s.add_clause([1, 2])
>>> sat = s.solve()
True
>>> print(s.get_model())
[1, 2]
```

The `solve()` method optionally takes an argument `assumptions` that
allows the user to set values to specific variables in the solver in a temporary
fashion. This means that in case the problem is satisfiable but e.g it's
unsatisfiable if variable 2 is FALSE, then `solve([-2])` will return
UNSAT. However, a subsequent call to `solve()` will still return a solution.
If instead of an assumption `add_clause()` would have been used, subsequent
`solve()` calls would have returned unsatisfiable.

`Solver` takes the following keyword arguments:
  * `time_limit`: the time limit (integer)
  * `confl_limit`: the propagation limit (integer)
  * `verbose`: the verbosity level (integer)

Both `time_limit` and `confl_limit` set a budget to the solver. The former is based on time elapsed while the former is based on number of conflicts met during search. If the solver runs out of budget, it returns with `(None, None)`. If both limits are used, the solver will terminate whenever one of the limits are hit (whichever first). Warning: Results from `time_limit` may differ from run to run, depending on compute load, etc. Use `confl_limit` for more reproducible runs.
