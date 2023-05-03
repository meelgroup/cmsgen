#!/bin/python3

import pycmsgen
import sys

fname = sys.argv[1]
print("using file: ", fname)
solver = pycmsgen.Solver()

num_cls = 0
with open(fname, "r") as f:
    for line in f:
        line = line.strip()
        if len(line) == 0: continue
        if line[0] == "c": continue
        if line[0] == "p": continue
        # print(line)
        clause = []
        for l in line.split(" "):
            # print(l)
            l = int(l)
            if l == 0: break
            clause.append(l)
        solver.add_clause(clause)
        num_cls+=1
print("num cls: ", num_cls)
print("num vars: ", solver.nb_vars())

for _ in range(10):
    ret = solver.solve()
    print("ret: ", ret)
    m = solver.get_model()
