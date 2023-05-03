#!/bin/python3
import pycmsgen

a = pycmsgen.Solver()
a.add_clause([1,2,3,4])
a.set_var_weight(1, 0.05)
a.set_var_weight(2, 0.05)
a.set_var_weight(3, 0.05)

num = [0]*4
for at in range(100):
    x = a.solve()
    m = a.get_model()
    if at==0: print("Example model: ", m)
    for i in range(1, 5) :
        if m[i-1] == i:
            num[i-1]+=1

print("num is: ", num)
