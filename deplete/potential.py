#!/usr/bin/env python3
from .curve import *
from .instance import *


def binding_energy(*arg):
    pass


c1, c2 = curve7()
rp = 0.03
side = 2
# x,y=center0(side,21,[1.6,-0.4])
x, y = center0(side, 11, [0, 0])
m = len(x)
n = len(y)
z = zeros([n, m])
p = 0
for i in range(n):
    # print(i)
    for j in range(m):
        z[i, j] = binding_energy(c1, c2, y[i], x[j], rp)
        print(x[j] / uni, y[i] / uni, z[i, j])
x /= uni
y /= uni
CS = contourf(x, y, z)
colorbar(CS)
plt.axis('equal')
# drawc(c)
savefig('binding.svg')
savefig('binding.pdf')
