#!/usr/bin/env python3
from curve import *
c1,c2=curve1()
rp=0.03
side=1
x,y=center1(side,21)
m=len(x)
n=len(y)
z=zeros([n,m])
p=0
for i in range(n):
    #print(i)
    for j in range(m):
        z[i,j]=binde(c1,c2,y[i],x[j],rp)
        print(x[j]/uni,y[i]/uni,z[i,j])
x/=uni
y/=uni
CS=contourf(x,y,z)
colorbar(CS)
plt.axis('equal')
#drawc(c)
savefig('binding.pdf')
