#! /usr/bin/env python3
from curve import *
plt.axis('equal')
hold(True)
rp=0.03
#测试旋转平移等功能的正确性
c1,c2=curve3()
t1=rand()*4*uni
t2=rand()*4*uni
ce1=rotc(c1,-2.5*uni)
ce2=rotc(c2,uni)
r1=closest(ce1,ce2,rp)
#drawc(ce2)
#for i in ce2:
    #i.pri()
draw2(ce1,ce2,r1,rp)
savefig('shape.pdf')
print('Tasks completed!')
