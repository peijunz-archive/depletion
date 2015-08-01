#! /usr/bin/env python3
from curve import *
from instance import *
plt.axis('equal')
hold(True)
rp=0.03
#测试旋转平移等功能的正确性
c1,c2=curve7()
t1=0
t2=0
ce1=rotc(c1,t1*uni)
c1c=extent(ce1,rp)
ce2=shiftc(rotc(c2,t2*uni),[3,0])
c2c=extent(ce2,rp)
plt.axis('off')
drawc(ce1,'k')
drawc(c1c,'r')
drawc(ce2,'k')
drawc(c2c,'g')
#plt.clf()
savefig('shape.svg')
print('Terminated Normally!')
