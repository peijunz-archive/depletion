#! /usr/bin/env python3
from pylab import *
from segm import *
from curve import *
#Test function
plt.axis('equal')
hold(True)
#测试旋转平移等功能的正确性
c1=conc(1,0.1,2,0.2)
c2=shiftc(conc(1,0.1,2,-0.2),[1.2,0.5])
ce1=extent(c1,.05)
ce2=extent(c2,.05)
#s=segm(0,[0,0],[-0.1,1])
#crm=rminsect(ce)
drawc(c1,'k')
drawc(c2,'k')
drawc(ce1,'g')
drawc(ce2,'r')
print(intsecc(ce1,ce2))
savefig('verge.pdf')
print('Tasks completed!')
