#! /usr/bin/env python3
from pylab import *
#from numpy.linalg import norm
uni=pi/6
pi2=2*pi
infs=1e-10
#@profile
def atanv(v):
    return arctan2(v[1],v[0])
#@profile
def rect(r,theta):
    return array([r*cos(theta),r*sin(theta)])
def rotp(p,theta):
    '''rotation a point'''
    c=cos(theta)
    s=sin(theta)
    M=array([[c,-s],[s,c]])
    return M.dot(p)
#@profile
def modi(base,theta):
    """Get a theta in the (0,2*pi) region by +2*k*pi
    ??inf
    """
    return theta-floor((theta-base-infs)/pi2)*pi2
#@profile
def betw(q,theta):
    if q[0]<q[1]:
        mi,ma=q
    else:
        ma,mi=q
    t1=modi(mi,theta+2*infs)
    if t1<ma:
        t2=modi(mi,theta)
        if(t2<ma):
            return (t1+t2)/2
    return False
#def cs(v1,v2):
    #return sign(cross(v1,v2))
#def intri(a,b,c,o):
    #'''利用三个叉乘矢量方向相同来判断'''
    #oa=o-a
    #ob=o-b
    #oc=o-c
    #C=cross(oa,ob)
    #A=cross(ob,oc)
    #B=cross(oc,oa)
    #if A*B>0 and B*C>0 and A*C>0:
        #return True
    #else:
        #return False
#def inbow(r,cen,ang,o):
    #'''Judge whether a point o is in a bow
    #利用圆周角大小关系以及矢量叉乘'''
    #t0=pi-(ang[1]-ang[0])/2
    #l=cen+rect(r,z)
    #r=cen+rect(r,w)
    #v1=o-l
    #v2=o-r
    #if dot(v1,v2) < norm2(v1)*norm2(v2)*cos(t0) and cross(v1,v2)<0:
        #return True
    #else:
        #return False
