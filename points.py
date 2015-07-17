#! /usr/bin/env python3
from pylab import *
infs=1e-10
def cs(v1,v2):
    return sign(cross(v1,v2))
def atanv(v):
    return arctan2(v[1],v[0])
def rect(r,theta):
    return array([r*cos(theta),r*sin(theta)])
def rotp(p,theta):
    '''rotation a point'''
    c=cos(theta)
    s=sin(theta)
    M=array([[c,-s],[s,c]])
    return M.dot(p)
def modi(base,theta):
    """Get a theta in the (0,2*pi) region by +2*k*pi"""
    return theta-floor((theta-base)/2/pi)*2*pi
def betw(q,theta):
    t1=modi(min(q),theta+infs)
    t2=modi(min(q),theta-infs)
    if t1<max(q) and t2<max(q):
        return (t1+t2)/2
    else:
        return False
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
    #if dot(v1,v2) < norm(v1)*norm(v2)*cos(t0) and cross(v1,v2)<0:
        #return True
    #else:
        #return False
