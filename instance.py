#!/usr/bin/env python3
from curve import *

def tricur(r):
    p0=rect(r,-uni)
    p1=rect(r,3*uni)
    s1=segm(0,p0,p1)
    s2=s1.rots(4*uni)
    s3=s1.rots(8*uni)
    return [s1,s2,s3]
def tricir(r):
    cen=rect(r,pi+uni)
    p0=rect(1,-uni)
    dr=p0-cen
    t1=atanv(dr)
    t2=2*uni-t1
    R=norm2(dr)
    s1=segm(R,cen,[t1,t2])
    s2=s1.rots(4*uni)
    s3=s1.rots(8*uni)
    return [s1,s2,s3]
def tricave(r,theta):
    #theta should be less than uni
    sc=segm(r,[0,0],[-uni+theta,3*uni-theta])
    p1=rect(r,3*uni-theta)
    p2=array([0,p1[1]-sqrt(3)*p1[0]])
    p3=rect(r,3*uni+theta)
    s1=segm(0,p1,p2)
    s2=segm(0,p2,p3)
    c1=[sc,s1,s2]
    c2=rotc(c1,4*uni)
    c3=rotc(c1,8*uni)
    return c1+c2+c3
def conc(H,b,r,h):
    '''Constrution of basic curve'''
    cur=[]
    for i in range(3):
        rt=i*4*uni
        p=rotp([-b,H],rt)
        q=rotp([b,H],rt)
        top=rotp([0,H+h],rt)
        cen=rect(r,rt+pi+uni)
        p1=rotp(p,-4*uni)
        z=arctan2(p1[1]-cen[1],p1[0]-cen[0])
        w=arctan2(q[1]-cen[1],q[0]-cen[0])
        R=norm2(q-cen)
        cur.append(segm(R,cen,[z,w]))
        cur.append(segm(0,q,top))
        cur.append(segm(0,top,p))
    return cur
def center0(delta,pt,cen=[0,0]):
    cen=array(cen)*uni
    x=linspace(cen[0]-delta*uni,cen[0]+delta*uni,pt)
    y=linspace(cen[1]-delta*uni,cen[1]+delta*uni,pt)
    return x,y
def center1(delta,pt):
    x=linspace(-uni-delta*uni,-uni+delta*uni,pt)
    y=linspace(uni-delta*uni,uni+delta*uni,pt)
    return x,y
def center2(delta,pt):
    x=linspace(uni-delta*uni,uni+delta*uni,pt)
    y=linspace(-uni-delta*uni,-uni+delta*uni,pt)
    return x,y
def curve1():
    c1=conc(1,0.2,0.5,0.4)
    c2=conc(1,0.2,0.5,-0.4)
    return c1,c2
def curve2(r=1,theta=uni/2):
    c1=tricur(r)
    c2=tricave(r,theta)
    return c1,c2
def curve3():
    s1=segm(1,[0,0],[-uni,2*uni])
    s2=segm(-1,[0,sqrt(3)],[-2*uni,-3*uni])
    s3=segm(0,[0,sqrt(3)-1],[0,1])
    c0=[s1,s2,s3]
    c1=rotc(c0,4*uni)
    c2=rotc(c0,8*uni)
    c=c0+c1+c2
    return c,c
def curve4():
    s1=segm(1,[0,0],[0,2*uni])
    s2=segm(0,rect(1,2*uni),rect(1,4*uni))
    c0=[s1,s2]
    c1=rotc(c0,4*uni)
    c2=rotc(c0,8*uni)
    c=c0+c1+c2
    return c,c
def curve5(r=0.8):
    c1=tricir(r)
    c2=tricir(r)
    return c1,c2
def symcave(r,h):
    cen=rect(r,pi+uni)
    p0=rect(1,-uni)
    dr=p0-cen
    t1=atanv(dr)
    t2=2*uni-t1
    R=norm2(dr)
    s1=segm(R,cen,[t1,t2])
    s2=segm(0,[0,1],[0,1-2*h])
    cen1=2*array([0,1-h])-cen
    s3t=segm(-R,cen1,[t2+pi,t1+pi])
    s4t=s1.rots(4*uni)
    s3t.pri()
    s4t.pri()
    p=intsecs(s3t,s4t)[0]
    s3=s3t.sep(p)[0]
    s4=s4t.sep(p)[1]
    c1=[s2,s3,s4]
    return c1+rotc(c1,4*uni)+rotc(c1,8*uni)
def curve6(r=1,h=0.15):
    c1=symcave(r,h)
    return [c1,c1]
def trify(c):
    return c+rotc(c,4*uni)+rotc(c,8*uni)
def curve7(dr=.7):
    cen=rect(1,pi+uni)+array([-dr,0])
    p0=[0,sqrt(9/4+sqrt(3)*dr)-0.5]
    R=sqrt(3)+dr
    t1=arccos((sqrt(3)/2+dr)/(sqrt(3)+dr))
    s1=segm(R,cen,[0,t1])
    s2=segm(0,p0,[0,1])
    c1=[s1,s2]
    cen2=cen+rect(sqrt(3)-dr,4*uni)
    s3=segm(R,cen2,[2*uni-t1,2*uni])
    p1=intsecs(s1,s3)[0]
    s11=s1.sep(p1)[1]
    s31=s3.sep(p1)[1]
    c2=shiftc([s2.rev(),s11.rev(),s31],rect(sqrt(3),-2*uni))
    return [trify(c1),trify(c2)]
