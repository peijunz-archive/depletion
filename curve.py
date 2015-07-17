#! /usr/bin/env python3
from pylab import *
from segm import *
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
        R=norm(q-cen)
        cur.append(segm(R,cen,[z,w]))
        cur.append(segm(0,q,top))
        cur.append(segm(0,top,p))
    return cur
def drawc(cur,prop):
    '''draw the curve in matplotlib'''
    for seg in cur:
        seg.dra(prop)
def shiftc(cur,delta):
    '''shift of a curve'''
    p=[]
    for seg in cur:
        p.append(seg.shifts(delta))
    return p
def rotc(cur,t):
    '''rotation of a curve'''
    p=[]
    for seg in cur:
        p.append(seg.rots(t))
    return p
def areal(poi):
    '''area of a polygon consists of line segments'''
    l=len(poi)
    area=0
    if l<2:
        return 0
    for i in range(l):
        area+=cross(poi[i],poi[(i+1)%l])/2
    return area
def areac(cur):
    '''area of a curve consists of line/circle segments'''
    poly=[]
    area=0
    for seg in cur:
        if seg.r==0:
            poly.append(seg.p)
        else:
            dt=seg.q[1]-seg.q[0]
            p1=rect(seg.r,seg.q[0])#TODO, midify it, use the end point function
            p2=rect(seg.r,seg.q[1])
            area+=sign(cross(p1,p2-p1))*seg.r**2*(dt-sin(dt))/2
            poly.append(seg.p+p1)
    #print('area of bow:', area)
    ap=areal(poly)
    #print('area of poly:', ap)
    area+=areal(poly)
    return area
def cur2pcur(cur):
    '''Only applicable when cur is close.
    p means the data structure is based on points
    convert the segment form curve to points form curve'''
    pcur=[]
    l=len(cur)
    for i in range(l):
        seg=cur[i]
        start=seg.endp()[0]
        pcur.append([start,seg,(i+1)%l])
    return pcur
def addpcur(pc1,pc2):
    '''add 2 pcur together, set the pointer properly'''
    l=len(pc1)
    ap=pc1+pc2
    for p in ap[l:]:
        p[2]+=l
    return ap
def cutpcur(pc):
    '''find all joint points and reconnect the pcurve according to it'''
    i=1
    while i<len(pc):
        for j in range(i):
            lis=intsecs(pc[j][1],pc[i][1])
            if len(lis)==0:
               continue
            if len(lis)==1:
                p=lis[0]
            else:
                dis=[pc[j][1].dist(k) for k in lis]
                p=lis[dis.index(min(dis))]
            [s1,s2]=pc[j][1].sep(p)
            [s3,s4]=pc[i][1].sep(p)
            jnex=pc[j][2]
            inex=pc[i][2]
            pc[i][1:]=[s3,len(pc)]
            pc.append([p,s2,jnex])
            pc[j][1:]=[s1,len(pc)]
            pc.append([p,s4,inex])
        i+=1
    return pc
def sepcur(pc):
    '''separate the different close simple curves'''
    cg=[]#Groups of curve
    l=len(pc)
    sig=zeros(l)
    for i in range(l):
        if sig[i]==0:
            t=pc[i]
            sig[i]=1
            tmp=[t[1]]
            while True:
                j=t[2]
                t=pc[j]
                if j==i:
                    break
                sig[j]=1
                tmp.append(t[1])
            cg.append(tmp)
    return cg
def extent(cur,d):
    '''cur is the curve to extent and d the distance of the pen
    first to get piles of smaller simple curves of its raw extension
    then get the final curve that has biggest area!
    '''
    l=len(cur)
    ext=[]
    for i in range(l):
        former=cur[i]
        latter=cur[(i+1)%l]
        ext.append(former.exts(d))
        sig=sign(cross(former.enddir()[1], latter.enddir()[0]))
        if sig==0:
            continue    #when sig<0, the curve break temporarily
        elif sig<0:    #the condition sig<0 can be simplified to return a null
            ff=former.exts(d)
            ll=latter.exts(d)
            s=segm(0,ff.endp()[1],ll.endp()[0])
            #s.dra('b')
            ext.append(s)
        else:
            c=former.endp()[1]
            q=[former.endth()[1],latter.endth()[0]]
            ext.append(segm(d,c,q))
    #Cut the redundant part of the curve
    t=cutpcur(cur2pcur(ext))
    sepc=sepcur(t)
    ar=[areac(cu) for cu in sepc]
    i=ar.index(max(ar))
    fcur=sepc[i]
    return fcur
def intsecc(cur1,cur2):
    pc1=cur2pcur(cur1)
    pc2=cur2pcur(cur2)
    pc=addpcur(pc1,pc2)
    l0=len(pc)
    m=cutpcur(pc)
    if len(m)==l0:
        return 0
    t=sepcur(m)
    ar=array([areac(cu) for cu in t])
    pos=(abs(ar)+ar)/2
    neg=ar-pos
    #print(pos)
    #print(neg)
    #print('area of regions',ar)
    #print('sum of area',sum(ar))
    #print('area of two region')
    #print(areac(cur1))
    #print(areac(cur2))
    #print('Intersection area:',sum(pos)-max(pos))
    return sum(pos)-max(pos)
