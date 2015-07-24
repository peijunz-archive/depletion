#! /usr/bin/env python3
from pylab import *
from segm import *
def tricur(r):
    p0=rect(r,-uni)
    p1=rect(r,3*uni)
    s1=segm(0,p0,p1)
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
def tritri(r,d):
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
def drawc(cur,prop='b'):
    '''draw the curve in matplotlib'''
    for seg in cur:
        seg.dra(prop)
def shiftc(cur,delta):
    '''shift of a curve'''
    p=[]
    for seg in cur:
        p.append(seg.shifts(delta))
    return p
#@profile
def rotc(cur,t):
    '''rotation of a curve'''
    p=[]
    for seg in cur:
        p.append(seg.rots(t))
    return p
#@profile
def areal(poi):
    '''area of a polygon consists of line segments'''
    l=len(poi)
    area=0
    if l<2:
        return 0
    for i in range(l):
        area+=cross2(poi[i],poi[(i+1)%l])/2
    return area
#@profile
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
            area+=sign(cross2(p1,p2-p1))*seg.r**2*(dt-sin(dt))/2
            poly.append(seg.p+p1)
    #print('area of poly:', ap)
    area+=areal(poly)
    return area
#@profile
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
#@profile
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
#@profile
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
#@profile
def extent_raw(cur,d):
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
        sig=sign(cross2(former.enddir()[1], latter.enddir()[0]))
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
    return ext
def extent(cur,d):
    return max(sepcur(cutpcur(cur2pcur(extent_raw(cur,d)))),key=areac)
#@profile
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
    return sum(pos)-max(pos)
#@profile
def issec(cur1,cur2):
    for si in cur1:
        for sj in cur2:
            if len(intsecs(si,sj))!=0:
                return True
    return False
#@profile
def poten(cur1,cur2,distance,rp):
    c2=shiftc(cur2,[distance,0])
    #if issec(cur1,c2):
        #return 0.1
    cr1=extent(cur1,rp)
    cr2=extent(c2,rp)
    return -intsecc(cr1,cr2)
#@profile
def closest(t1,t2,rp):
    #TODO:提高效率
    lef=1
    rig=3
    num=12+int(floor(log((rig-lef)/rp)))
    for i in range(num):
        cen=(lef+rig)/2
        #print(issec(t1,shiftc(t2,[cen,0])))
        if issec(t1,shiftc(t2,array([cen,0]))):
            lef=cen
        else:
            rig=cen
    return rig
def draw2_raw(c1,cur2,distance,rp):
    c2=shiftc(cur2,[distance,0])
    drawc(c1,'k')
    drawc(c2,'k')
    t1=extent_raw(c1,rp)
    t2=extent_raw(c2,rp)
    drawc(t1,'g')
    drawc(t2,'r')
    print(intsecc(t1,t2))
    return 0
def draw2(c1,cur2,distance,rp):
    c2=shiftc(cur2,[distance,0])
    drawc(c1,'k')
    drawc(c2,'k')
    t1=extent(c1,rp)
    t2=extent(c2,rp)
    drawc(t1,'g')
    drawc(t2,'r')
    print(intsecc(t1,t2))
    return 0
def binde(a,b,t1,t2,rp):

    c1=rotc(a,t1)
    c2=rotc(b,t2)
    r1=closest(c1,c2,rp)
    return poten(c1,c2,r1,rp)
def center1(delta,pt):
    x=linspace(-uni-delta*uni,-uni+delta*uni,pt)
    y=linspace(uni-delta*uni,uni+delta*uni,pt)
    return x,y
def center2(delta,pt):
    x=linspace(uni-delta*uni,uni+delta*uni,pt)
    y=linspace(-uni-delta*uni,-uni+delta*uni,pt)
    return x,y
def curve1():
    c1=conc(1,0.2,1,0.5)
    c2=conc(1,0.2,1,-0.5)
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
