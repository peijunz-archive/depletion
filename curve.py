#! /usr/bin/env python3
from pylab import *
from segm import *
uni=pi/6
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
def verge(cur):
    '''print the points'''
    for seg in cur:
        seg.pri()
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
def cur2pcur(cur):
    '''Only applicable when cur is close'''
    pcur=[]
    l=len(cur)
    for i in range(l):
        seg=cur[i]
        start=seg.endp()[0]
        pcur.append([start,seg,(i+1)%l])
    return pcur
def addpcur(pc1,pc2):
    l=len(pc1)
    ap=pc1+pc2
    for p in ap[l:]:
        p[2]+=l
    return ap
def cutpcur(pc):
    i=1
    while i<len(pc):
        for j in range(i):
            lis=intsecs(pc[j][1],pc[i][1])
            if len(lis)==1:
                p=lis[0]
                [s1,s2]=pc[j][1].sep(p)
                [s3,s4]=pc[i][1].sep(p)
                jnex=pc[j][2]
                inex=pc[i][2]
                pc[i][1:]=[s3,len(pc)]
                pc.append([p,s2,jnex])
                pc[j][1:]=[s1,len(pc)]
                pc.append([p,s4,inex])

            elif len(lis)>1:
                print('Consider more! FIX ME')
        i+=1
    return pc
def sepcur(pc):
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
    print('area of bow:', area)
    ap=areal(poly)
    print('area of poly:', ap)
    area+=areal(poly)
    return area
def extent(cur,d):
    '''cur is the curve to extent and d the distance of the pen'''
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
    return ext
def jointc(cn,seg):
    '''multi points condition not concerned'''
    tmp=[]
    if len(cn)<=1:
        cn.append(seg)
        return cn
    for s in cn:
        cro=intsecs(s,seg)
        if cro:
            i=cn.index(s)
            close=cn[i:]
            cn=cn[0:i]
            po=cro[0]
            s1=s.sep(po)[0]
            s2=seg.sep(po)[1]
            s1.pri()
            s2.pri()
            cn.append(s1)
            cn.append(s2)
            print('newnew')
            verge(cn)
            return cn
    cn.append(seg)
    return cn
def rminsect(cur):
    '''remove the inner intersect segments of a curve. DIRTY'''
    #cgroup=[]
    cn=[]
    l=len(cur)
    for i in range(l):
        seg=cur[i]
        #print('CN1:')
        #verge(cn)
        #print('seg:')
        seg.pri()
        cn=jointc(cn,seg)
        #tmp=jointc(cn,seg)
        #cn=tmp[0]
        #cgroup.append(tmp[1])
        print('CN2:')
        verge(cn)
    #cgroup.append(cn)
    #carea=[areac(j) for j in cgroup]
    return cn
if __name__=='__main__':
    #Test function
    plt.axis('equal')
    hold(True)
    #测试旋转平移等功能的正确性
    c=conc(1,0.1,2,0.2)
    ce=extent(c,.3)
    s=segm(0,[0,0],[-0.1,1])
    #crm=rminsect(ce)
    drawc(ce,'b')
    drawc(c,'r')
    t=cutpcur(cur2pcur(ce))
    sepc=sepcur(t)
    ar=[areac(cu) for cu in sepc]
    i=ar.index(max(ar))
    fcur=sepc[i]
    drawc(fcur,'g')
    #print(len(crm))
    #verge(crm)
    #drawc(crm,'g')
    #d=conc(1.1,0.11,2.2,+0.3)
    #drawc(c,'r')
    #drawc(shiftc(rotc(d,2*uni),[3,4]),'b')
    #测试判断相交性的函数
    #intri([0,1],[1,0],[0,0],[0.2,0.3])
    #p1=array([0,1])
    #p2=array([2,0])
    #q1=array([0.3,0.3])
    #q2=array([0.7,0.7])
    #print(interll(p1,p2,q1,q2))
    #q=array([-20*uni,-3*uni])
    #print(modi(q))
    #测试面积函数的正确性
    #l=segm(0,[-0.6,-0.81],[1,0.2])
    #c1=segm(1,[-0.3,3],[2*uni,11*uni])
    #c2=segm(1.2,[0,0.1],[3*uni,9*uni])
    #print(intercc(c1,c2))
    #c1.dra('b')
    #c2.dra('r')
    savefig('verge.pdf')
    ##print(interlc(l,c1))
    #cur=[segm(0,[0,0],[sqrt(3)/2,0.5]),segm(1,[0,0],[uni,2*uni]),segm(0,[0.5,sqrt(3)/2],[0,0])]
    #print(areac(cur))
    #print(pi/12)
    #测试去除内自交
    print('Tasks completed!')
#Area function should me modified for the clockwise case?
