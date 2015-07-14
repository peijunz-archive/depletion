#! /usr/bin/env python3
from pylab import *
def modi(q):
    """Get a theta in the (0,2*pi) region by +2*k*pi"""
    theta=q[1]
    base=q[0]
    return array([base,theta-floor((theta-base)/2/pi)*2*pi])
class segm:
    '''Line/Circle segment'''
    def __init__(self,r,p,q):
        '''The meaning of parameters:
        r   radius
        For line, p(x,y) is the start point and q(z,w) the end point
        For circle, p(x,y) is the center and q(z,w) the start/end angle'''
        self.r=r
        self.p=array(p)
        if r==0:
            self.q=array(q)
        else:
            self.q=modi(q)
    def pri(self):
        """Print the info of the segmant"""
        if self.r == 0:
            print('Line\t',self.p,'____',self.q)
        else:
            print('Circle\tr=',self.r,'center',self.p,'range',self.q)
    def dra(self,prop):
        [x,y]=self.p
        [z,w]=self.q
        """Draw the segment by matplotlib"""
        if self.r == 0:
            plot((x,z),(y,w),prop)
        else:
            xx=[]
            yy=[]
            for t in arange(z,w,0.01):
                xx.append(x+self.r*cos(t))
                yy.append(y+self.r*sin(t))
            plot(xx,yy,prop)
    def rots(self,t):
        '''rotate the segment'''
        [x,y]=self.p
        [z,w]=self.q
        p1=rotp(self.p,t)
        if self.r==0:
            q1=rotp(self.q,t)
            return segm(0,p1,q1)
        else:
            return segm(self.r,p1,self.q+t)
    def shifts(self,delta):
        '''add a delta(dx,dy) shift to the segment'''
        if self.r==0:
            return segm(0, self.p+delta, self.q+delta)
        else:
            return segm(self.r, self.p+delta, self.q)
        
uni=pi/6
def rotp(p,theta):
    '''rotation a point'''
    c=cos(theta)
    s=sin(theta)
    M=array([[c,-s],[s,c]])
    return M.dot(p)
def conc(H,b,r,h):
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
        R=mo(q-cen)
        cur.append(segm(R,cen,[z,w]))
        cur.append(segm(0,q,top))
        cur.append(segm(0,top,p))
    return cur
def verge(cur):
    for seg in cur:
        seg.pri()
def drawc(cur,prop):
    for seg in cur:
        seg.dra(prop)
def shiftc(cur,delta):
    p=[]
    for seg in cur:
        p.append(seg.shifts(delta))
    return p
def rotc(cur,t):
    p=[]
    for seg in cur:
        p.append(seg.rots(t))
    return p
def intri(a,b,c,o):
    '''利用三个叉乘矢量方向相同来判断'''
    oa=o-a
    ob=o-b
    oc=o-c
    C=cross(oa,ob)
    A=cross(ob,oc)
    B=cross(oc,oa)
    if A*B>0 and B*C>0 and A*C>0:
        return True
    else:
        return False
def inbow(r,cen,ang,o):
    '''Judge whether a point o is in a bow
    利用圆周角大小关系以及矢量叉乘
    '''
    t0=pi-(ang[1]-ang[0])/2
    l=cen+rect(r,z)
    r=cen+rect(r,w)
    v1=o-l
    v2=o-r
    if dot(v1,v2) < mo(v1)*mo(v2)*cos(t0) and cross(v1,v2)<0:
        return True
    else:
        return False
def atanv(v):
    return arctan2(v[1],v[0])
def mo(v):
    return sqrt(sum(array(v)**2))
def rect(r,theta):
    return array([r*cos(theta),r*sin(theta)])
def interll(s1,s2):
    '''lm is lambda and mu'''
    [l,m]=inv(array([s1.p-s1.q,s2.p-s2.q])).T.dot(s2.p-s1.q)
    if 0<l and l<1 and 0<m and m<1:
        return [l*s1.p+(1-l)*s1.q]
    else:
        return array([])
def betw(q,theta):
    return modi([q[0],theta])[1]<q[1]
def interlc(l,c):
    '''line---circle
    '''
    #三角形的三条边
    a=l.p-c.p
    b=l.q-c.p
    delta=l.p-l.q
    direc=delta/mo(delta)#直线的方向矢量
    vert=a-dot(a,direc)*direc#垂线矢量，减去它相当于做投影变换
    h=mo(vert)#垂线长度height
    ah=a-vert
    bh=b-vert
    ret=[]
    if h<abs(c.r):
        con=sqrt(c.r**2-h**2)
        pm=[con*direc,-con*direc]
        for i in pm:
            theta=atanv(vert+i)
            #print(theta,c.q)s
            if dot(ah-i,bh-i)<0 and betw(c.q,theta):
                ret.append(c.p+vert+i)
    return ret
def intercc(c1,c2):
    delta=c2.p-c1.p
    d=mo(delta)
    ret=[]
    if abs(c1.r-c2.r)<d and d<c1.r+c2.r:
        t1=atanv(delta)
        t2=t1+pi
        phi1=arccos((c1.r**2+d**2-c2.r**2)/(2*c1.r*d))
        phi2=-arccos((c2.r**2+d**2-c1.r**2)/(2*c2.r*d))
        print('theta',t1,t2)
        print('phi:',phi1,phi2) 
        for i in [1,-1]:
            if betw(c1.q,t1+i*phi1) and betw(c2.q,t2+i*phi2):
                ret.append(c1.p+rect(c1.r,t1+i*phi1))
    return ret
def intsecs(s1,s2):
    '''Judge and get the intersect point(if exist)'''
    if s1.r==0:
        if s2.r==0:
            return interll(s1,s2)
        else:
            return interlc(s1,s2)
    else:
        if s2.r==0:
            return interlc(s2,s1)
        else:
            return intercc(s1,s2)
        
if __name__=='__main__':
    #Test function
    plt.axis('equal')
    hold(True)
    #c=conc(1,0.1,2,-0.2)
    #d=conc(1.1,0.11,2.2,+0.3)
    #drawc(c,'r')
    #drawc(shiftc(rotc(d,2*uni),[3,4]),'b')
    #intri([0,1],[1,0],[0,0],[0.2,0.3])
    #p1=array([0,1])
    #p2=array([2,0])
    #q1=array([0.3,0.3])
    #q2=array([0.7,0.7])
    #print(interll(p1,p2,q1,q2))
    #q=array([-20*uni,-3*uni])
    #print(modi(q))
    l=segm(0,[-0.6,-0.81],[1,0.2])
    c1=segm(1,[-0.3,3],[2*uni,11*uni])
    c2=segm(1.2,[0,0.1],[3*uni,9*uni])
    print(intercc(c1,c2))
    c1.dra('b')
    c2.dra('r')
    savefig('verge.pdf')
    #print(interlc(l,c1))
    print('Tasks completed!')