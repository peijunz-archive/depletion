#! /usr/bin/env python3
from pylab import *
def modi(theta,base):
    """Get a theta in the (0,2*pi) region by +2*k*pi"""
    return theta-floor((theta-base)/2/pi)*2*pi
class segm:
    '''Line/Circle segment'''
    def __init__(self,r,x,y,z,w):
        '''The meaning of parameters:
        r   radius
        For line, (x,y) is the start point and (z,w) the end point
        For circle, (x,y) is the center and (z,w) the start/end angle
        '''
        self.r=r
        self.x=x
        self.y=y
        self.z=z
        if r==0:
            self.w=w
        else:
            self.w=modi(w,z)
    def pri(self):
        """Print the info of the segmant"""
        if self.r == 0:
            print('Line\t ({0:.3f},{1:.3f})--({2:.3f},{3:.3f})'.format(self.x,self.y,self.z,self.w))
        else:
            print('Circle\tr={0:.3f}, O({1:.3f},{2:.3f}), {3:.3f}--{4:.3f}'.format(self.r,self.x,self.y,self.z,self.w))
    def dra(self,prop):
        """Draw the segment by matplotlib"""
        if self.r == 0:
            plot((self.x,self.z),(self.y,self.w),prop)
        else:
            xx=[]
            yy=[]
            for t in arange(self.z,self.w,0.01):
                xx.append(self.x+self.r*cos(t))
                yy.append(self.y+self.r*sin(t))
            plot(xx,yy,prop)
    def rots(self,t):
        '''rotate the segment'''
        p=rotp(self.x,self.y,t)
        if self.r==0:
            q=rotp(self.z,self.w,t)
            return segm(0,p[0],p[1],q[0],q[1])
        else:
            return segm(self.r,p[0],p[1],self.z+t,self.w+t)
    def shifts(self,dx,dy):
        '''add a (dx,dy) shift to the segment'''
        if self.r==0:
            return segm(0, self.x+dx, self.y+dy, self.z+dx, self.w+dy)
        else:
            return segm(self.r, self.x+dx, self.y+dy, self.z, self.w)
        
uni=pi/6
def rotp(x,y,theta):
    '''rotation a point'''
    c=cos(theta)
    s=sin(theta)
#    M=array([[1,2],[3,4]])
    x1=x*c-y*s
    y1=x*s+y*c
    return array([x1,y1])
def conc(H,b,r,h):
    cur=[]
    for i in range(3):
        rt=i*4*uni
        p=rotp(-b,H,rt)
        q=rotp(b,H,rt)
        top=rotp(0,H+h,rt)
        cen=rotp(r,0,rt+pi+uni)
        p1=rotp(p[0],p[1],-4*uni)
        z=arctan2(p1[1]-cen[1],p1[0]-cen[0])
        w=arctan2(q[1]-cen[1],q[0]-cen[0])
        R=hypot(q[0]-cen[0],q[1]-cen[1])
        cur.append(segm(R,cen[0],cen[1],z,w))
        cur.append(segm(0,q[0],q[1],top[0],top[1]))
        cur.append(segm(0,top[0],top[1],p[0],p[1]))
    return cur
def verge(cur):
    for seg in cur:
        seg.pri()
def drawc(cur,prop):
    for seg in cur:
        seg.dra(prop)
def shiftc(cur,dx,dy):
    p=[]
    for seg in cur:
        p.append(seg.shifts(dx,dy))
    return p
def rotc(cur,t):
    p=[]
    for seg in cur:
        p.append(seg.rots(t))
    return p
def intri(xa,ya,xb,yb,xc,yc,x,y):
    a=[x-xa,y-ya]
    b=[x-xb,y-yb]
    c=[x-xc,y-yc]
    C=cross(a,b)
    A=cross(b,c)
    B=cross(c,a)
    if A*B>0 and B*C>0 and A*C>0:
        return True
    else:
        return False
def inbow(r,cx,cy,z,w,x,y):
    t0=pi+(z-w)/2
    x1=cx+r*cos(z)
    y1=cy+r*sin(z)
    x2=cx+r*cos(w)
    y2=cy+r*sin(w)
    v1=[x-x1,y-y1]
    v2=[x-x2,y-y2]
    if dot(v1,v2) < hypot(*v1)*hypot(*v2)*cos(t0) and cross(v1,v2)<0:
        return True
    else:
        return False
def interll(p1,p2,q1,q2):
    '''lm is lambda and mu'''
    #print(array([p1-p2,q1-q2]))
    #print(inv(array([p1-p2,q1-q2])))
    lm=inv(array([p1-p2,q1-q2])).T.dot(q1-p2)
    #print(lm)
    #print(lm[0]*p1+(1-lm[0])*p2)
    if 0<lm[0] and lm[0]<1 and 0<lm[1] and lm[1]<1:
        return lm[0]*p1+(1-lm[0])*p2
    else:
        return array([])
def interlc()
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
    #c=conc(1,0.1,2,-0.2)
    #d=conc(1.1,0.11,2.2,+0.2)
    #plt.axis('equal')
    #hold(True)
    #drawc(c,'r')
    #drawc(shiftc(rotc(d,2*uni),3,4),'b')
    #savefig('verge.pdf')
    #intri(0,1,1,0,0,0,0.2,0.3)
    p1=array([0,1])
    p2=array([2,0])
    q1=array([0.3,0.3])
    q2=array([0.7,0.7])
    print(interll(p1,p2,q1,q2))
    print('Tasks completed!')
