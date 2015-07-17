#! /usr/bin/env python3
from pylab import *
from points import *
class segm:
    '''Line/Circle segment'''
    def __init__(self,r,p,q):
        '''The meaning of parameters:
        r   radius
        For line, p(x,y) is the start point and q(z,w) the end point
        For circle, p(x,y) is the center and q(z,w) the start/end angle'''
        self.p=array(p)
        if r==0:
            self.q=array(q)
        elif r>0:
            self.q=array([q[0],modi(*q)])
        else:
            self.q=array([q[0],modi(*q)-2*pi])
        self.r=abs(r)
    def pri(self):
        """Print the info of the segmant"""
        print(self.r,self.p,self.q)
    def __str__(self):
        '''print end points infomations'''
        return str(list(self.endp()[0]))
    def dra(self,prop):
        [x,y]=self.p
        [z,w]=self.q
        """Draw the segment by matplotlib"""
        if self.r == 0:
            plot((x,z),(y,w),prop)
        else:
            xx=[]
            yy=[]
            for t in linspace(z,w,100):
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
    def exts(self,dis):
        '''Get extension of a SINGLE segment'''
        if self.r==0:
            direc=self.q-self.p
            n=direc/norm(direc)
            return self.shifts(dis*rotp(n,-pi/2))
        if self.r>0:
            if(self.q[1]>self.q[0]):
                return segm(self.r+dis,self.p,self.q)
            else:
                if self.r<dis:
                    print('ERROR: r should be bigger than dis')
                    return self
                return segm(self.r-dis,self.p,self.q)
    def enddir(self):
        '''Calculate Tangent direction vector of start point and end point'''
        if self.r==0:
            delta=self.q-self.p
            t=delta/norm(delta)
            return [t,t]
        else:
            t1=rect(1,self.q[0]+pi/2)
            t2=rect(1,self.q[1]+pi/2)
            sig=sign(self.q[1]-self.q[0])
            return [sig*t1,sig*t2]
    def endth(self):
        '''Calculate the theta of end points(anti-clockwise)
        !!!!not efficient, Need improvement
        '''
        td=self.enddir()
        return array([atanv(td[0])-pi/2,atanv(td[1])-pi/2])
    def endp(self):
        '''calculate the end points'''
        if self.r==0:
            return [self.p,self.q]
        else:
            return [self.p+rect(self.r,self.q[0]),self.p+rect(self.r,self.q[1])]
    def sep(self,pt):
        '''use the node pt to separate the self into two segs'''
        if(self.r==0):
                return [segm(0,self.p,pt),segm(0,pt,self.q)]
        else:
            phi=betw(self.q,atanv(pt-self.p))
            s1=segm(self.r,self.p,[self.q[0],phi])
            s2=segm(self.r,self.p,[phi,self.q[1]])
            return [s1,s2]
    def dist(self,pt):
        '''calculate the distance between pt and start point'''
        if self.r==0:
            return norm(pt-selp.p)
        else:
            phi=betw(self.q,atanv(pt-self.p))
            return abs(phi-self.q[0])
    __repr__ = __str__
def interll(s1,s2):
    '''intersection of two line segment. lm is lambda and mu'''
    if abs(cross(s1.p-s1.q,s2.p-s2.q))<infs:
        return []
    [l,m]=inv(array([s1.p-s1.q,s2.p-s2.q])).T.dot(s2.p-s1.q)
    if infs<l and l<1-infs and infs<m and m<1-infs:
        return [l*s1.p+(1-l)*s1.q]
    else:
        return []
def interlc(l,c):
    '''intersection of two line--circle'''
    #三角形的三条边
    a=l.p-c.p
    b=l.q-c.p
    delta=l.p-l.q
    direc=delta/norm(delta)#直线的方向矢量
    vert=a-dot(a,direc)*direc#垂线矢量，减去它相当于做投影变换
    h=norm(vert)#垂线长度height
    ah=a-vert
    bh=b-vert
    ret=[]
    if h<abs(c.r):
        con=sqrt(c.r**2-h**2)
        pm=[con*direc,-con*direc]
        for i in pm:
            theta=atanv(vert+i)
            #print(theta,c.q)
            if dot(ah-i,bh-i)<-infs and betw(c.q,theta+infs) and betw(c.q,theta-infs):
                ret.append(c.p+vert+i)
    return ret
def intercc(c1,c2):
    '''intersection of two circle--circle'''
    delta=c2.p-c1.p
    d=norm(delta)
    ret=[]
    if abs(c1.r-c2.r)+infs<d and d+infs<c1.r+c2.r:
        #print(c1.r,c2.r,d)
        t1=atanv(delta)
        t2=t1+pi
        phi1=arccos((c1.r**2+d**2-c2.r**2)/(2*c1.r*d))
        phi2=-arccos((c2.r**2+d**2-c1.r**2)/(2*c2.r*d))
        #print('theta',t1,t2)
        #print('phi:',phi1,phi2)
        for i in [1,-1]:
            if betw(c1.q,t1+i*phi1) and betw(c2.q,t2+i*phi2):
                ret.append(c1.p+rect(c1.r,t1+i*phi1))
    return ret
def intsecs(s1,s2):
    '''Judge and get the intersect point(if exist)
    Intersect with end point will not be considered intersect'''
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
