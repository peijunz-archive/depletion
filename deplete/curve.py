from .segm import *
class Curve(list):
    '''Curve class for boundaries of colloids'''
    def __init__(self, segm_list):  
        '''Initialize by segment list'''
        list.__init__([])
        self.extend(segm_list)
    def __str__(self):
        return '\n'.join([str(i) for i in self])
    def __add__(self, rhs):
        '''Use + to add multiple curves together while preserving Curve class'''
        return Curve(list.__add__(self, rhs))
    def draw(self, *args, **kargs):
        '''Draw the curve in matplotlib
        Use uniform properties for all'''
        points=np.concatenate([seg.draw_raw() for seg in self], axis=0).transpose()
        plt.plot(*points, *args, **kargs)
    def draw_each(self, *args, **kargs):
        '''draw the curve in matplotlib
        use new properties for each'''
        for seg in self:
            self.draw(*args, **kargs)
    def shift(self, delta):
        '''shift of a curve'''
        return Curve(seg.shift(delta) for seg in self)
    def rotate(self, t):
        '''rotation of a curve'''
        return Curve(seg.rotate(t) for seg in self)
    def area(self):
        '''area of a curve consists of line/circle segments'''
        poly=[]
        area=0
        for seg in self:
            area+=cross2(seg.start(), seg.end())/2
            if not np.isinf(seg.r):
                area+=(seg.q[1]-np.sin(seg.q[1]))*seg.r**2/2
        return area
    def offset_raw(self, dis):
        '''offset curve and dis is boldness of offset
        1. Get piles of smaller simple curves of its raw extension
        2. Get the final curve that has biggest area!
        '''
        l=len(self)
        ext=[]
        for i in range(l):
            former=self[i]
            latter=self[(i+1)%l]
            ff=former.extend(dis)
            ext.append(ff)
            sig=cross2(former.end_dir(), latter.start_dir())
            if abs(sig)<err:
                continue
            elif sig*dis<0:
                ext.append(Line(ff.end(), latter.extend(dis).start()))
            elif sig*dis>0:
                # Angles?
                delta=normalize_angle(former.end_angle(), latter.start_angle())
                angle=former.end_angle()-np.pi/2
                if dis<0:
                    delta-=2*np.pi
                    angle+=np.pi
                ext.append(Circle(abs(dis), former.end(), [angle, delta]))
        #Cut the redundant part of the curve
        return Curve(ext)
    def _split(self):
        '''Find the self intersecting points
        and use them to split curve into curves.
        It is designed for finding outmost boundaries'''
        return PCurve(self).cut().separate()
    def clean(self):
        return Curve(seg for seg in self if seg.scale()>err)
    def offset(self, d):
        '''Offset the out boundary.
        TODO add function for inner boundary by "negative" offset? '''
        cur=max(self.offset_raw(d)._split(), key=Curve.area)
        return cur.clean()
    @staticmethod
    def intersect(cur1, cur2):
        '''Intersect curves, and return simple curve components
        of the new curve'''
        pc=PCurve.join(PCurve(cur1), PCurve(cur2))
        l=len(pc)
        pc.cut()
        if len(pc)==l:
            return False
        return pc.separate()
    @staticmethod
    def interarea(cur1, cur2):
        '''Intersection area between two curves'''
        curs=Curve.intersect(cur1, cur2)
        if curs:
            areas=[cur.area() for cur in curs]
            pareas=[a for a in areas if a>0]
            return sum(pareas)-max(pareas)
        else:
            return False
class PCurve(list):
    '''Utility for Curve based on (start + curves + ending)'''
    def __init__(self, cur, shift=0):
        list.__init__([])
        l=len(cur)
        self.extend([seg.start(), seg, (i+1)%l] for i, seg in enumerate(cur))
    def cut(self):
        '''find all joint points and reconnect the pcurve according to it
        #TODO Improve efficiency using sweep-like method?'''
        i=1
        while i<len(self):
            for j in range(i):
                lis=Segm.intersect(self[j][1],self[i][1])
                if len(lis)==0:
                    continue
                if len(lis)==1:
                    p=lis[0]
                else:
                    #???
                    dis=[self[j][1].dist(k) for k in lis]
                    p=lis[dis.index(min(dis))]
                s1,s2=self[j][1].split(p)
                s3,s4=self[i][1].split(p)
                jnex=self[j][2]
                inex=self[i][2]
                self[i][1:]=[s3,len(self)]
                self.append([p,s2,jnex])
                self[j][1:]=[s1,len(self)]
                self.append([p,s4,inex])
            i+=1
        return self
    def separate(self):
        '''separate the different close simple curves'''
        l=len(self)
        sig=np.zeros(l)
        for i in range(l):
            if sig[i]==0:
                t=self[i]
                sig[i]=1
                tmp=[t[1]]
                while True:
                    j=t[2]
                    t=self[j]
                    if j==i:
                        break
                    sig[j]=1
                    tmp.append(t[1])
                yield Curve(tmp)
    @staticmethod
    def join(pc1, pc2):
        '''Join them into one'''
        l=len(pc1)
        for p in pc2:
            p[2]+=l
        pc=PCurve([])
        pc.extend(pc1+pc2)
        return pc
