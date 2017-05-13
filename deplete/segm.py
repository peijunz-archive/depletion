#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
err = 1e-10


def norm2(l):
    '''Norm of 2D vector
    >>> norm2([3,4])
    5.0
    '''
    return np.sqrt(l[0] * l[0] + l[1] * l[1])


def cross2(a, b):
    '''Cross product
    >>> cross2([1.,2],[2,1])
    -3.0
    '''
    return a[0] * b[1] - b[0] * a[1]


def pol2rect(theta):
    '''
    >>> pol2rect(np.pi/2)[1]==1.0
    True
    '''
    return np.array([np.cos(theta), np.sin(theta)])


def rect2pol(n):
    '''
    >>> rect2pol(pol2rect(np.pi/3))-np.pi/3
    0.0
    '''
    return np.arctan2(n[1], n[0])


def rot_point(p, theta):
    '''Rotation a point about origin
    >>> rect2pol(rot_point([0, 1.0], np.pi/3))-5*np.pi/6
    0.0
    '''
    c = np.cos(theta)
    s = np.sin(theta)
    M = np.array([[c, -s], [s, c]])
    return M.dot(p)


def normalize_angle(base, theta):
    """Get a angle in the (0, 2*np.pi) relative to base
    # TODO Edge cases?
    >>> normalize_angle(0, 0.01)
    0.01
    >>> normalize_angle(1, 2*np.pi)==2*np.pi-1
    True
    >>> normalize_angle(1, 2*np.pi+2)==1.0
    True
    """
    delta = theta - base
    return delta - np.floor(delta / (2 * np.pi)) * 2 * np.pi


def angle_between(q, theta):
    '''Find a theta+2*k*np.pi between range q
    >>> angle_between([np.pi/6, np.pi/3], np.pi/3)==np.pi/6
    True
    '''
    delta = normalize_angle(q[0], theta)
    if q[1] > 0:
        if err < delta < q[1] - err:
            return delta
    else:
        delta -= 2 * np.pi
        if -err > delta > q[1] + err:
            return delta
    return False


class Segm:
    '''Segment class, the base of Lines and Circles
    Also implemented intersection'''

    def __str__(self):
        return '{}{}{}'.format(self.start(), self.symbol(), self.end())

    def draw(self, *args, **kargs):
        """Draw the segment by matplotlib"""
        plt.plot(*(self.draw_raw().transpose()), *args, **kargs)

    @staticmethod
    def interll(s1, s2):
        '''intersection of two line segment'''
        a = s1.p - s1.q
        b = s2.p - s2.q
        c = s2.p - s1.q
        ab = cross2(a, b)
        if abs(ab) < err:
            return []
        bc = cross2(c, b)
        ca = cross2(a, c)
        l = bc / ab
        m = ca / ab
        if err < l < 1 - err and err < m < 1 - err:
            return [l * s1.p + (1 - l) * s1.q]
        else:
            return []

    @staticmethod
    def interlc(l, c):
        '''intersection of line--circle'''
        a = l.p - c.p
        b = l.q - c.p
        delta = l.p - l.q
        direc = delta / norm2(delta)  # 直线的方向矢量
        vert = a - np.dot(a, direc) * direc  # 垂线矢量，减去它相当于做投影变换
        h = norm2(vert)  # 垂线长度height
        ah = a - vert
        bh = b - vert
        ret = []
        if h < abs(c.r):
            con = np.sqrt(c.r**2 - h**2)
            pm = [con * direc, -con * direc]
            for i in pm:
                theta = rect2pol(vert + i)
                # print(theta,c.q)
                if np.dot(ah - i, bh - i) < -err and angle_between(c.q, theta):
                    ret.append(c.p + vert + i)
        return ret

    @staticmethod
    def intercc(c1, c2):
        '''intersection of circle--circle'''
        delta = c2.p - c1.p
        d = norm2(delta)
        ret = []
        if abs(c1.r - c2.r) + err < d < c1.r + c2.r - err:
            # print(c1.r,c2.r,d)
            t1 = rect2pol(delta)
            t2 = t1 + np.pi
            phi1 = np.arccos((c1.r**2 + d**2 - c2.r**2) / (2 * c1.r * d))
            phi2 = -np.arccos((c2.r**2 + d**2 - c1.r**2) / (2 * c2.r * d))
            for i in [1, -1]:
                if angle_between(c1.q, t1 + i * phi1) and\
                        angle_between(c2.q, t2 + i * phi2):
                    ret.append(c1.p + c1.r * pol2rect(t1 + i * phi1))
        return ret

    def intersect(s1, s2):
        '''Judge and get the intersect point(if exist)
        Intersect with end point will not be considered intersect
        Consider that under most condition they will not cross with others
        Calculate and compare the bound maybe useful
        This is the key function affect the performance
        '''
        if np.isinf(s1.r):
            if np.isinf(s2.r):
                return Segm.interll(s1, s2)
            else:
                return Segm.interlc(s1, s2)
        else:
            if np.isinf(s2.r):
                return Segm.interlc(s2, s1)
            else:
                return Segm.intercc(s1, s2)

    def scale(self):
        '''Full circle?'''
        return norm2(self.start() - self.end())


class Line(Segm):

    def __init__(self, start_p, end_p):
        '''Line Segment defined by starting point and ending point'''
        self.r = np.inf
        self.p = np.array(start_p)
        self.q = np.array(end_p)

    def __repr__(self):
        return 'Line({}, {})'.format(list(self.p), list(self.q))

    def _angle(self):
        '''Angle between line and x-axis'''
        v = self.q - self.p
        return np.arctan2(v[1], v[0])

    def _dir(self):
        '''Calculate direction vector of line'''
        v = self.q - self.p
        return v / norm2(v)

    def symbol(self):
        return '⟶'

    def start(self):
        '''Starting point'''
        return self.p

    def end(self):
        '''Ending point'''
        return self.q

    def start_angle(self):
        '''Starting tangent angle'''
        return self._angle()

    def end_angle(self):
        '''Ending tangent angle'''
        return self._angle()

    def start_dir(self):
        '''Starting tangent direction vector'''
        return self._dir()

    def end_dir(self):
        '''Ending tangent direction vector'''
        return self._dir()

    def draw_raw(self):
        return np.array([self.p, self.q])

    def shift(self, delta):
        '''add a delta=(dx,dy) shift to the segment'''
        return Line(self.p + delta, self.q + delta)

    def reverse(self):
        '''Reverse direction'''
        return Line(self.q, self.p)

    def rotate(self, t):
        '''Rotate the segment for t rad about origin=(0, 0)'''
        return Line(rot_point(self.p, t), rot_point(self.q, t))

    def split(self, pt):
        '''Use the node pt to separate the self into two segs'''
        return Line(self.p, pt), Line(pt, self.q)

    def extend(self, distance):
        '''Extend the line in the positive direction'''
        return self.shift(distance * rot_point(self.end_dir(), -np.pi / 2))

    def length(self):
        return norm2(self.p - self.q)


class Circle(Segm):
    _dots = 15

    def __init__(self, radius, center, angles):
        '''Use relative angle, r>0'''
        self.r = radius
        self.p = np.array(center)
        self.q = angles
        self.qend = self.q[0] + self.q[1]

    def set_dots(d=15):
        '''Set dots each radian'''
        Circle._dots = max(int(d), 6)

    def symbol(self):
        '''Mark direction of circles'''
        return '↺' if self.q[1] > 0 else '↻'

    def __repr__(self):
        return 'Circle({}, {}, {})'.format(self.r, list(self.p), list(self.q))

    def draw_raw(self):
        """Draw the segment by matplotlib, one point for 0.1 rad"""
        x, y = self.p
        z, w = self.q
        angles = z + np.linspace(0, w, 3 + Circle._dots * abs(w))
        A = np.empty([len(angles), 2])
        A[:, 0] = x + self.r * np.cos(angles)
        A[:, 1] = y + self.r * np.sin(angles)
        return A

    def start_angle(self):
        '''Starting tangent angle'''
        if self.q[1] > 0:
            return self.q[0] + np.pi / 2
        else:
            return self.q[0] - np.pi / 2

    def end_angle(self):
        '''Ending tangent angle'''
        if self.q[1] > 0:
            return self.qend + np.pi / 2
        else:
            return self.qend - np.pi / 2

    def start(self):
        '''Starting point'''
        return self.p + self.r * pol2rect(self.q[0])

    def end(self):
        '''Ending point'''
        return self.p + self.r * pol2rect(self.qend)

    def start_dir(self):
        '''Starting tangent direction vector'''
        return pol2rect(self.start_angle())

    def end_dir(self):
        '''Ending tangent direction vector'''
        return pol2rect(self.end_angle())

    def rotate(self, theta):
        '''Rotate the segment for t rad about origin=(0, 0)'''
        return Circle(self.r,
                      rot_point(self.p, theta),
                      [self.q[0] + theta, self.q[1]])

    def shift(self, delta):
        '''add a delta=(dx,dy) shift to the segment'''
        return Circle(self.r, self.p + delta, self.q)

    def split(self, pt):
        '''Use the node pt to separate the self into two segs'''
        phi = angle_between(self.q, rect2pol(pt - self.p))
        s1 = Circle(self.r, self.p, [self.q[0], phi])
        s2 = Circle(self.r, self.p, [self.q[0] + phi, self.q[1] - phi])
        return s1, s2

    def reverse(self):
        '''Reverse direction'''
        return Circle(self.r, self.p, [self.qend, -self.q[1]])

    def extend(self, dis):
        '''Extend the line in the positive direction
        TODO: If radius smaller than zero, simply connect?'''
        if self.q[1] > 0:
            return Circle(self.r + dis, self.p, self.q)
        else:
            return Circle(self.r - dis, self.p, self.q)

    def length(self):
        return abs(self.r * self.q[1])


if __name__ == "__main__":
    import doctest
    doctest.testmod()
