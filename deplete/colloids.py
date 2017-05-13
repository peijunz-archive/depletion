from .segm import (pol2rect, rect2pol, Segm, Circle, Line, norm2,
                   np, plt, rot_point, err)
from .curve import Curve
pie = np.pi / 6


def Trify(c):
    '''Construct curve with 3 copies'''
    return c + c.rotate(4 * pie) + c.rotate(8 * pie)
# Single colloid


def dartoid(r, h):
    '''Constrution a family of symmetric curves concave/convex'''
    cen = r * pol2rect(7 * pie)
    p0 = pol2rect(-pie)
    dr = p0 - cen
    t1 = rect2pol(dr)
    t2 = 2 * pie - 2 * t1
    R = norm2(dr)
    s1 = Circle(R, cen, [t1, t2])
    s2 = Line([0, 1], [0, 1 - 2 * h])
    cen1 = 2 * np.array([0, 1 - h]) - cen
    s3t = Circle(R, cen1, [t1 + t2 + 6 * pie, -t2])
    s4t = s1.rotate(4 * pie)
    p = Segm.intersect(s3t, s4t)[0]
    s3 = s3t.split(p)[0]
    s4 = s4t.split(p)[1]
    return Trify(Curve([s2, s3, s4]))


def triangloid(H=1, r=0, b=0, h=0):
    '''Constrution a family of curves with concave/convex'''
    top_left = [-b, H]
    top_right = [b, H]
    top = [0, H + h]
    center = r * pol2rect(7 * pie)
    z = rect2pol(rot_point(top_left, -4 * pie) - center)
    w = rect2pol(top_right - center)
    R = norm2(top_right - center)
    if abs(h) < err or abs(b) < err:
        L = [Circle(R, center, [z, w - z]), Line(top_right, top_left)]
    else:
        L = [Circle(R, center, [z, w - z]),
             Line(top_right, top), Line(top, top_left)]
    return Trify(Curve(L).clean())


def twistoid(dr=.7):
    '''Twisted lattice'''
    cen = pol2rect(7 * pie) + np.array([-dr, 0])
    p0 = [0, sqrt(9 / 4 + sqrt(3) * dr) - 0.5]
    R = sqrt(3) + dr
    t1 = arccos((sqrt(3) / 2 + dr) / (sqrt(3) + dr))
    s1 = Circle(R, cen, [0, t1])
    s2 = Line(p0, [0, 1])
    c1 = Curve([s1, s2])
    cen2 = cen + (sqrt(3) - dr) * pol2rect(4 * pie)
    s3 = Circle(R, cen2, [2 * pie - t1, 2 * pie])
    p1 = Segm.intersect(s1, s3)[0]
    s11 = s1.split(p1)[1]
    s31 = s3.split(p1)[1]
    c2 = Curve([s2.reverse(), s11.reverse(), s31]).shift(
        sqrt(3) * pol2rect(-2 * pie))
    return Trify(c1), Trify(c2)


def mark_intersect(l):
    '''mark intersection areas'''
    A = [i.area() for i in l]
    m = max(A)
    c = [(a > 0 and a != m) for a in A]
    return c


def show_colloids(c1, c2, off, inter=False):
    '''Show two colloids with their offset and intersection'''
    m, n = c1.offset(off), c2.offset(off)
    plt.axis('equal')
    c1.draw()
    c2.draw()
    m.draw('--')
    n.draw('--')
    if not inter:
        return 0
    sect = Curve.intersect(m, n)
    if sect:
        cs = list(sect)
        col = mark_intersect(cs)
        for co, seg in zip(col, cs):
            if co:
                seg.draw('r', linewidth=2)


def binding_energy(c1, c2, off):
    '''Binding energy is propto inter area of offset'''
    m, n = c1.offset(off), c2.offset(off)
    return Curve.interarea(m, n)
