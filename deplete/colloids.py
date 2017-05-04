from .curve import *
uni = np.pi / 6


def Trify(c):
    return c + c.rotate(4 * uni) + c.rotate(8 * uni)
# Single colloid


def dartoid(r, h):
    cen = r * pol2rect(np.pi + uni)
    p0 = pol2rect(-uni)
    dr = p0 - cen
    t1 = rect2pol(dr)
    t2 = np.pi / 3 - 2 * t1
    R = norm2(dr)
    s1 = Circle(R, cen, [t1, t2])
    s2 = Line([0, 1], [0, 1 - 2 * h])
    cen1 = 2 * np.array([0, 1 - h]) - cen
    s3t = Circle(R, cen1, [t1 + t2 + np.pi, -t2])
    s4t = s1.rotate(4 * uni)
    p = Segm.intersect(s3t, s4t)[0]
    s3 = s3t.split(p)[0]
    s4 = s4t.split(p)[1]
    return Trify(Curve([s2, s3, s4]))


def triangloid(H=1, r=0, b=0, h=0):
    '''Constrution a family of curves with concave/convex'''
    top_left = [-b, H]
    top_right = [b, H]
    top = [0, H + h]
    center = r * pol2rect(np.pi + uni)
    z = rect2pol(rot_point(top_left, -4 * uni) - center)
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
    cen = pol2rect(np.pi + uni) + np.array([-dr, 0])
    p0 = [0, sqrt(9 / 4 + sqrt(3) * dr) - 0.5]
    R = sqrt(3) + dr
    t1 = arccos((sqrt(3) / 2 + dr) / (sqrt(3) + dr))
    s1 = Circle(R, cen, [0, t1])
    s2 = Line(p0, [0, 1])
    c1 = Curve([s1, s2])
    cen2 = cen + (sqrt(3) - dr) * pol2rect(4 * uni)
    s3 = Circle(R, cen2, [2 * uni - t1, 2 * uni])
    p1 = Segm.intersect(s1, s3)[0]
    s11 = s1.split(p1)[1]
    s31 = s3.split(p1)[1]
    c2 = Curve([s2.reverse(), s11.reverse(), s31]).shift(
        sqrt(3) * pol2rect(-2 * uni))
    return Trify(c1), Trify(c2)


def mark_intersect(l):
    A = [i.area() for i in l]
    m = max(A)
    c = [(a > 0 and a != m) for a in A]
    return c


def show_colloids(c1, c2, off, inter=False):
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
    return Curve.interarea(m, n)
