"""
Hybrid C-mesh for the S1223 + NACA 0009 tandem configuration.

  * airfoil surfaces discretised with a sine (both-ends) clustering
  * 4 cm thick STRUCTURED quad layer around every airfoil, built here as a
    local C-grid (wall + short wake cut), first cell from the y+ estimate
  * outer region meshed by gmsh (quad-dominant, frontal-Delaunay + blossom)
    with a size that grows gradually from the quad-layer edge to >= 1 chord
    on the C-shaped far boundary
  * flat wake line, 30 chords long, with a cosine node distribution
  * 1-cell extrusion in z, written as msh2 for gmshToFoam

usage: python3 mesh.py <output.msh>
"""
import math
import sys
from collections import Counter

import gmsh
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree

import params as P


# ------------------------------------------------------------ airfoil data
def read_selig(path):
    pts = []
    with open(path) as f:
        for line in f:
            s = line.split()
            if len(s) == 2:
                try:
                    pts.append((float(s[0]), float(s[1])))
                except ValueError:
                    pass
    return np.array(pts)


def naca4_sym(t, n=400):
    """Closed-TE NACA 00xx, Selig order (TE -> upper -> LE -> lower -> TE)."""
    b = np.linspace(0.0, math.pi, n)
    x = 0.5 * (1.0 - np.cos(b))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2
                  + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    up = np.c_[x[::-1], yt[::-1]]
    lo = np.c_[x[1:], -yt[1:]]
    return np.vstack([up, lo])


def sine_dist(n, w=P.SINE_WEIGHT):
    """n elements on [0,1], small at both ends (sine/cosine law, blended
    with a uniform part so the LE/TE cells do not become vanishingly small)."""
    th = np.linspace(0.0, math.pi, n + 1)
    return w * 0.5 * (1.0 - np.cos(th)) + (1.0 - w) * th / math.pi


def resample_surface(xy, n):
    """n sine-clustered elements on one surface (LE->TE), cubic spline in arc length."""
    d = np.r_[0.0, np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]
    sx, sy = CubicSpline(d, xy[:, 0]), CubicSpline(d, xy[:, 1])
    s = sine_dist(n) * d[-1]
    return np.c_[sx(s), sy(s)]


def airfoil_surfaces(raw, chord, x0, n):
    """Upper and lower surfaces (LE->TE, n elements each, common LE/TE)."""
    raw = raw.copy()
    raw[:, 0] -= raw[:, 0].min()
    raw /= raw[:, 0].max()
    ile = int(np.argmin(raw[:, 0]))
    upper = raw[: ile + 1][::-1]
    lower = raw[ile:].copy()
    te = 0.5 * (upper[-1] + lower[-1])
    upper[-1] = te
    lower[-1] = te
    up = resample_surface(upper, n) * chord + [x0, 0.0]
    lo = resample_surface(lower, n) * chord + [x0, 0.0]
    lo[0] = up[0]
    lo[-1] = up[-1]
    return up, lo


def rotate(xy, deg, piv):
    c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    return (xy - piv) @ np.array([[c, s], [-s, c]]) + piv


# ------------------------------------------------------------ structured layer
def geometric(h1, total, n):
    """n cells, first h1, summing to total -> cumulative distances (n+1)."""
    if h1 * n >= total:
        return np.linspace(0.0, total, n + 1)
    lo_, hi_ = 1.0 + 1e-9, 2.0
    for _ in range(200):
        r = 0.5 * (lo_ + hi_)
        if h1 * (r ** n - 1) / (r - 1) > total:
            hi_ = r
        else:
            lo_ = r
    d = h1 * (r ** np.arange(n + 1) - 1) / (r - 1)
    return d * total / d[-1]


def c_block(up, lo, h1):
    """Local C-grid (wall + wake cut) of thickness BL_THICK.

    Path (fluid on the left): wake end -> TE -> lower -> LE -> upper -> TE -> wake end.
    Returns node array X[i, j] (i along path, j wall-normal) and info.
    """
    T = P.BL_THICK
    nl = int(math.ceil(math.log(1 + T * (P.BL_RATIO - 1) / h1) / math.log(P.BL_RATIO)))
    te = up[-1]
    # wake cut along the trailing-edge bisector, geometric spacing from TE
    w = (te - up[-2]) / np.linalg.norm(te - up[-2]) + (te - lo[-2]) / np.linalg.norm(te - lo[-2])
    w /= np.linalg.norm(w)
    a = np.linalg.norm(up[-1] - up[-2])
    wk = [0.0]
    while wk[-1] < T:
        wk.append(wk[-1] + a * 1.15 ** (len(wk) - 1))
    wk = np.array(wk) * T / wk[-1]
    wake = te + wk[1:, None] * w
    m, n = len(wake), len(up) - 1
    path = np.vstack([wake[::-1], lo[::-1], up[1:], wake])
    i_te_up = m + 2 * n
    # marching distances: wall columns geometric from h1; along the wake cut
    # the first cell blends towards uniform so the block exit is not tiny
    dist = np.empty((len(path), nl + 1))
    for i in range(len(path)):
        if i < m or i > i_te_up:
            s = wk[m - i] if i < m else wk[i - i_te_up]
            f = (s / T) ** 1.5
            hh = h1 * (1 - f) + (T / nl) * f
        else:
            hh = h1
        dist[i] = geometric(hh, T, nl)
    # vertex normals (left of the path direction), averaged at corners
    tang = np.diff(path, axis=0)
    tang /= np.linalg.norm(tang, axis=1)[:, None]
    segn = np.c_[-tang[:, 1], tang[:, 0]]
    nrm = np.empty_like(path)
    nrm[0], nrm[-1] = segn[0], segn[-1]
    nrm[1:-1] = segn[:-1] + segn[1:]
    nrm /= np.linalg.norm(nrm, axis=1)[:, None]
    # heavily smoothed directions for the outer part of the layer
    sm = nrm.copy()
    for _ in range(400):
        sm[1:-1] = 0.5 * sm[1:-1] + 0.25 * (sm[:-2] + sm[2:])
        sm /= np.linalg.norm(sm, axis=1)[:, None]
    # outer row: offset along the smoothed directions, then the nodes are
    # locally re-spaced (Gaussian-smoothed arc-length parameter) so the dense
    # LE/TE columns spread out and cannot cross
    raw = path + T * sm
    s_raw = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(raw, axis=0), axis=1))]
    s_new = gaussian_filter1d(s_raw, sigma=P.BL_SPREAD, mode="nearest")
    s_new = (s_new - s_new[0]) / (s_new[-1] - s_new[0]) * s_raw[-1]
    outer = np.c_[np.interp(s_new, s_raw, raw[:, 0]), np.interp(s_new, s_raw, raw[:, 1])]
    # columns leave the wall along the (lightly smoothed) local normal and
    # turn quickly onto the straight line wall -> outer node
    nw = nrm.copy()
    for _ in range(5):
        nw[1:-1] = 0.5 * nw[1:-1] + 0.25 * (nw[:-2] + nw[2:])
    nw /= np.linalg.norm(nw, axis=1)[:, None]
    X = np.empty((len(path), nl + 1, 2))
    for j in range(nl + 1):
        d = dist[:, j][:, None]
        w = (1.0 - d / T) ** 4
        X[:, j] = path + d * (w * nw + (1.0 - w) * (outer - path) / T)
    return X, dict(nl=nl, m=m, n=n, i_le=m + n)


def quads_of(X, offset):
    ni, nj = X.shape[:2]
    idx = offset + np.arange(ni * nj).reshape(ni, nj)
    q = np.stack([idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]], axis=-1)
    return X.reshape(-1, 2), q.reshape(-1, 4)


def block_outline(X):
    """Closed outer polygon of a C-block: exit column (low), outer row, exit column (up)."""
    nl = X.shape[1] - 1
    a = [X[0, j] for j in range(nl + 1)]
    b = [X[i, nl] for i in range(1, X.shape[0])]
    c = [X[-1, j] for j in range(nl - 1, 0, -1)]
    return np.array(a + b + c)


# ------------------------------------------------------------ build
def build(out):
    up1, lo1 = airfoil_surfaces(read_selig(P.ROOT + "/s1223.dat"), P.C1, P.X_LE1, P.N_SURF1)
    up2, lo2 = airfoil_surfaces(naca4_sym(0.09), P.C2, P.X_LE2, P.N_SURF2)
    if abs(P.AOA_DEG) > 1e-12:     # geometric AoA: rotate both about the S1223 TE
        piv = np.array([P.X_TE1, 0.0])
        up1, lo1, up2, lo2 = (rotate(a, P.AOA_DEG, piv) for a in (up1, lo1, up2, lo2))

    X1, i1 = c_block(up1, lo1, P.H1_1)
    X2, i2 = c_block(up2, lo2, P.H1_2)
    for X, name in ((X1, "S1223"), (X2, "NACA0009")):
        _, q = quads_of(X, 0)
        p = X.reshape(-1, 2)
        a = 0.5 * ((p[q[:, 1], 0] - p[q[:, 0], 0]) * (p[q[:, 3], 1] - p[q[:, 0], 1])
                   - (p[q[:, 3], 0] - p[q[:, 0], 0]) * (p[q[:, 1], 1] - p[q[:, 0], 1]))
        print("%s quad layer: %d x %d cells, min area %.3e" %
              (name, X.shape[0] - 1, X.shape[1] - 1, a.min()))
        if a.min() <= 0:
            raise RuntimeError("inverted cell in the %s quad layer" % name)
    out1, out2 = block_outline(X1), block_outline(X2)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("tandem")
    geo = gmsh.model.geo

    def polygon(poly):
        pts = [geo.addPoint(x, y, 0.0) for x, y in poly]
        lines = [geo.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
        return pts, lines

    pts1, lines1 = polygon(out1)
    pts2, lines2 = polygon(out2)

    xc, R = P.X_TE1, P.R_FAR
    y_w = out2[0, 1]                     # wake line height (NACA block exit)
    pc = geo.addPoint(xc, 0, 0)
    pa_top = geo.addPoint(xc, R, 0)
    pa_front = geo.addPoint(xc - R, 0, 0)
    pa_bot = geo.addPoint(xc, -R, 0)
    po_top = geo.addPoint(P.X_OUT, R, 0)
    po_mid = geo.addPoint(P.X_OUT, y_w, 0)
    po_bot = geo.addPoint(P.X_OUT, -R, 0)
    arc1 = geo.addCircleArc(pa_top, pc, pa_front)
    arc2 = geo.addCircleArc(pa_front, pc, pa_bot)
    bot = geo.addLine(pa_bot, po_bot)
    o1 = geo.addLine(po_bot, po_mid)
    o2 = geo.addLine(po_mid, po_top)
    top = geo.addLine(po_top, pa_top)
    outer = geo.addCurveLoop([arc1, arc2, bot, o1, o2, top])
    surf = geo.addPlaneSurface([outer, geo.addCurveLoop(lines1), geo.addCurveLoop(lines2)])

    # ---------------- size function
    def edge_h(poly):
        d = np.hypot(*np.diff(np.vstack([poly, poly[:1]]), axis=0).T)
        return 0.5 * (d + np.roll(d, 1))

    dense, dense_h = [], []
    for poly in (out1, out2):
        q = np.vstack([poly, poly[:1]])
        h = edge_h(poly)
        hh = np.r_[h, h[:1]]
        for t in np.linspace(0, 1, 4, endpoint=False):
            dense.append(q[:-1] * (1 - t) + q[1:] * t)
            dense_h.append(hh[:-1] * (1 - t) + hh[1:] * t)
    dense, dense_h = np.vstack(dense), np.concatenate(dense_h)
    tree = cKDTree(dense)
    dth = 0.5 * math.pi / P.N_WAKE
    ds_min = P.L_WAKE * (1.0 - math.cos(dth))

    def wake_size(x, y):
        """Cosine law along the 30c wake: x' = L (1 - cos theta), theta uniform."""
        xp = x - P.X_TE1
        if xp < 0.0:
            return 1e9
        th = math.acos(1.0 - min(xp, P.L_WAKE) / P.L_WAKE)
        return max(P.L_WAKE * math.sin(th) * dth, ds_min) + P.GROWTH * abs(y - y_w)

    def size(x, y):
        d, i = tree.query((x, y))
        return min(P.FAR_SIZE, dense_h[i] + P.GROWTH * d, wake_size(x, y))

    def size_cb(dim, tag, x, y, z, lc):
        return float(size(x, y))

    # ---------------- embedded flat wake line (30c wake behind the tandem)
    def march(x0, x1, y):
        xs = [x0]
        while xs[-1] < x1:
            xs.append(xs[-1] + size(xs[-1], y))
        if len(xs) <= 2:
            return np.array([x0, x1])
        xs = np.array(xs)
        return x0 + (xs - x0) * (x1 - x0) / (xs[-1] - x0)

    def wake_line(p_start, x_start, x_end, p_end, y):
        xs = march(x_start, x_end, y)
        wp = [p_start] + [geo.addPoint(x, y, 0) for x in xs[1:-1]] + [p_end]
        return [geo.addLine(wp[k], wp[k + 1]) for k in range(len(wp) - 1)]

    wake_lines = wake_line(pts2[0], out2[0, 0], P.X_OUT, po_mid, y_w)

    geo.synchronize()
    for c in lines1 + lines2 + wake_lines:
        gmsh.model.mesh.setTransfiniteCurve(c, 2)
    gmsh.model.mesh.embed(1, wake_lines, 2, surf)
    gmsh.model.mesh.setSizeCallback(size_cb)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 8)          # frontal-Delaunay for quads
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
    gmsh.option.setNumber("Mesh.Smoothing", 10)
    gmsh.model.mesh.generate(2)

    ntag, ncoord, _ = gmsh.model.mesh.getNodes()
    gpts = ncoord.reshape(-1, 3)[:, :2]
    tag2i = {t: k for k, t in enumerate(ntag)}
    gel = []
    for ty, _, en in zip(*gmsh.model.mesh.getElements(2, surf)):
        nn = 3 if ty == 2 else 4
        for e in en.reshape(-1, nn):
            gel.append([tag2i[t] for t in e])
    gmsh.finalize()

    # ---------------- merge structured layers + gmsh region
    p1, q1 = quads_of(X1, 0)
    p2, q2 = quads_of(X2, len(p1))
    allp = np.vstack([p1, p2, gpts])
    cells = [list(c) for c in q1] + [list(c) for c in q2] + \
        [[k + len(p1) + len(p2) for k in c] for c in gel]
    key = np.round(allp / 1e-9).astype(np.int64)
    _, uniq, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    pts = allp[uniq]
    cells = [[int(inv[k]) for k in c] for c in cells]
    for c in cells:                      # counter-clockwise
        xy = pts[c]
        a = np.sum(xy[:, 0] * np.roll(xy[:, 1], -1) - np.roll(xy[:, 0], -1) * xy[:, 1])
        if a < 0:
            c.reverse()
        if len(set(c)) != len(c):
            raise RuntimeError("degenerate cell %s" % c)

    # boundary edges
    ec = Counter()
    for c in cells:
        for k in range(len(c)):
            ec[tuple(sorted((c[k], c[(k + 1) % len(c)])))] += 1
    if max(ec.values()) > 2:
        raise RuntimeError("non-manifold edge in the merged mesh")
    bedges = [e for e, v in ec.items() if v == 1]

    def wall_nodes(X, info, off):
        col0 = off + np.arange(X.shape[0]) * X.shape[1]
        return set(inv[col0[info["m"]: info["m"] + 2 * info["n"] + 1]].tolist())

    wall1, wall2 = wall_nodes(X1, i1, 0), wall_nodes(X2, i2, len(p1))
    patches = {"airfoil1": [], "airfoil2": [], "outlet": [], "farfield": []}
    rr = np.hypot(pts[:, 0] - xc, pts[:, 1])
    for a, b in bedges:
        if a in wall1 and b in wall1:
            patches["airfoil1"].append((a, b))
        elif a in wall2 and b in wall2:
            patches["airfoil2"].append((a, b))
        elif abs(pts[a, 0] - P.X_OUT) < 1e-6 and abs(pts[b, 0] - P.X_OUT) < 1e-6:
            patches["outlet"].append((a, b))
        elif abs(abs(pts[a, 1]) - R) < 1e-6 or abs(rr[a] - R) < 1e-4:
            patches["farfield"].append((a, b))
        else:
            raise RuntimeError("unclassified boundary edge at %s" % pts[a])
    for k, v in patches.items():
        print("  patch %-9s %6d faces" % (k, len(v)))
    nq = sum(1 for c in cells if len(c) == 4)
    print("2-D cells: %d  (quads %d, triangles %d), points %d"
          % (len(cells), nq, len(cells) - nq, len(pts)))
    if len(patches["airfoil1"]) != 2 * P.N_SURF1 or len(patches["airfoil2"]) != 2 * P.N_SURF2:
        raise RuntimeError("wall face count mismatch")

    write_msh2(out, pts, cells, patches)
    np.savez(out.replace(".msh", "") + "_geom.npz", up1=up1, lo1=lo1, up2=up2, lo2=lo2,
             ncells=len(cells), nquad=nq, nl1=i1["nl"], nl2=i2["nl"])
    return len(cells)


def write_msh2(path, pts, cells, patches):
    """Extrude the 2-D mesh one cell in z and write a gmsh 2.2 ascii file."""
    n = len(pts)
    names = list(patches) + ["frontAndBack", "fluid"]
    tag = {k: i + 1 for i, k in enumerate(names)}
    lines = []
    eid = 0
    for name, edges in patches.items():
        for a, b in edges:
            eid += 1
            lines.append("%d 3 2 %d %d %d %d %d %d" % (eid, tag[name], tag[name],
                                                       a + 1, b + 1, b + 1 + n, a + 1 + n))
    for c in cells:
        for zoff in (0, n):
            eid += 1
            ty = 2 if len(c) == 3 else 3
            lines.append("%d %d 2 %d %d %s" % (eid, ty, tag["frontAndBack"], tag["frontAndBack"],
                                              " ".join(str(k + 1 + zoff) for k in c)))
    for c in cells:
        eid += 1
        ty = 6 if len(c) == 3 else 5
        nodes = [k + 1 for k in c] + [k + 1 + n for k in c]
        lines.append("%d %d 2 %d %d %s" % (eid, ty, tag["fluid"], tag["fluid"],
                                          " ".join(map(str, nodes))))
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n%d\n" % len(names))
        for k in names:
            f.write("%d %d \"%s\"\n" % (3 if k == "fluid" else 2, tag[k], k))
        f.write("$EndPhysicalNames\n$Nodes\n%d\n" % (2 * n))
        for z in (0.0, P.SPAN):
            for i, (x, y) in enumerate(pts):
                f.write("%d %.12g %.12g %.12g\n" % (i + 1 + (n if z else 0), x, y, z))
        f.write("$EndNodes\n$Elements\n%d\n" % len(lines))
        f.write("\n".join(lines))
        f.write("\n$EndElements\n")


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "tandem.msh")
