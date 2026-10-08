"""2D hybrid C-mesh template of the S1223 section (mean chord, x-z plane).

* structured quad strip (QUAD_THICK thick) marched off the airfoil + wake cut
  (C-topology), first cell = H1 for y+ ~ 0.7
* unstructured triangles (gmsh) from the strip to the C-shaped far field
* triangles inside the airfoil (only used outboard of the tip)

Writes template2d.npz consumed by extrude.py.
"""
import math
import sys

import gmsh
import numpy as np
from scipy.interpolate import splev, splprep
from scipy.optimize import brentq
from scipy.spatial import cKDTree

import params as P


# ------------------------------------------------------------------ airfoil
def airfoil_surfaces():
    pts = np.loadtxt(P.AIRFOIL_FILE, skiprows=1) * P.C
    tck, u = splprep([pts[:, 0], pts[:, 1]], s=0, k=3, per=0)
    uu = np.linspace(0, 1, 200001)
    xx, zz = splev(uu, tck)
    i_le = int(np.argmin(xx))
    u_le = uu[i_le]

    def side(u0, u1):
        us = np.linspace(u0, u1, 100001)
        x, z = splev(us, tck)
        ds = np.hypot(np.diff(x), np.diff(z))
        s = np.concatenate([[0], np.cumsum(ds)])
        t = np.linspace(0, 1, P.N_SURF + 1)
        f = P.SINE_BLEND * (1 - np.cos(np.pi * t)) / 2 + (1 - P.SINE_BLEND) * t
        ui = np.interp(f * s[-1], s, us)
        x, z = splev(ui, tck)
        return np.column_stack([x, z])

    upper = side(u_le, 0.0)          # LE -> TE
    lower = side(u_le, 1.0)          # LE -> TE
    upper[0] = lower[0] = [xx[i_le], zz[i_le]]
    upper[-1] = lower[-1] = [P.C, 0.0]
    return upper, lower


def wake_x():
    k = np.arange(P.N_WAKE + 1)
    return P.C + P.WAKE_LEN * (1 - np.cos(np.pi * k / (2 * P.N_WAKE)))


def geometric_dist(h1, total, n):
    """Cumulative distances of n layers, first h1, summing to total."""
    if h1 * n >= total * 0.999:
        return np.linspace(0, total, n + 1)
    g = lambda r: h1 * (r ** n - 1) / (r - 1) - total
    r = brentq(g, 1.0 + 1e-9, 2.0)
    d = h1 * (r ** np.arange(n) )
    return np.concatenate([[0], np.cumsum(d)]) * total / d.sum()


# ------------------------------------------------------------- quad strip
def build_strip():
    upper, lower = airfoil_surfaces()
    xw = wake_x()
    wl = np.column_stack([xw[::-1], np.zeros_like(xw)])[:-1]   # outlet -> TE)
    wu = np.column_stack([xw, np.zeros_like(xw)])[1:]          # (TE -> outlet
    foil = np.vstack([lower[::-1], upper[1:]])                 # TE -> LE -> TE
    curve = np.vstack([wl, foil, wu])
    nwl, nf = len(wl), len(foil)
    kind = np.array([1] * nwl + [0] * nf + [1] * len(wu))      # 0 wall 1 wake
    i_te_lo, i_te_up = nwl, nwl + nf - 1

    # first-cell height: H1 on the wall, growing log-linearly along the
    # wake (in sqrt(s)) to a uniform distribution at the outlet
    s = np.where(kind == 1, np.abs(curve[:, 0] - P.C), 0.0)
    h_out = P.QUAD_THICK / P.N_LAYERS
    frac = np.sqrt(np.clip(s / P.WAKE_LEN, 0, 1))
    h1 = np.exp(np.log(P.H1) + frac * (np.log(h_out) - np.log(P.H1)))
    D = np.array([geometric_dist(h, P.QUAD_THICK, P.N_LAYERS) for h in h1])
    dh = np.diff(D, axis=1)

    M = len(curve)
    X = np.zeros((P.N_LAYERS + 1, M, 2))
    X[0] = curve
    for k in range(P.N_LAYERS):
        c = X[k]
        t = np.zeros_like(c)
        t[1:-1] = c[2:] - c[:-2]
        t[0] = c[1] - c[0]
        t[-1] = c[-1] - c[-2]
        t /= np.linalg.norm(t, axis=1)[:, None]
        n = np.column_stack([-t[:, 1], t[:, 0]])
        # distance-weighted normal smoothing: at height d the normal is the
        # Gaussian average over an arc-length window ~ d, so lines from
        # concave corners (upper TE) fan out instead of crossing
        sc = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(c, axis=0), axis=1))])
        sig = np.maximum(1.5 * D[:, k], 1e-6)
        wgt = np.exp(-0.5 * ((sc[None, :] - sc[:, None]) / sig[:, None]) ** 2)
        n = wgt @ n
        n /= np.linalg.norm(n, axis=1)[:, None]
        n[0] = [0, -1]
        n[-1] = [0, 1]
        nxt = c + n * dh[:, k:k + 1]
        if k > 3:
            # mild tangential (Laplacian) smoothing away from the wall
            w = min(0.5, 0.05 * (k - 3))
            lap = np.zeros_like(nxt)
            lap[1:-1] = 0.5 * (nxt[:-2] + nxt[2:]) - nxt[1:-1]
            # keep only the tangential part of the correction
            tt = np.zeros_like(nxt)
            tt[1:-1] = nxt[2:] - nxt[:-2]
            tt[0], tt[-1] = [1, 0], [1, 0]
            tt /= np.linalg.norm(tt, axis=1)[:, None]
            lap = (lap * tt).sum(1)[:, None] * tt
            # cap the sideways shift so thin cells are not sheared
            lm = np.linalg.norm(lap, axis=1)[:, None] * w
            cap = 0.25 * dh[:, k:k + 1]
            lap = np.where(lm > cap, lap * cap / np.maximum(lm, 1e-30), lap)
            # restore the prescribed layer height after the shift
            cand = nxt + w * lap
            dv = cand - c
            dist = np.linalg.norm(dv, axis=1)[:, None]
            cand = c + dv / dist * dh[:, k:k + 1]
            cand[0], cand[-1] = nxt[0], nxt[-1]
            nxt = cand
        X[k + 1] = nxt
    return X, kind, i_te_lo, i_te_up, upper, lower, h1


def quad_quality(X):
    a = X[:-1, :-1]
    b = X[:-1, 1:]
    c = X[1:, 1:]
    d = X[1:, :-1]
    def cr(p, q, r):
        return (q[..., 0] - p[..., 0]) * (r[..., 1] - p[..., 1]) - \
               (q[..., 1] - p[..., 1]) * (r[..., 0] - p[..., 0])
    return np.minimum.reduce([cr(a, b, d), cr(b, c, a), cr(c, d, b), cr(d, a, c)])


# ---------------------------------------------------------- gmsh triangles
def wake_spacing_expr():
    # local spacing of the half-cosine wake distribution at distance s
    L, N = P.WAKE_LEN, P.N_WAKE
    s = f"Max(x-{P.C},0)"
    return f"({math.pi/(2*N)})*Sqrt({s}*(2*{L}-Min({s},{L})))"


def triangulate(outer_strip, foil_closed):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("c")
    occ = gmsh.model.geo
    xo = P.C + P.WAKE_LEN
    R = P.R_FAR

    # strip outer curve (outlet-bottom -> around -> outlet-top), one element
    # per segment so its nodes are exactly the strip nodes
    sp = [occ.addPoint(x, z, 0) for x, z in outer_strip]
    sl = [occ.addLine(sp[i], sp[i + 1]) for i in range(len(sp) - 1)]
    pt = occ.addPoint(xo, R, 0, P.FAR_SIZE)
    pb = occ.addPoint(xo, -R, 0, P.FAR_SIZE)
    pc = occ.addPoint(P.C, 0, 0)
    pft = occ.addPoint(P.C, R, 0, P.FAR_SIZE)
    pfb = occ.addPoint(P.C, -R, 0, P.FAR_SIZE)
    pfl = occ.addPoint(P.C - R, 0, 0, P.FAR_SIZE)
    l_out_t = occ.addLine(sp[-1], pt)
    l_top = occ.addLine(pt, pft)
    a1 = occ.addCircleArc(pft, pc, pfl)
    a2 = occ.addCircleArc(pfl, pc, pfb)
    l_bot = occ.addLine(pfb, pb)
    l_out_b = occ.addLine(pb, sp[0])
    loop = occ.addCurveLoop(sl + [l_out_t, l_top, a1, a2, l_bot, l_out_b])
    s_fluid = occ.addPlaneSurface([loop])

    fp = [occ.addPoint(x, z, 0) for x, z in foil_closed[:-1]]
    fl = [occ.addLine(fp[i], fp[(i + 1) % len(fp)]) for i in range(len(fp))]
    s_foil = occ.addPlaneSurface([occ.addCurveLoop(fl)])
    occ.synchronize()
    for l in sl + fl:
        gmsh.model.mesh.setTransfiniteCurve(l, 2)
    # outlet lines: grow from the strip spacing to the far-field size
    h_wake_end = math.pi / (2 * P.N_WAKE) * P.WAKE_LEN
    n_out = 14
    for l in (l_out_t, l_out_b):
        gmsh.model.mesh.setTransfiniteCurve(l, n_out + 1, "Progression",
                                            1.25 if l == l_out_t else 1 / 1.25)

    f = gmsh.model.mesh.field
    fd = f.add("Distance")
    f.setNumbers(fd, "CurvesList", sl)
    f.setNumber(fd, "Sampling", 4)
    fm = f.add("MathEval")
    base = f"Max(0.0022, Min({h_wake_end}, {wake_spacing_expr()}))"
    f.setString(fm, "F", f"Min({P.FAR_SIZE}, {base} + {P.GROWTH}*F{fd})")
    f.setAsBackgroundMesh(fm)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Laplace2D")

    out = {}
    tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = coords.reshape(-1, 3)[:, :2]
    tmap = {t: i for i, t in enumerate(tags)}
    for name, surf in (("fluid", s_fluid), ("foil", s_foil)):
        et, etags, enodes = gmsh.model.mesh.getElements(2, surf)
        tri = np.array([tmap[n] for n in enodes[0]]).reshape(-1, 3)
        out[name] = tri
    edges = {}
    for name, curves in (("outlet", [l_out_t, l_out_b]),
                         ("farfield", [l_top, a1, a2, l_bot])):
        e = []
        for cv in curves:
            _, _, en = gmsh.model.mesh.getElements(1, cv)
            e.append(np.array([tmap[n] for n in en[0]]).reshape(-1, 2))
        edges[name] = np.vstack(e)
    gmsh.finalize()
    return coords, out, edges


# --------------------------------------------------------------- assemble
def main():
    X, kind, i_lo, i_up, upper, lower, h1col = build_strip()
    q = quad_quality(X)
    print(f"strip: {X.shape[1]} columns x {P.N_LAYERS} layers, "
          f"min corner cross = {q.min():.3e} (neg -> folded)")
    if q.min() <= 0:
        bad = np.argwhere(q <= 0)
        print("folded quads at (layer, column):", bad[:10])
        sys.exit(1)
    foil = np.vstack([lower[::-1], upper[1:]])   # TE->lower->LE->upper->TE
    tc, tris, tedges = triangulate(X[-1], foil)

    # ---- global node numbering: strip nodes, merging coincident ones
    nl1, M = X.shape[0], X.shape[1]
    allpts = [X.reshape(-1, 2), tc]
    pts = np.vstack(allpts)
    tree = cKDTree(pts)
    tol = 1e-9
    groups = tree.query_ball_point(pts, tol)
    rep = np.array([min(g) for g in groups])
    uniq, inv = np.unique(rep, return_inverse=True)
    nodes = pts[uniq]
    sid = inv[: nl1 * M].reshape(nl1, M)      # strip node ids
    gid = inv[nl1 * M:]                       # gmsh node ids

    # ---- cells (CCW in x-z)
    quads = np.stack([sid[:-1, :-1], sid[:-1, 1:], sid[1:, 1:], sid[1:, :-1]],
                     axis=-1).reshape(-1, 4)
    ftri = gid[tris["fluid"]]
    stri = gid[tris["foil"]]

    def ccw(cells):
        p = nodes[cells]
        a = 0.5 * np.sum(p[:, :, 0] * np.roll(p[:, :, 1], -1, 1)
                         - np.roll(p[:, :, 0], -1, 1) * p[:, :, 1], axis=1)
        cells = cells.copy()
        cells[a < 0] = cells[a < 0][:, ::-1]
        return cells, np.abs(a)

    quads, aq = ccw(quads)
    ftri, at = ccw(ftri)
    stri, ast = ccw(stri)
    print(f"quads {len(quads)}, fluid tris {len(ftri)}, foil tris {len(stri)}")
    print(f"min areas: quad {aq.min():.2e} tri {at.min():.2e} foil {ast.min():.2e}")

    # ---- tagged boundary edges
    wall = np.column_stack([sid[0, i_lo:i_up], sid[0, i_lo + 1:i_up + 1]])
    outlet = np.vstack([np.column_stack([sid[:-1, 0], sid[1:, 0]]),
                        np.column_stack([sid[:-1, -1], sid[1:, -1]]),
                        gid[tedges["outlet"]]])
    farfield = gid[tedges["farfield"]]

    # distance of every node from the airfoil (for the LE-wave morph)
    fd = cKDTree(np.vstack([foil[i] + (foil[i + 1] - foil[i]) * t
                            for i in range(len(foil) - 1)
                            for t in np.linspace(0, 1, 5)]))
    dist, _ = fd.query(nodes)
    inside = np.zeros(len(nodes), bool)
    inside[np.unique(stri)] = True
    dist[inside] = 0.0

    np.savez_compressed(
        "template2d.npz", nodes=nodes, quads=quads, ftri=ftri, stri=stri,
        wall=wall, outlet=outlet, farfield=farfield, dist=dist, strip=X, sid=sid, h1col=h1col,
        foil=foil, upper=upper, lower=lower)
    print(f"2D nodes {len(nodes)}, 2D fluid cells {len(quads) + len(ftri)}")
    return nodes, quads, ftri, stri


if __name__ == "__main__":
    main()
