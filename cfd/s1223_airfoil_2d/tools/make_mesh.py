#!/usr/bin/env python3
"""
2D hybrid C-mesh for the S1223 airfoil with a blunt trailing edge, written
directly as an OpenFOAM polyMesh (one cell deep, front/back = empty).

  * Airfoil cut straight at x/c = teCut (blunt base, treated as wall).
  * Structured QUAD band (C-topology), `quadThickness` thick, around the
    airfoil and along the wake to `wakeChords` chords behind the TE:
      - `nSurface` elements on each surface, sine-type (cosine) spacing:
        small at the LE and TE, largest at mid-chord;
      - wake: one-sided cosine spacing, small at the TE, growing downstream;
      - wall-normal: geometric growth from a first cell sized for
        y+ = yPlusTarget at the first cell centre (flat-plate Cf estimate).
  * TRIANGLES everywhere else (gmsh Frontal-Delaunay), growing from the band
    edge size to `triMaxSize` (>= chord) at the far field; they share the
    band's outer nodes (conformal interface).
  * Mesh quality is measured per cell (equiangle and orthogonal quality,
    1 = ideal) and written to constant/meshQuality.txt.

Usage:  python3 tools/make_mesh.py [case_dir]
Needs:  numpy, scipy, gmsh  (pip install numpy scipy gmsh)
"""
import math
import os
import re
import sys

import numpy as np
from scipy.interpolate import splev, splprep
from scipy.optimize import brentq, minimize_scalar
from scipy.spatial import cKDTree

CASE = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else
                       os.path.join(os.path.dirname(__file__), ".."))
AIRFOIL = os.environ.get("AIRFOIL_DAT", os.path.join(CASE, "..", "..", "s1223.dat"))


# ----------------------------------------------------------------- parameters
def read_parameters(path):
    """Plain `name number;` entries of an OpenFOAM dictionary."""
    num = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
    out = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(rf"^\s*(\w+)\s+({num})\s*;", line)
            if m:
                out[m.group(1)] = float(m.group(2))
    return out


P = read_parameters(os.path.join(CASE, "caseParameters"))
C = P["chord"]
NU = P["mu"] / P["rho"]
UINF = P["Re"] * NU / C
N_SURF = int(P["nSurface"])
N_BASE = int(P["nBase"])
N_LAY = int(P["nLayers"])
N_WAKE = int(P["nWake"])
T_BAND = P["quadThickness"]
L_WAKE = P["wakeChords"] * C
R_FAR = P["farfieldChords"] * C
DZ = P["spanDepth"]
WAKE_BEND = 0.25 * C      # length over which the wake cut turns to horizontal [m]
TRI_ALGO = int(os.environ.get("MESH_ALGO", 6))           # gmsh 2D algorithm
TRI_OPT = [m for m in os.environ.get("MESH_OPT", "Laplace2D").split(",") if m]
CUT_ANGLE = float(os.environ.get("MESH_CUT", 1.0))   # wake-cut start angle / TE bisector angle
SIG = float(os.environ.get("MESH_SIG", 1.5))     # normal smoothing width / wall distance
W_MAX = float(os.environ.get("MESH_W", 0.25))     # tangential relaxation weight
SWEEPS = int(os.environ.get("MESH_SWEEPS", 2))   # relaxation sweeps per layer


def first_cell_height():
    """Schlichting flat-plate Cf -> u_tau -> wall distance for the y+ target.
    OpenFOAM evaluates y+ at the first cell CENTRE, so the cell is twice that."""
    cf = (2.0 * math.log10(P["Re"]) - 0.65) ** -2.3
    utau = UINF * math.sqrt(cf / 2.0)
    y = P["yPlusTarget"] * NU / utau
    return 2.0 * y * P.get("firstCellFactor", 1.0), cf, utau, y


H1, CF, UTAU, Y_WALL = first_cell_height()


# ------------------------------------------------------------------- airfoil
def airfoil_surfaces():
    """Sine-spaced points on the airfoil cut at x/c = teCut.
    Returns lower (TE corner -> LE), upper (LE -> TE corner)."""
    xy = np.loadtxt(AIRFOIL, skiprows=1) * C        # Selig: TE upper..LE..TE lower
    tck, _ = splprep(xy.T, s=0, k=3)
    ud = np.linspace(0.0, 1.0, 400001)
    xd, yd = splev(ud, tck)
    i = int(np.argmin(xd))
    u_le = minimize_scalar(lambda u: float(splev(u, tck)[0]),
                           bounds=(ud[i - 5], ud[i + 5]), method="bounded",
                           options=dict(xatol=1e-12)).x
    x_cut = P["teCut"] * C
    fx = lambda u: float(splev(u, tck)[0]) - x_cut
    u_up = brentq(fx, 0.0, 0.05)              # upper surface: u from 0 (TE)
    u_lo = brentq(fx, 0.95, 1.0)              # lower surface: u up to 1 (TE)
    seg = np.hypot(np.diff(xd), np.diff(yd))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    s_of = lambda u: float(np.interp(u, ud, s))

    def surface(s0, s1):
        k = np.arange(N_SURF + 1)
        sk = s0 + (s1 - s0) * 0.5 * (1.0 - np.cos(np.pi * k / N_SURF))
        return np.column_stack(splev(np.interp(sk, s, ud), tck))

    upper = surface(s_of(u_le), s_of(u_up))   # LE -> upper TE corner
    lower = surface(s_of(u_lo), s_of(u_le))   # lower TE corner -> LE
    upper[0] = lower[-1]
    lower[0, 0] = upper[-1, 0] = x_cut        # exactly vertical base
    return lower, upper


# ---------------------------------------------------------------- quad band
def growth_ratio(h1, total, n):
    """r such that h1 (r^n - 1)/(r - 1) = total."""
    f = lambda r: h1 * (r ** n - 1.0) / (r - 1.0) - total
    return brentq(f, 1.0 + 1e-12, 3.0)


def layer_distances(h1, total, n):
    r = growth_ratio(h1, total, n)
    j = np.arange(n + 1)
    return h1 * (r ** j - 1.0) / (r - 1.0), r


def arc(F):
    return np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(F, axis=0).T))])


def front_normals(F):
    d = np.diff(F, axis=0)
    t = d / np.linalg.norm(d, axis=1)[:, None]
    T = np.empty_like(F)
    T[1:-1] = t[:-1] + t[1:]
    T[0], T[-1] = t[0], t[-1]
    T /= np.linalg.norm(T, axis=1)[:, None]
    n = np.column_stack([-T[:, 1], T[:, 0]])  # left normal = away from the wall
    n[0], n[-1] = (0.0, -1.0), (0.0, 1.0)     # stay on the outlet plane
    return n


def smooth_normals(F, n, sigma):
    """Gaussian smoothing of front normals along the front arc length."""
    s = arc(F)
    out = n.copy()
    for i in range(1, len(F) - 1):
        sg = sigma[i]
        if sg <= 0.0:
            continue
        a, b = np.searchsorted(s, [s[i] - 4 * sg, s[i] + 4 * sg])
        w = np.exp(-0.5 * ((s[a:b] - s[i]) / sg) ** 2)
        v = (w[:, None] * n[a:b]).sum(axis=0)
        out[i] = v / np.linalg.norm(v)
    return out


def redistribute(G, w, sweeps):
    """Relax node spacing along the front (Laplace on arc length), ends fixed."""
    s = arc(G)
    s0 = s.copy()
    for _ in range(sweeps):
        s[1:-1] += w[1:-1] * (0.5 * (s[:-2] + s[2:]) - s[1:-1])
    return np.column_stack([np.interp(s, s0, G[:, 0]), np.interp(s, s0, G[:, 1])])


def wake_h1(dx):
    """Wake-centreline cell height: kept equal to the wall value. Wake cells
    are very long and thin (aspect ratio 100-1000); any growth of their
    height along the wake tapers them and moves the centroid off the face
    normals (poor orthogonal quality), so the height stays constant."""
    return np.full_like(np.asarray(dx, dtype=float), H1)


def build_band(lower, upper):
    """C-line: lower wake branch -> lower surface -> upper surface -> upper
    wake branch. The blunt base is NOT on the C-line: the two wake branches
    start at the TE corners and enclose a thin slit block (see build_slit),
    so the C-line has no sharp corners."""
    x_te = lower[0, 0]
    y_m = 0.5 * (lower[0, 1] + upper[-1, 1])
    base = upper[-1, 1] - lower[0, 1]
    k = np.arange(0, N_WAKE + 1)
    xw = x_te + L_WAKE * (1.0 - np.cos(0.5 * np.pi * k / N_WAKE))
    half = 0.5 * (base + N_BASE * (wake_h1(xw - x_te) - H1))   # slit half width
    # wake cut: leaves the TE along the TE bisector, then relaxes to the
    # freestream direction (a kinked cut would shear the cells at the TE)
    t_up = upper[-1] - upper[-2]
    t_lo = lower[0] - lower[1]
    slope = math.tan(CUT_ANGLE * 0.5 * (math.atan2(t_up[1], t_up[0]) + math.atan2(t_lo[1], t_lo[0])))
    yc = y_m + slope * WAKE_BEND * (1.0 - np.exp(-(xw - x_te) / WAKE_BEND))
    wake_lo = np.column_stack([xw, yc - half])[1:]
    wake_up = np.column_stack([xw, yc + half])[1:]
    line = np.vstack([wake_lo[::-1], lower, upper[1:], wake_up])
    ni = len(line)
    wall = (N_WAKE, N_WAKE + 2 * N_SURF)                      # TE corners

    h1 = np.full(ni, H1)
    for idx in list(range(0, wall[0])) + list(range(wall[1] + 1, ni)):
        h1[idx] = wake_h1(line[idx, 0] - x_te)
    dist = np.empty((ni, N_LAY + 1))
    ratios = np.empty(ni)
    for i in range(ni):
        dist[i], ratios[i] = layer_distances(h1[i], T_BAND, N_LAY)

    Pb = np.empty((N_LAY + 1, ni, 2))
    Pb[0] = line
    for j in range(N_LAY):
        F = Pb[j]
        n = smooth_normals(F, front_normals(F), SIG * dist[:, j])
        G = F + (dist[:, j + 1] - dist[:, j])[:, None] * n
        dj = dist[:, j + 1]
        w = W_MAX * np.clip((dj - 0.0005) / 0.01, 0.0, 1.0)
        Pb[j + 1] = redistribute(G, w, sweeps=SWEEPS)
    return Pb, line, wall, ratios


# ------------------------------------------------------------- triangles
def triangulate(ring, x_out):
    import gmsh

    seg = np.hypot(*np.diff(ring, axis=0).T)
    h_ring = np.empty(len(ring))
    h_ring[1:-1] = 0.5 * (seg[:-1] + seg[1:])
    h_ring[0], h_ring[-1] = seg[0], seg[-1]
    tree = cKDTree(ring)
    g, hmax = P["triGrowth"], P["triMaxSize"]
    kq = min(64, len(ring))

    def size(dim, tag, x, y, z, lc):
        d, idx = tree.query((x, y), k=kq)
        return float(min(hmax, np.min(h_ring[idx] + g * d)))

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("tri")
    geo = gmsh.model.geo
    rp = [geo.addPoint(x, y, 0.0) for x, y in ring]
    pB = geo.addPoint(x_out, -R_FAR, 0.0)
    pC = geo.addPoint(0.0, -R_FAR, 0.0)
    pO = geo.addPoint(0.0, 0.0, 0.0)
    pW = geo.addPoint(-R_FAR, 0.0, 0.0)
    pD = geo.addPoint(0.0, R_FAR, 0.0)
    pE = geo.addPoint(x_out, R_FAR, 0.0)
    ring_lines = [geo.addLine(rp[k], rp[k + 1]) for k in range(len(rp) - 1)]
    outer = [geo.addLine(rp[-1], pE), geo.addLine(pE, pD),
             geo.addCircleArc(pD, pO, pW), geo.addCircleArc(pW, pO, pC),
             geo.addLine(pC, pB), geo.addLine(pB, rp[0])]
    geo.addPlaneSurface([geo.addCurveLoop(ring_lines + outer)])
    geo.synchronize()
    for t in ring_lines:
        gmsh.model.mesh.setTransfiniteCurve(t, 2)
    for opt, val in [("Mesh.MeshSizeExtendFromBoundary", 0),
                     ("Mesh.MeshSizeFromPoints", 0),
                     ("Mesh.MeshSizeFromCurvature", 0),
                     ("Mesh.Algorithm", TRI_ALGO),
                     ("Mesh.Smoothing", 10)]:
        gmsh.option.setNumber(opt, val)
    gmsh.model.mesh.setSizeCallback(size)
    gmsh.model.mesh.generate(2)
    for method in TRI_OPT:
        gmsh.model.mesh.optimize(method, niter=5)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    _, _, enodes = gmsh.model.mesh.getElements(2)
    tri_nodes = enodes[0].reshape(-1, 3)
    gmsh.finalize()

    order = {int(t): i for i, t in enumerate(tags)}
    V = xyz.reshape(-1, 3)[:, :2].copy()
    T = np.vectorize(order.get)(tri_nodes.astype(np.int64))
    return V, T


# --------------------------------------------------------------- assembly
def assemble(Pb, wall, V, T):
    nj, ni = Pb.shape[0] - 1, Pb.shape[1]
    i_lc, i_uc = wall                         # lower / upper TE corner
    pts = []

    def node(xy):
        pts.append((float(xy[0]), float(xy[1])))
        return len(pts) - 1

    nid = np.array([[node(Pb[j, i]) for i in range(ni)] for j in range(nj + 1)])
    wall_nodes = set(int(nid[0, i]) for i in range(i_lc, i_uc + 1))

    cells = []
    for j in range(nj):
        for i in range(ni - 1):
            cells.append([int(nid[j, i]), int(nid[j, i + 1]),
                          int(nid[j + 1, i + 1]), int(nid[j + 1, i])])

    # slit block between the two wake branches; its first column face is
    # the blunt base (wall)
    for c in range(N_WAKE + 1):
        lo, up = nid[0, i_lc - c], nid[0, i_uc + c]
        col = [int(lo)]
        for r in range(1, N_BASE):
            col.append(node(Pb[0, i_lc - c] + r / N_BASE * (Pb[0, i_uc + c] - Pb[0, i_lc - c])))
        col.append(int(up))
        if c == 0:
            wall_nodes.update(col)
        else:
            for r in range(N_BASE):
                cells.append([prev[r], col[r], col[r + 1], prev[r + 1]])
        prev = col
    n_quad = len(cells)

    ring = Pb[nj]
    d, ridx = cKDTree(ring).query(V)
    vmap = np.empty(len(V), dtype=np.int64)
    used = np.zeros(len(V), bool)
    used[T.ravel()] = True
    for v in range(len(V)):
        if d[v] < 1e-9:
            vmap[v] = nid[nj, ridx[v]]
        elif used[v]:
            vmap[v] = node(V[v])
        else:
            vmap[v] = -1                      # construction point (arc centre)
    for tri in T:
        cells.append([int(vmap[v]) for v in tri])

    pts = np.array(pts)
    for k, poly in enumerate(cells):          # every cell CCW
        xy = pts[poly]
        a = np.sum(xy[:, 0] * np.roll(xy[:, 1], -1) - np.roll(xy[:, 0], -1) * xy[:, 1])
        if a < 0:
            cells[k] = poly[::-1]
    return pts, cells, wall_nodes, n_quad


# --------------------------------------------------------------- quality
def polygon_centroid(xy):
    x, y = xy[:, 0], xy[:, 1]
    xn, yn = np.roll(x, -1), np.roll(y, -1)
    cr = x * yn - xn * y
    a = cr.sum() / 2.0
    return np.array([((x + xn) * cr).sum(), ((y + yn) * cr).sum()]) / (6.0 * a), a


def quality(pts, cells, edges):
    """Per-cell equiangle quality (1 - equiangle skewness) and orthogonal
    quality (min cosine between each face normal and the vectors from the
    cell centre to the face centre and to the neighbour centre)."""
    cen = np.empty((len(cells), 2))
    area = np.empty(len(cells))
    eq = np.empty(len(cells))
    for c, poly in enumerate(cells):
        xy = pts[poly]
        cen[c], area[c] = polygon_centroid(xy)
        n = len(poly)
        th_e = 180.0 * (n - 2) / n
        a = xy - np.roll(xy, 1, axis=0)
        b = np.roll(xy, -1, axis=0) - xy
        ang = 180.0 - np.degrees(np.arctan2(a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0],
                                            (a * b).sum(axis=1)))
        eq[c] = 1.0 - max((ang.max() - th_e) / (180.0 - th_e), (th_e - ang.min()) / th_e)
    oq = np.ones(len(cells))
    for (a, b), lst in edges.items():
        fc = 0.5 * (pts[a] + pts[b])
        for c, p, q in lst:
            e = pts[q] - pts[p]
            nrm = np.array([e[1], -e[0]]) / np.linalg.norm(e)      # outward of c
            v = fc - cen[c]
            oq[c] = min(oq[c], nrm @ v / np.linalg.norm(v))
        if len(lst) == 2:
            (c0, p, q), (c1, _, _) = lst
            e = pts[q] - pts[p]
            nrm = np.array([e[1], -e[0]]) / np.linalg.norm(e)
            v = cen[c1] - cen[c0]
            cosv = nrm @ v / np.linalg.norm(v)
            oq[c0] = min(oq[c0], cosv)
            oq[c1] = min(oq[c1], cosv)
    return eq, oq, cen, area


# ------------------------------------------------------------ polyMesh I/O
HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "constant/polyMesh";
    object      {obj};
}}
"""


def edge_table(cells):
    edges = {}
    for c, poly in enumerate(cells):
        for k in range(len(poly)):
            a, b = poly[k], poly[(k + 1) % len(poly)]
            edges.setdefault((min(a, b), max(a, b)), []).append((c, a, b))
    return edges


def write_polymesh(case, pts, cells, edges, wall_nodes, x_out):
    np2 = len(pts)
    internal, patches = [], {"airfoil": [], "farfield": [], "outlet": []}
    for key, lst in edges.items():
        if len(lst) == 2:
            (c0, a, b), (c1, _, _) = sorted(lst)
            internal.append((c0, c1, a, b))
        elif len(lst) == 1:
            c0, a, b = lst[0]
            if a in wall_nodes and b in wall_nodes:
                patches["airfoil"].append((c0, a, b))
            elif abs(pts[a, 0] - x_out) < 1e-9 and abs(pts[b, 0] - x_out) < 1e-9:
                patches["outlet"].append((c0, a, b))
            else:
                patches["farfield"].append((c0, a, b))
        else:
            raise RuntimeError(f"edge {key} shared by {len(lst)} cells")
    internal.sort()

    faces, owner, neigh = [], [], []
    for c0, c1, a, b in internal:
        faces.append((a, b, b + np2, a + np2)); owner.append(c0); neigh.append(c1)
    bounds = []
    for name in ("airfoil", "farfield", "outlet"):
        start = len(faces)
        for c0, a, b in sorted(patches[name]):
            faces.append((a, b, b + np2, a + np2)); owner.append(c0)
        bounds.append((name, "wall" if name == "airfoil" else "patch", start, len(faces) - start))
    start = len(faces)
    for c, poly in enumerate(cells):
        faces.append(tuple(poly[::-1])); owner.append(c)                 # back, -z
        faces.append(tuple(v + np2 for v in poly)); owner.append(c)      # front, +z
    bounds.append(("frontAndBack", "empty", start, len(faces) - start))

    d = os.path.join(case, "constant", "polyMesh")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "points"), "w") as f:
        f.write(HEADER.format(cls="vectorField", obj="points"))
        f.write(f"\n{2*np2}\n(\n")
        for z in (0.0, DZ):
            f.writelines(f"({x:.12g} {y:.12g} {z:.12g})\n" for x, y in pts)
        f.write(")\n")
    with open(os.path.join(d, "faces"), "w") as f:
        f.write(HEADER.format(cls="faceList", obj="faces"))
        f.write(f"\n{len(faces)}\n(\n")
        f.writelines(f"{len(fc)}({' '.join(map(str, fc))})\n" for fc in faces)
        f.write(")\n")
    for name, data in (("owner", owner), ("neighbour", neigh)):
        with open(os.path.join(d, name), "w") as f:
            f.write(HEADER.format(cls="labelList", obj=name))
            f.write(f"\n{len(data)}\n(\n")
            f.writelines(f"{v}\n" for v in data)
            f.write(")\n")
    with open(os.path.join(d, "boundary"), "w") as f:
        f.write(HEADER.format(cls="polyBoundaryMesh", obj="boundary"))
        f.write(f"\n{len(bounds)}\n(\n")
        for name, typ, st, n in bounds:
            extra = "        inGroups        1(wall);\n" if typ == "wall" else ""
            f.write(f"    {name}\n    {{\n        type            {typ};\n{extra}"
                    f"        nFaces          {n};\n        startFace       {st};\n    }}\n")
        f.write(")\n")
    return {k: len(v) for k, v in patches.items()}


# ------------------------------------------------------------------- main
def main():
    lower, upper = airfoil_surfaces()
    sp = np.hypot(*np.diff(upper, axis=0).T)
    base = upper[-1, 1] - lower[0, 1]
    print(f"Freestream U = {UINF:.4f} m/s, nu = {NU:.5e} m2/s, Re = {P['Re']:.3g}")
    print(f"y+ estimate: Cf = {CF:.5f}, u_tau = {UTAU:.4f} m/s, "
          f"y(y+={P['yPlusTarget']}) = {Y_WALL:.3e} m -> cell height {2*Y_WALL:.3e} m; "
          f"x firstCellFactor {P.get('firstCellFactor', 1.0)} -> first cell height {H1:.3e} m")
    print(f"Blunt TE at x/c = {P['teCut']}: base thickness {base*1e3:.4f} mm, "
          f"{N_BASE} elements")
    print(f"Surfaces: {N_SURF} upper + {N_SURF} lower elements (sine spacing), "
          f"size {sp.min():.2e} .. {sp.max():.2e} m")

    Pb, line, wall, ratios = build_band(lower, upper)
    x_out = line[0, 0]
    print(f"Quad band: {Pb.shape[1]-1} x {N_LAY} cells, {T_BAND} m thick, wall growth "
          f"ratio {ratios[wall[0]+5]:.4f}; wake {N_WAKE} cells (cosine), first "
          f"{line[N_WAKE-1,0]-line[N_WAKE,0]:.2e} m, last {line[0,0]-line[1,0]:.3f} m; "
          f"outlet x = {x_out:.3f} m")

    V, T = triangulate(Pb[-1], x_out)
    pts, cells, wall_nodes, n_quad = assemble(Pb, wall, V, T)
    edges = edge_table(cells)
    eq, oq, cen, area = quality(pts, cells, edges)
    if area.min() <= 0:
        sys.exit(f"Inverted cells: {int((area <= 0).sum())}")
    counts = write_polymesh(CASE, pts, cells, edges, wall_nodes, x_out)
    n_tri = len(cells) - n_quad

    far = [e for e, l in edges.items() if len(l) == 1
           and not (e[0] in wall_nodes and e[1] in wall_nodes)
           and not (abs(pts[e[0], 0] - x_out) < 1e-9 and abs(pts[e[1], 0] - x_out) < 1e-9)]
    far_len = np.array([np.linalg.norm(pts[a] - pts[b]) for a, b in far])
    with open(os.path.join(CASE, "constant", "meshQuality.txt"), "w") as fh:
        rows = [("cells", f"{len(cells)}  (quads {n_quad}, triangles {n_tri})"),
                ("airfoil faces", f"{counts['airfoil']}  ({N_SURF} upper, {N_SURF} lower, "
                                  f"{N_BASE} base)"),
                ("first cell height", f"{H1:.3e} m"),
                ("equiangle quality min / 1st pct / mean",
                 f"{eq.min():.3f} / {np.percentile(eq, 1):.3f} / {eq.mean():.3f}"),
                ("  quads min, triangles min", f"{eq[:n_quad].min():.3f}, {eq[n_quad:].min():.3f}"),
                ("orthogonal quality min / 1st pct / mean",
                 f"{oq.min():.3f} / {np.percentile(oq, 1):.3f} / {oq.mean():.3f}"),
                ("  quads min, triangles min", f"{oq[:n_quad].min():.3f}, {oq[n_quad:].min():.3f}"),
                ("far-field (C-boundary) edge length min / max", f"{far_len.min():.3f} / {far_len.max():.3f} m")]
        for k, v in rows:
            fh.write(f"{k:42s} {v}\n")
            print(f"{k:42s} {v}")
    np.savez(os.path.join(CASE, "constant", "meshInfo.npz"),
             lower=lower, upper=upper, h1=H1, n_quad=n_quad, n_tri=n_tri,
             y_wall=Y_WALL, cf=CF, utau=UTAU, uinf=UINF, nu=NU,
             eq=eq, oq=oq, cen=cen)


if __name__ == "__main__":
    main()
