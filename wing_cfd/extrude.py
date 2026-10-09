"""Extrude the 2D template spanwise into an OpenFOAM polyMesh.

* layers 0 .. N_WING-1 span the wing (y = 0 .. B): airfoil interior absent,
  its edges are the 'wing' wall
* layers outboard of the tip include the airfoil-interior triangles; their
  bottom faces at y = B close the tip ('wing' wall)
* spanwise spacing is clustered to TIP_DY_FIRST on both sides of y = B so
  the flat tip cap gets y+ < 1 (the section's wall-normal clustering is kept
  outboard: relaxing it makes thin wedge cells ~90 deg non-orthogonal)
* every spanwise plane is morphed so the section matches the local chord of
  the sine leading edge (S1223 scaled about the fixed TE); the morph decays
  to zero 0.05 -> 0.5 m from the wall
Patches: wing (wall), root (symmetryPlane y=0), side (y=ymax), farfield, outlet
"""
import math
import os

import numpy as np

import params as P


def span_planes():
    """Spanwise planes: uniform DY_WING on the wing, clustered geometrically
    to TIP_DY_FIRST on both sides of the tip plane y = B, then growing by
    TIP_RATIO outboard up to TIP_DY_MAX."""
    clus = []
    h = P.TIP_DY_FIRST
    while h < P.DY_WING:
        clus.append(h)
        h *= P.TIP_RATIO
    n_uni = int(math.ceil((P.B - sum(clus)) / P.DY_WING))
    dy_uni = (P.B - sum(clus)) / n_uni
    dys = [dy_uni] * n_uni + clus[::-1]
    h = P.TIP_DY_FIRST
    out = []
    while h <= P.TIP_DY_MAX:
        out.append(h)
        h *= P.TIP_RATIO
    y = np.concatenate([[0.0], np.cumsum(dys)])
    y[-1] = P.B
    y = np.concatenate([y, P.B + np.cumsum(out)])
    return y, len(dys)


def morph_weight(d, d0=P.QUAD_THICK, d1=0.5):
    t = np.clip((d - d0) / (d1 - d0), 0, 1)
    return 1 - t * t * (3 - 2 * t)


def header(cls, obj, note=None):
    s = ("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
         f"    class       {cls};\n")
    if note:
        s += f'    note        "{note}";\n'
    s += f"    location    \"constant/polyMesh\";\n    object      {obj};\n}}\n\n"
    return s


def newell(pts):
    """Area vectors of polygons pts[nf, nv, 3]."""
    a = pts
    b = np.roll(pts, -1, axis=1)
    return 0.5 * np.cross(a, b).sum(axis=1)


def main(case="case"):
    d = np.load("template2d.npz")
    nodes2 = d["nodes"]
    quads, ftri, stri = d["quads"], d["ftri"], d["stri"]
    n2 = len(nodes2)
    yp, n_wing = span_planes()
    NL = len(yp) - 1

    # ---- 2D cells: fluid (quads, tris) then solid (foil tris)
    cells2 = [list(c) for c in quads] + [list(c) for c in ftri] + \
             [list(c) for c in stri]
    nfl = len(quads) + len(ftri)
    nsol = len(stri)
    nc2 = nfl + nsol
    # active per layer
    active = np.ones((NL, nc2), bool)
    active[:n_wing, nfl:] = False
    cnt = active.sum(1)
    off = np.concatenate([[0], np.cumsum(cnt)])
    loc = np.cumsum(active, axis=1) - 1            # local index in layer
    id3 = np.where(active, off[:-1, None] + loc, -1)
    ncells = int(off[-1])

    # ---- points (morphed per plane)
    xte = np.array([P.C, 0.0])
    w = morph_weight(d["dist"])[:, None]
    pts = np.zeros((NL + 1, n2, 3))
    for j, y in enumerate(yp):
        s = P.local_chord(y) / P.C
        xz = nodes2 + w * ((xte + (nodes2 - xte) * s) - nodes2)
        pts[j, :, 0] = xz[:, 0]
        pts[j, :, 1] = y
        pts[j, :, 2] = xz[:, 1]
    pid = lambda node, j: j * n2 + node

    # ---- 2D edges -> cells
    e_a, e_b, e_c = [], [], []
    for ci, c in enumerate(cells2):
        for k in range(len(c)):
            e_a.append(c[k]); e_b.append(c[(k + 1) % len(c)]); e_c.append(ci)
    e_a, e_b, e_c = map(np.array, (e_a, e_b, e_c))
    key = np.minimum(e_a, e_b).astype(np.int64) * n2 + np.maximum(e_a, e_b)
    order = np.argsort(key, kind="stable")
    key_s = key[order]
    uk, first, counts = np.unique(key_s, return_index=True, return_counts=True)
    assert counts.max() <= 2
    ia = order[first]
    c1 = e_c[ia]
    c2 = np.where(counts == 2, e_c[order[np.minimum(first + 1, len(order) - 1)]], -1)
    ea, eb = e_a[ia], e_b[ia]
    # boundary tags of 2D boundary edges
    def keys(E):
        return np.minimum(E[:, 0], E[:, 1]).astype(np.int64) * n2 + np.maximum(E[:, 0], E[:, 1])
    tag = np.full(len(uk), -1)
    for t, name in enumerate(("outlet", "farfield")):
        tag[np.isin(uk, keys(d[name]))] = t
    is_wall = np.isin(uk, keys(d["wall"]))

    patches = ["wing", "root", "side", "farfield", "outlet"]
    PW, PR, PS, PF, PO = range(5)
    tag2patch = {0: PO, 1: PF}

    faces, owner, nbr, patch = [], [], [], []

    # ---- side faces
    for j in range(NL):
        a1 = id3[j][c1]
        a2 = np.where(c2 >= 0, id3[j][np.maximum(c2, 0)], -1)
        quad = np.column_stack([pid(ea, j), pid(eb, j), pid(eb, j + 1), pid(ea, j + 1)])
        both = (a1 >= 0) & (a2 >= 0)
        faces.append(quad[both]); owner.append(np.minimum(a1, a2)[both])
        nbr.append(np.maximum(a1, a2)[both]); patch.append(np.full(both.sum(), -1))
        one = (a1 >= 0) ^ (a2 >= 0)
        own = np.where(a1 >= 0, a1, a2)
        pt = np.where(is_wall, PW, np.vectorize(lambda t: tag2patch.get(t, -9))(tag))
        sel = one
        assert (pt[sel] >= 0).all(), "untagged boundary edge"
        faces.append(quad[sel]); owner.append(own[sel]); nbr.append(np.full(sel.sum(), -1))
        patch.append(pt[sel])

    # ---- cap faces (tri caps padded with -1)
    cap = np.full((nc2, 4), -1)
    for ci, c in enumerate(cells2):
        cap[ci, :len(c)] = c
    for j in range(NL + 1):
        below = id3[j - 1] if j > 0 else np.full(nc2, -1)
        above = id3[j] if j < NL else np.full(nc2, -1)
        f = np.where(cap >= 0, pid(cap, j), -1)
        both = (below >= 0) & (above >= 0)
        faces.append(f[both]); owner.append(below[both]); nbr.append(above[both])
        patch.append(np.full(both.sum(), -1))
        if j == 0:
            sel, own, pt = above >= 0, above, PR
        elif j == NL:
            sel, own, pt = below >= 0, below, PS
        else:
            sel = (above >= 0) & (below < 0)
            own, pt = above, PW
            assert not ((below >= 0) & (above < 0)).any()
        faces.append(f[sel]); owner.append(own[sel]); nbr.append(np.full(sel.sum(), -1))
        patch.append(np.full(sel.sum(), pt))

    F = np.vstack([np.pad(f, ((0, 0), (0, 4 - f.shape[1])), constant_values=-1)
                   if f.shape[1] < 4 else f for f in faces])
    own = np.concatenate(owner)
    nb = np.concatenate(nbr)
    pa = np.concatenate(patch)

    # ---- compact points
    used = np.unique(F[F >= 0])
    newid = np.full((NL + 1) * n2, -1)
    newid[used] = np.arange(len(used))
    allp = pts.reshape(-1, 3)
    points = allp[used]

    # ---- orientation: normal must point out of the owner
    # cell centres = mean of the face centres of its faces (approx.)
    nv = (F >= 0).sum(1)
    Fp = allp[np.where(F >= 0, F, F[:, :1])]          # pad with first vertex
    fc = Fp.sum(1) / 4.0
    fc = np.where(nv[:, None] == 3, (Fp[:, :3].sum(1)) / 3.0, fc)
    cc = np.zeros((ncells, 3)); cn = np.zeros(ncells)
    np.add.at(cc, own, fc); np.add.at(cn, own, 1)
    m = nb >= 0
    np.add.at(cc, nb[m], fc[m]); np.add.at(cn, nb[m], 1)
    cc /= cn[:, None]
    tri = nv == 3
    area = np.zeros((len(F), 3))
    area[~tri] = newell(Fp[~tri])
    area[tri] = newell(Fp[tri][:, :3])
    flip = (area * (fc - cc[own])).sum(1) < 0
    Fr = F.copy()
    Fr[flip & ~tri] = F[flip & ~tri][:, ::-1]
    Fr[flip & tri, :3] = F[flip & tri][:, 2::-1]
    F = Fr
    print(f"flipped {flip.sum()} faces")

    # ---- order: internal (owner, neighbour), then patches
    internal = nb >= 0
    oi = np.where(internal)[0]
    oi = oi[np.lexsort((nb[oi], own[oi]))]
    ob = [np.where(pa == p)[0] for p in range(len(patches))]
    ordr = np.concatenate([oi] + ob)
    F, own, nb = F[ordr], own[ordr], nb[ordr]
    nint = len(oi)

    out = os.path.join(case, "constant", "polyMesh")
    os.makedirs(out, exist_ok=True)
    note = f"nPoints:{len(points)} nCells:{ncells} nFaces:{len(F)} nInternalFaces:{nint}"
    with open(os.path.join(out, "points"), "w") as fh:
        fh.write(header("vectorField", "points"))
        fh.write(f"{len(points)}\n(\n")
        np.savetxt(fh, points, fmt="(%.9g %.9g %.9g)")
        fh.write(")\n")
    Fn = np.where(F >= 0, newid[np.maximum(F, 0)], -1)
    with open(os.path.join(out, "faces"), "w") as fh:
        fh.write(header("faceList", "faces"))
        fh.write(f"{len(Fn)}\n(\n")
        q = Fn[:, 3] >= 0
        lines = np.empty(len(Fn), dtype=object)
        lines[q] = [f"4({a} {b} {c} {e})" for a, b, c, e in Fn[q]]
        lines[~q] = [f"3({a} {b} {c})" for a, b, c in Fn[~q, :3]]
        fh.write("\n".join(lines))
        fh.write("\n)\n")
    for name, arr in (("owner", own), ("neighbour", nb[:nint])):
        with open(os.path.join(out, name), "w") as fh:
            fh.write(header("labelList", name, note))
            fh.write(f"{len(arr)}\n(\n")
            np.savetxt(fh, arr, fmt="%d")
            fh.write(")\n")
    types = {"wing": "wall", "root": "symmetryPlane", "side": "patch",
             "farfield": "patch", "outlet": "patch"}
    start = nint
    with open(os.path.join(out, "boundary"), "w") as fh:
        fh.write(header("polyBoundaryMesh", "boundary"))
        fh.write(f"{len(patches)}\n(\n")
        for p, name in enumerate(patches):
            n = len(ob[p])
            ig = "        inGroups        1(wall);\n" if types[name] == "wall" else ""
            fh.write(f"    {name}\n    {{\n        type            {types[name]};\n{ig}"
                     f"        nFaces          {n};\n        startFace       {start};\n    }}\n")
            start += n
        fh.write(")\n")
    print(note)
    print("span planes:", len(yp), "wing layers:", n_wing, "ymax:", yp[-1])
    np.save(os.path.join(case, "span_planes.npy"), yp)


if __name__ == "__main__":
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else "case")
