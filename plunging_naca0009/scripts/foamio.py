"""
Minimal readers for ascii OpenFOAM meshes and fields (2-D, one cell deep),
used by post.py and animate.py.
"""
import os
import re

import matplotlib.tri as mtri
import numpy as np


# =================================================================== readers
def _body(path):
    with open(path) as f:
        txt = f.read()
    i = txt.find("FoamFile")
    if i >= 0:
        txt = txt[txt.index("}", i) + 1:]
    return re.sub(r"//.*", "", txt)


def _list_block(txt, start=0):
    """Parse 'N ( ... )' starting at the first count after `start`."""
    m = re.compile(r"(\d+)\s*\(").search(txt, start)
    n, i = int(m.group(1)), m.end()
    depth, j = 1, i
    while depth:
        c = txt[j]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        j += 1
    return n, txt[i:j - 1], j


def read_points(path):
    n, s, _ = _list_block(_body(path))
    return np.array(s.replace("(", " ").replace(")", " ").split(), float).reshape(n, 3)


def read_labels(path):
    n, s, _ = _list_block(_body(path))
    return np.array(s.split(), int)


def read_faces(path):
    txt = _body(path)
    with open(path) as f:
        head = f.read(2000)
    if "faceCompactList" in head:
        n1, s1, j = _list_block(txt)
        n2, s2, _ = _list_block(txt, j)
        off, lab = np.array(s1.split(), int), np.array(s2.split(), int)
        return [lab[off[k]:off[k + 1]] for k in range(len(off) - 1)]
    _, s, _ = _list_block(txt)
    return [np.array(g.split(), int) for g in re.findall(r"\d+\s*\(([^)]*)\)", s)]


def read_boundary(path):
    txt = _body(path)
    out = {}
    for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", txt):
        d = dict(re.findall(r"(\w+)\s+([^;]+);", m.group(2)))
        if "nFaces" in d:
            out[m.group(1)] = (int(d["startFace"]), int(d["nFaces"]))
    return out


def _parse_values(txt, ncomp):
    s = txt.replace("(", " ").replace(")", " ")
    a = np.array(s.split(), float)
    return a.reshape(-1, ncomp) if ncomp > 1 else a


def read_field(path, ncells):
    """Return (internal values, {patch: values}) for an ascii vol field."""
    txt = _body(path)
    vec = "volVectorField" in open(path).read(1500)
    nc = 3 if vec else 1
    i = txt.index("internalField")
    seg = txt[i: i + 200]
    if "nonuniform" in seg:
        _, s, _ = _list_block(txt, i)
        internal = _parse_values(s, nc)
    else:
        v = re.search(r"uniform\s+(\([^)]*\)|[-+.\deE]+)", seg).group(1)
        internal = np.tile(_parse_values(v, nc), (ncells, 1)) if vec else np.full(ncells, float(v))
    patches = {}
    bf = txt.index("boundaryField")
    for m in re.finditer(r"\n\s{4}(\w+)\s*\n\s{4}\{", txt[bf:]):
        name = m.group(1)
        k0 = bf + m.end()
        k1 = txt.index("\n    }", k0)
        blk = txt[k0:k1]
        vm = re.search(r"value\s+nonuniform", blk)
        if vm:
            _, s, _ = _list_block(blk, vm.end())
            patches[name] = _parse_values(s, nc)
        else:
            vm = re.search(r"value\s+uniform\s+(\([^)]*\)|[-+.\deE]+)", blk)
            if vm:
                patches[name] = _parse_values(vm.group(1), nc)
    return internal, patches


# =================================================================== mesh
class Mesh2D:
    def __init__(self, case):
        pm = os.path.join(case, "constant", "polyMesh")
        self.points = read_points(os.path.join(pm, "points"))
        self.faces = read_faces(os.path.join(pm, "faces"))
        self.owner = read_labels(os.path.join(pm, "owner"))
        self.boundary = read_boundary(os.path.join(pm, "boundary"))
        self.ncells = self.owner.max() + 1
        nb = os.path.join(pm, "neighbour")
        if os.path.exists(nb):
            self.ncells = max(self.ncells, read_labels(nb).max() + 1)
        s, n = self.boundary["frontAndBack"]
        z0 = self.points[:, 2].min()
        polys, cells = [], []
        for f in range(s, s + n):
            fv = self.faces[f]
            if np.all(np.abs(self.points[fv, 2] - z0) < 1e-9):
                polys.append(fv)
                cells.append(self.owner[f])
        self.poly_cell = np.array(cells)
        used = np.unique(np.concatenate(polys))
        remap = -np.ones(len(self.points), int)
        remap[used] = np.arange(len(used))
        self.xy = self.points[used, :2]
        self.polys = [remap[p] for p in polys]
        tris, tcell = [], []
        for p, c in zip(self.polys, self.poly_cell):
            for k in range(1, len(p) - 1):
                tris.append((p[0], p[k], p[k + 1]))
                tcell.append(c)
        self.tri = mtri.Triangulation(self.xy[:, 0], self.xy[:, 1], np.array(tris))
        # cell centres & cell->point averaging operator
        cc = np.zeros((self.ncells, 2))
        cnt = np.zeros(self.ncells)
        for p, c in zip(self.polys, self.poly_cell):
            cc[c] += self.xy[p].mean(axis=0)
            cnt[c] += 1
        self.cc = cc / np.maximum(cnt, 1)[:, None]
        rows = np.concatenate(self.polys)
        cols = np.concatenate([[c] * len(p) for p, c in zip(self.polys, self.poly_cell)])
        self._pr, self._pc = rows, cols
        self._pn = np.bincount(rows, minlength=len(self.xy))

    def to_points(self, cellvals):
        return np.bincount(self._pr, weights=cellvals[self._pc], minlength=len(self.xy)) / self._pn

    def patch(self, name):
        """Wall-face centres (x, y), outward unit normals (into the wall), owner cells."""
        s, n = self.boundary[name]
        fc, nrm, own = [], [], []
        for f in range(s, s + n):
            fv = self.points[self.faces[f]]
            c = fv.mean(axis=0)
            nv = np.cross(fv[1] - fv[0], fv[2] - fv[0])
            fc.append(c[:2])
            nrm.append(nv[:2] / np.linalg.norm(nv[:2]))
            own.append(self.owner[f])
        return np.array(fc), np.array(nrm), np.array(own)


def latest_time(case):
    ts = []
    for d in os.listdir(case):
        try:
            ts.append((float(d), d))
        except ValueError:
            pass
    ts = [t for t in ts if t[0] > 0]
    return max(ts)[1] if ts else None


def time_dirs(case):
    """All time directories > 0 as (value, name), sorted."""
    ts = []
    for d in os.listdir(case):
        try:
            v = float(d)
        except ValueError:
            continue
        if v > 0 and os.path.isdir(os.path.join(case, d)):
            ts.append((v, d))
    return sorted(ts)
