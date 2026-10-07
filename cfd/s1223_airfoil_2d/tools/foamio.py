"""Minimal readers for an ASCII 2D OpenFOAM case (one cell deep)."""
import os
import re

import numpy as np


def _body(path):
    with open(path) as fh:
        s = fh.read()
    s = re.sub(r"/\*.*?\*/", "", s, flags=re.S)
    s = re.sub(r"//[^\n]*", "", s)
    return s[s.index("}", s.index("FoamFile")) + 1:]


def _list(s, start=0):
    """Parse an OpenFOAM list `N ( ... )` starting at or after `start`."""
    m = re.compile(r"(\d+)\s*\(").search(s, start)
    n, i = int(m.group(1)), m.end()
    depth, j = 1, i
    while depth:
        c = s[j]
        depth += (c == "(") - (c == ")")
        j += 1
    return n, s[i:j - 1], j


def read_points(case):
    n, body, _ = _list(_body(os.path.join(case, "constant/polyMesh/points")))
    return np.array([[float(v) for v in t.split()] for t in re.findall(r"\(([^()]*)\)", body)])


def read_faces(case):
    n, body, _ = _list(_body(os.path.join(case, "constant/polyMesh/faces")))
    return [[int(v) for v in t.split()] for t in re.findall(r"\d+\(([^()]*)\)", body)]


def read_labels(case, name):
    n, body, _ = _list(_body(os.path.join(case, "constant/polyMesh", name)))
    return np.array(body.split(), dtype=np.int64)


def read_boundary(case):
    s = _body(os.path.join(case, "constant/polyMesh/boundary"))
    out = {}
    for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", s):
        d = dict(re.findall(r"(\w+)\s+([^;]+);", m.group(2)))
        out[m.group(1)] = dict(type=d["type"].strip(), nFaces=int(d["nFaces"]),
                               startFace=int(d["startFace"]))
    return out


class Mesh2D:
    """Cell polygons (z = 0 plane), centres and wall-patch faces."""

    def __init__(self, case):
        pts = read_points(case)
        faces = read_faces(case)
        owner = read_labels(case, "owner")
        bnd = read_boundary(case)
        self.ncells = int(owner.max()) + 1
        zmin = pts[:, 2].min()
        poly = [None] * self.ncells
        fb = bnd["frontAndBack"]
        for f in range(fb["startFace"], fb["startFace"] + fb["nFaces"]):
            nodes = faces[f]
            if np.all(np.abs(pts[nodes, 2] - zmin) < 1e-12):
                poly[owner[f]] = pts[nodes, :2]
        self.polys = poly
        self.centres = np.array([p.mean(axis=0) for p in poly])
        self.pts, self.faces, self.owner, self.boundary = pts, faces, owner, bnd

    def patch(self, name):
        """Face centres (x, y), outward normals and owner cells of a 2D patch."""
        b = self.boundary[name]
        idx = range(b["startFace"], b["startFace"] + b["nFaces"])
        xy, nrm, own = [], [], []
        for f in idx:
            p = self.pts[self.faces[f]]
            xy.append(p[:, :2].mean(axis=0))
            n = np.cross(p[1] - p[0], p[2] - p[1])[:2]
            nrm.append(n / np.linalg.norm(n))
            own.append(self.owner[f])
        return np.array(xy), np.array(nrm), np.array(own)


def read_field(case, time, name, patch=None):
    """internalField (or a patch's boundary values) of an ASCII field."""
    path = os.path.join(case, str(time), name)
    with open(path) as fh:
        vector = "vector" in re.search(r"class\s+(\w+)", fh.read(2000)).group(1).lower()
    s = _body(path)
    if patch is None:
        m = re.search(r"internalField\s+(uniform|nonuniform)", s)
        start = m.end()
    else:
        m = re.search(rf"\b{patch}\s*\{{", s[s.index("boundaryField"):])
        start = s.index("boundaryField") + m.end()
        m = re.search(r"value\s+(uniform|nonuniform)", s[start:])
        start += m.end()
    if m.group(1) == "nonuniform":
        n, body, _ = _list(s, start)
        if vector:
            return np.array([[float(v) for v in t.split()] for t in re.findall(r"\(([^()]*)\)", body)])
        return np.array(body.split(), dtype=float)
    raw = s[start:s.index(";", start)].strip()
    return np.array([float(v) for v in raw.strip("()").split()])
