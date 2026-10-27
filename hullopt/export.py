"""STL export of a design: front pod (with windows, thruster and pylon) and the
whole boat (hulls, skids, pods, thrusters, struts) for viewing.

Viewing meshes only: the pod is its outer skin with the side windows cut out,
so the thruster shows through. Body frame, metres.

  python -m hullopt.export --design results/stage2/designs/finalist_1.json --out finalist_1
  -> finalist_1_pod.stl, finalist_1_assembly.stl
"""
from __future__ import annotations

import argparse
import json
import math

import numpy as np

from .geometry import Mesh, catamaran_mesh, hull_centre_y, keel_z, loft, write_stl
from .layout import Pod, Strut, Thruster
from .params import Config, Design
from .stage1 import evaluate


def _ellipse(n: int, a: float, b: float):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return a * np.cos(t), b * np.sin(t)   # CCW in (y, z)


def _to_world(m: Mesh, R: np.ndarray, t) -> Mesh:
    return Mesh(m.vertices @ R.T + np.asarray(t, float), m.faces)


def _bay_len(pod: Pod) -> float:
    """Fore-aft length of the thruster bay (inside the walls)."""
    return pod.bay_side_area / (2 * (pod.half_h - pod.wall))


def pod_profile(pod: Pod, x: np.ndarray) -> np.ndarray:
    """Section scale (0..1) along the pod: tail fairing, full-section bay, elliptic nose."""
    nose0 = pod.x1 - pod.nose_len
    tail = nose0 - pod.x0 - 2 * pod.wall - _bay_len(pod)   # = design pod_tail_len
    s = np.ones_like(x)
    if tail > 1e-6:
        m = x < pod.x0 + tail
        u = (x[m] - pod.x0) / tail
        s[m] = 0.6 + 0.4 * np.sqrt(np.clip(1 - (1 - u) ** 2, 0, 1))
    m = x > nose0
    u = (x[m] - nose0) / pod.nose_len
    s[m] = np.sqrt(np.clip(1 - u**2, 0, 1))
    return np.maximum(s, 0.03)


def pod_mesh(pod: Pod, n_x: int = 60, n_c: int = 48, windows: bool = True) -> Mesh:
    """Pod outer skin; side windows (area pod.window_area each) cut out of the bay."""
    x = np.linspace(pod.x0, pod.x1, n_x)
    s = pod_profile(pod, x)
    yy, zz = _ellipse(n_c, 1.0, 1.0)
    Y = pod.y + s[:, None] * pod.half_w * yy[None]
    Z = pod.z + s[:, None] * pod.half_h * zz[None]
    m = loft(x, Y, Z)
    if not windows:
        return m
    # window: centred on the bay, on each side, rectangle in the side projection
    bay_len = _bay_len(pod)
    nose0 = pod.x1 - pod.nose_len
    xc = nose0 - pod.wall - bay_len / 2
    wl = 0.8 * bay_len
    wh = min(pod.window_area / wl, 1.6 * pod.half_h)
    c = m.triangles.mean(1)
    side = np.abs(c[:, 1] - pod.y) > 0.5 * pod.half_w
    inside = (np.abs(c[:, 0] - xc) < wl / 2) & (np.abs(c[:, 2] - pod.z) < wh / 2) & side
    return Mesh(m.vertices, m.faces[~inside])


def thruster_mesh(t: Thruster, n: int = 32) -> Mesh:
    """T200 as a cylinder (duct outer radius, overall length) on its axis."""
    yy, zz = _ellipse(n, t.radius, t.radius)
    x = np.array([-t.length / 2, t.length / 2])
    m = loft(x, np.tile(yy, (2, 1)), np.tile(zz, (2, 1)))
    c, s_ = math.cos(t.heading), math.sin(t.heading)
    R = np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])
    return _to_world(m, R, t.centre)


def foil_mesh(s: Strut, n: int = 24) -> Mesh:
    """Vertical NACA 00xx strut from z_bot to z_top, chord along x, centred on (x, y)."""
    xc = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n // 2 + 1)))
    yt = 5 * s.tc * (0.2969 * np.sqrt(xc) - 0.1260 * xc - 0.3516 * xc**2 + 0.2843 * xc**3 - 0.1036 * xc**4)
    # upper surface LE->TE then lower TE->LE is clockwise in (x, y); reverse it to CCW
    px = (np.concatenate([xc, xc[-2:0:-1]]) - 0.5)[::-1]
    py = np.concatenate([yt, -yt[-2:0:-1]])[::-1]
    a = np.array([s.z_bot, s.z_top])
    m = loft(a, np.tile(px * s.chord, (2, 1)), np.tile(py * s.chord, (2, 1)))
    # loft frame (a, b, c) -> world (x, y, z) = (b, c, a): a proper rotation
    R = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float)
    return _to_world(m, R, [s.x, s.y, 0.0])


def skid_mesh(skid, y: float) -> Mesh:
    """Skid rail: rectangle-ish section from its flat bottom up to the keel line."""
    xs, h, w = skid.xs, skid.heights, skid.width
    zb = skid.z_bottom
    n = len(xs)
    ys = np.array([0, w / 2, w / 2, 0, -w / 2, -w / 2]) + y
    Y = np.tile(ys, (n, 1))
    Z = np.stack([np.full(n, zb), np.full(n, zb), zb + h, zb + h, zb + h, np.full(n, zb)], 1)
    return loft(xs, Y, Z)


def build(d: Design, cfg: Config):
    row, det = evaluate(d, cfg, detail=True)
    if det is None:
        raise ValueError(f"design does not float: {row.get('fail_reason')}")
    th = {t.name: t for t in det.thrusters}
    st = {s.name: s for s in det.struts}
    pod = det.pods[0]
    pod_parts = Mesh.merge(pod_mesh(pod), thruster_mesh(th["front_port"]), foil_mesh(st["front_port"]))
    yc = hull_centre_y(d)
    parts = [catamaran_mesh(d, cfg.n_stations)]
    parts += [skid_mesh(det.skid, yc), skid_mesh(det.skid, -yc)]
    parts += [pod_mesh(p) for p in det.pods]
    parts += [thruster_mesh(t) for t in det.thrusters]
    parts += [foil_mesh(s) for s in det.struts]
    return row, pod_parts, Mesh.merge(*parts)


def export(d: Design, cfg: Config, prefix: str) -> list[str]:
    _, pod, boat = build(d, cfg)
    paths = [f"{prefix}_pod.stl", f"{prefix}_assembly.stl"]
    write_stl(paths[0], pod.triangles, "front pod")
    write_stl(paths[1], boat.triangles, "assembly")
    return paths


def main(argv=None):
    p = argparse.ArgumentParser(prog="hullopt.export", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--design", help="design JSON (defaults if omitted)")
    p.add_argument("--config")
    p.add_argument("--out", default="design", help="output prefix")
    a = p.parse_args(argv)
    d = Design()
    if a.design:
        with open(a.design) as f:
            d = Design.from_dict(json.load(f))
    for path in export(d, Config.load(a.config), a.out):
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
