"""Printed-mass model and print-bed segmentation."""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .hydrostatics import area_props, volume_props
from .params import Config, Design


@dataclass
class Item:
    name: str
    mass: float
    pos: np.ndarray  # body frame CG


@dataclass
class MassBudget:
    items: list[Item] = field(default_factory=list)

    def add(self, name: str, mass: float, pos) -> None:
        self.items.append(Item(name, float(mass), np.asarray(pos, float)))

    @property
    def total(self) -> float:
        return sum(i.mass for i in self.items)

    @property
    def cg(self) -> np.ndarray:
        return sum(i.mass * i.pos for i in self.items) / self.total

    def by_name(self) -> dict[str, float]:
        out: dict[str, float] = {}
        for i in self.items:
            out[i.name] = out.get(i.name, 0.0) + i.mass
        return out


def print_pieces(d: Design, bed) -> int:
    """Fewest pieces one demihull splits into to fit the printer bed.

    The hull is cut along x into segments. The segment length uses one bed
    axis and the section (beam x depth) uses the other two, split further if
    it does not fit.
    """
    best = math.inf
    for k in range(3):
        lb = bed[k]
        c1, c2 = [bed[j] for j in range(3) if j != k]
        n_len = math.ceil(d.loa / lb)
        n_sec = min(
            math.ceil(d.beam / c1) * math.ceil(d.depth / c2),
            math.ceil(d.beam / c2) * math.ceil(d.depth / c1),
        )
        best = min(best, n_len * n_sec)
    return int(best)


@dataclass
class HullMass:
    mass: float
    cg: np.ndarray        # demihull centred on y=0
    shell: float
    ribs: float
    sealer: float
    joints: float
    pieces: int
    area: float
    enclosed_volume: float


def demihull_mass(tri: np.ndarray, d: Design, cfg: Config) -> HullMass:
    """Shell + ribs/infill + sealer + joint flanges for one demihull.

    shell = surface area x wall thickness x filament density
    ribs  = enclosed volume x rib volume fraction x filament density
    The printed part (shell + ribs + joints) is scaled by `mass_calibration`,
    which should come from weighing a test print.
    """
    A, ca = area_props(tri)
    V, cv = volume_props(tri)
    shell = A * cfg.wall_thickness * cfg.fil_density
    ribs = V * cfg.rib_vol_frac * cfg.fil_density
    pieces = print_pieces(d, cfg.bed)
    girth = 2 * d.depth + 2 * d.beam
    joints = (pieces - 1) * girth * cfg.joint_kg_per_m
    k = cfg.mass_calibration
    shell, ribs, joints = k * shell, k * ribs, k * joints
    sealer = A * cfg.sealer_kg_per_m2
    total = shell + ribs + joints + sealer
    cg = ((shell + sealer + joints) * ca + ribs * cv) / total
    return HullMass(total, cg, shell, ribs, sealer, joints, pieces, A, V)


def calibration_factor(measured_kg: float, predicted_printed_kg: float) -> float:
    """Factor for Config.mass_calibration from a weighed test print (before sealing)."""
    return measured_kg / predicted_printed_kg
