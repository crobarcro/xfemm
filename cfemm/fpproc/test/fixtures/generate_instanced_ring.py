#!/usr/bin/env python3
"""Generate the Milestone D rotational template and full-ring control fixtures.

The template is one 45-degree annular sector whose two radial edges share a
periodic boundary property. The control is the equivalent full ring with
Dirichlet inner and outer boundaries. Both use a uniform applied current
density so the field is axisymmetric and the two models are directly
comparable.

This script has no external dependencies and no moving revision: it is the
recorded generator for d_ring_tile.fem and d_ring_control.fem.
"""

import math
import os

SECTORS = 8
ANGLE = 360.0 / SECTORS
INNER_RADIUS = 0.05
OUTER_RADIUS = 0.10
CURRENT_DENSITY = 1.0
ARC_SEGMENT_ANGLE = 5.0


def coord(radius, degrees):
    radians = math.radians(degrees)
    return radius * math.cos(radians), radius * math.sin(radians)


HEADER = """[Format]      =  4.0
[Frequency]   =  0
[Precision]   =  1e-8
[MinAngle]    =  30
[Depth]       =  1
[LengthUnits] =  meters
[ProblemType] =  planar
[Coordinates] =  cartesian
[ACSolver]    =  0
[PrevSoln]    = ""
[PrevType]    =  0
[Comment]     =  "{comment}"
"""


def boundary(name, bdry_type):
    return f"""  <BeginBdry>
    <BdryName> = "{name}"
    <BdryType> = {bdry_type}
    <A_0> = 0
    <A_1> = 0
    <A_2> = 0
    <Phi> = 0
    <c0> = 0
    <c0i> = 0
    <c1> = 0
    <c1i> = 0
    <Mu_ssd> = 0
    <Sigma_ssd> = 0
    <innerangle> = 0
    <outerangle> = 0
  <EndBdry>
"""


BLOCK = """  <BeginBlock>
    <BlockName> = "Air"
    <Mu_x> = 1
    <Mu_y> = 1
    <H_c> = 0
    <H_cAngle> = 0
    <J_re> = {current}
    <J_im> = 0
    <Sigma> = 0
    <d_lam> = 0
    <Phi_h> = 0
    <Phi_hx> = 0
    <Phi_hy> = 0
    <LamType> = 0
    <LamFill> = 1
    <NStrands> = 0
    <WireD> = 0
    <BHPoints> = 0
  <EndBlock>
"""


def number(value):
    return f"{value:.17g}"


def write_tile(path):
    lines = [HEADER.format(comment="Milestone D rotational template (one 45 degree sector)")]
    lines.append("[PointProps]  = 0\n")
    lines.append("[BdryProps]   = 2\n")
    lines.append(boundary("Dirichlet", 0))
    lines.append(boundary("Periodic", 4))
    lines.append("[BlockProps]  = 1\n")
    lines.append(BLOCK.format(current=number(CURRENT_DENSITY)))
    lines.append("[CircuitProps]  = 0\n")

    inner0 = coord(INNER_RADIUS, 0.0)
    inner1 = coord(INNER_RADIUS, ANGLE)
    outer0 = coord(OUTER_RADIUS, 0.0)
    outer1 = coord(OUTER_RADIUS, ANGLE)
    nodes = [inner0, inner1, outer0, outer1]
    lines.append(f"[NumPoints] = {len(nodes)}\n")
    for x, y in nodes:
        lines.append(f"{number(x)} {number(y)} 0 0\n")

    # Radial edges carry the periodic boundary (1-based property 2).
    lines.append("[NumSegments] = 2\n")
    lines.append("0 2 -1 2 0 0\n")
    lines.append("1 3 -1 2 0 0\n")

    # Inner and outer arcs carry the Dirichlet boundary (1-based property 1).
    lines.append("[NumArcSegments] = 2\n")
    lines.append(f"0 1 {number(ANGLE)} {number(ARC_SEGMENT_ANGLE)} 1 0 0\n")
    lines.append(f"2 3 {number(ANGLE)} {number(ARC_SEGMENT_ANGLE)} 1 0 0\n")

    lines.append("[NumHoles] = 0\n")
    lines.append("[NumBlockLabels] = 1\n")
    label = coord((INNER_RADIUS + OUTER_RADIUS) / 2.0, ANGLE / 2.0)
    lines.append(f"{number(label[0])} {number(label[1])} 1 0.001 0 0 0 1 0\n")

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("".join(lines))


def write_control(path):
    lines = [HEADER.format(comment="Milestone D conventional full-ring control")]
    lines.append("[PointProps]  = 0\n")
    lines.append("[BdryProps]   = 1\n")
    lines.append(boundary("Dirichlet", 0))
    lines.append("[BlockProps]  = 1\n")
    lines.append(BLOCK.format(current=number(CURRENT_DENSITY)))
    lines.append("[CircuitProps]  = 0\n")

    inner = [coord(INNER_RADIUS, k * ANGLE) for k in range(SECTORS)]
    outer = [coord(OUTER_RADIUS, k * ANGLE) for k in range(SECTORS)]
    lines.append(f"[NumPoints] = {2 * SECTORS}\n")
    for x, y in inner + outer:
        lines.append(f"{number(x)} {number(y)} 0 0\n")

    lines.append("[NumSegments] = 0\n")
    lines.append(f"[NumArcSegments] = {2 * SECTORS}\n")
    for k in range(SECTORS):
        n0, n1 = k, (k + 1) % SECTORS
        lines.append(f"{n0} {n1} {number(ANGLE)} {number(ARC_SEGMENT_ANGLE)} 1 0 0\n")
    for k in range(SECTORS):
        n0, n1 = SECTORS + k, SECTORS + (k + 1) % SECTORS
        lines.append(f"{n0} {n1} {number(ANGLE)} {number(ARC_SEGMENT_ANGLE)} 1 0 0\n")

    # The enclosed inner disk is declared as a hole so the mesher removes it.
    # xfemm's hole line is "x y InGroup"; Tangle reads only x and y.
    lines.append("[NumHoles] = 1\n")
    lines.append("0 0 0\n")
    lines.append("[NumBlockLabels] = 1\n")
    label = coord((INNER_RADIUS + OUTER_RADIUS) / 2.0, 0.0)
    lines.append(f"{number(label[0])} {number(label[1])} 1 0.001 0 0 0 1 0\n")

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("".join(lines))


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    write_tile(os.path.join(here, "d_ring_tile.fem"))
    write_control(os.path.join(here, "d_ring_control.fem"))


if __name__ == "__main__":
    main()
