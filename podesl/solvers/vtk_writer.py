
from __future__ import annotations
import numpy as np

def write_vtk_any(path, nodes, elems, point_data=None, cell_data=None):
    nodes = np.asarray(nodes, dtype=float)
    elems = np.asarray(elems, dtype=int)
    npts = nodes.shape[0]
    ncells = elems.shape[0]
    with open(path, "w", encoding="utf-8") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("PODESL Mesh\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {npts} float\n")
        if nodes.shape[1] == 2:
            for x,y in nodes: f.write(f"{x:.9e} {y:.9e} 0.0\n")
        else:
            for x,y,z in nodes: f.write(f"{x:.9e} {y:.9e} {z:.9e}\n")
        nper = elems.shape[1]
        f.write(f"CELLS {ncells} {ncells*(1+nper)}\n")
        if nper == 2:
            for e in elems: f.write(f"2 {e[0]} {e[1]}\n"); cell_type = 3
        elif nper == 3:
            for e in elems: f.write(f"3 {e[0]} {e[1]} {e[2]}\n"); cell_type = 5
        elif nper == 4:
            if nodes.shape[1] >= 3:
                for e in elems: f.write(f"4 {e[0]} {e[1]} {e[2]} {e[3]}\n"); cell_type = 10
            else:
                for e in elems: f.write(f"4 {e[0]} {e[1]} {e[2]} {e[3]}\n"); cell_type = 9
        elif nper == 8:
            for e in elems:
                f.write(f"8 {e[0]} {e[1]} {e[2]} {e[3]} {e[4]} {e[5]} {e[6]} {e[7]}\n")
            cell_type = 12
        else:
            raise ValueError("Unsupported element node count for VTK")
        f.write(f"CELL_TYPES {ncells}\n")
        for _ in range(ncells): f.write(f"{cell_type}\n")
        if point_data:
            f.write(f"POINT_DATA {npts}\n")
            for name, arr in point_data.items():
                a = np.asarray(arr)
                if a.ndim == 1:
                    f.write(f"SCALARS {name} float 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    for v in a: f.write(f"{float(v):.9e}\n")
                elif a.ndim == 2 and a.shape[1] in (2,3):
                    f.write(f"VECTORS {name} float\n")
                    if a.shape[1] == 2:
                        for vx,vy in a: f.write(f"{float(vx):.9e} {float(vy):.9e} 0.0\n")
                    else:
                        for vx,vy,vz in a: f.write(f"{float(vx):.9e} {float(vy):.9e} {float(vz):.9e}\n")
        if cell_data:
            f.write(f"CELL_DATA {ncells}\n")
            for name, arr in cell_data.items():
                a = np.asarray(arr)
                if a.ndim == 1:
                    f.write(f"SCALARS {name} float 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    for v in a: f.write(f"{float(v):.9e}\n")
                elif a.ndim == 2 and a.shape[1] in (2,3):
                    f.write(f"VECTORS {name} float\n")
                    if a.shape[1] == 2:
                        for vx,vy in a: f.write(f"{float(vx):.9e} {float(vy):.9e} 0.0\n")
                    else:
                        for vx,vy,vz in a: f.write(f"{float(vx):.9e} {float(vy):.9e} {float(vz):.9e}\n")


def write_vtk_quad2d(
    path,
    nodes,
    elems,
    point_data=None,
    cell_data=None,
):
    write_vtk_any(path, nodes, elems, point_data=point_data, cell_data=cell_data)
