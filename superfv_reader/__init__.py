from dataclasses import dataclass

import numpy as np


@dataclass
class Mesh:
    nx: int
    ny: int
    nz: int
    xlim0: float
    xlim1: float
    ylim0: float
    ylim1: float
    zlim0: float
    zlim1: float


def read_snapshot(path, dtype=np.float64):
    """Return (t, mesh, u), with u shaped (5, nx, ny, nz).

    mesh is a Mesh dataclass with dimensions and domain bounds. dtype must match
    SUPERFV_REAL_TYPE; the binary format uses native byte order.
    """
    dtype = np.dtype(dtype)
    if dtype.kind != "f" or not dtype.isnative:
        raise ValueError("dtype must be a native floating-point type")

    header_size = 3 * np.dtype(np.uint64).itemsize + 7 * dtype.itemsize
    with open(path, "rb") as file:
        dimensions = np.fromfile(file, dtype=np.uint64, count=3)
        metadata = np.fromfile(file, dtype=dtype, count=7)
        if dimensions.size != 3 or metadata.size != 7:
            raise ValueError("incomplete snapshot header")
        nx, ny, nz = map(int, dimensions)

        size = 5 * nx * ny * nz
        if file.seek(0, 2) != header_size + size * dtype.itemsize:
            raise ValueError("snapshot size does not match its dimensions and dtype")
        file.seek(header_size)
        u = np.fromfile(file, dtype=dtype, count=size).reshape((5, nx, ny, nz))

    mesh = Mesh(nx, ny, nz, *map(float, metadata[:6]))
    return metadata[6], mesh, u
