#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple test file to load voxel data.
"""
import re
import os
import numpy


def read_itk_metaimage(filename, order=[2, 1, 0]):
    """
    Read voxel data saved as a pair of .mhd/.raw files.
    Assume the binary data is stored as uint8.

    Additional documentation:
    http://www.itk.org/Wiki/ITK/MetaIO/Documentation

    Args:
        filename (str): Name of the file to load (with or withour extension)
        order (list, optional): Storage order in the .raw file.
    """
    name = os.path.splitext(filename)[0]
    with open(name + ".mhd", "r") as mhdfile:
        s = mhdfile.read()
    m = re.search('DimSize = ([0-9]*) ([0-9]*) ([0-9]*)', s)
    shape = (int(m.group(order[0] + 1)), int(m.group(order[1] + 1)), int(m.group(order[2] + 1)))
    m = re.search('ElementDataFile = (.*).raw', s)
    my_dir = os.path.dirname(name)
    assert os.path.join(my_dir, m.group(1)) == name
    with open(name + ".raw", "rb") as rawfile:
        volume = numpy.frombuffer(bytearray(rawfile.read()), dtype=numpy.uint8)
    return volume.reshape(shape)


def plot_voxels(volume):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from skimage import measure

    # Use marching cubes to obtain the surface mesh
    verts, faces = measure.marching_cubes(volume, 0.5)

    # Display resulting triangular mesh using Matplotlib. This can also be done
    # with mayavi (see skimage.measure.marching_cubes docstring).
    fig = plt.figure(figsize=(10, 12))
    ax = fig.add_subplot(111, projection='3d')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(verts[faces])
    ax.add_collection3d(mesh)

    ax.set_xlim(0, volume.shape[0])
    ax.set_ylim(0, volume.shape[1])
    ax.set_zlim(0, volume.shape[2])

    ax.set_aspect('equal')
    plt.show()


if __name__ == "__main__":
    filename = 'build/output.mhd'
    volume = read_itk_metaimage(filename)
    print("volume = {0}%".format(100 * volume.sum() / volume.size))
    print(volume.shape)
    plot_voxels(volume)
