VoxMesh
=======

A simple voxelization program, which takes a surface mesh as input.


Usage
-----

Voxelize mesh:

    ./voxmesh ../bunny.stl numvox=128

Open result with Paraview:

    paraview output.mhd

Read binary data with Python:

    ./visualize.py

Description
-----------

The program launches rays in the Z direction, and computes intersections with the input mesh. The computation is pretty fast even with an output grid of 1024. Results can be visualized with Paraview. Patch welcome if you want to add another output format.
