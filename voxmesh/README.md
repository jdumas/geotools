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

Disclaimer
----------

Because `voxmesh` traces rays in the Z direction only, it might introduce some directional bias during the voxelization. This is especially important if you are trying to voxelize a mesh with cubic symmetry, because there is no guarantee that `voxmesh` will produce a result with cubic symmetric. A more sensible approach would be to combine raytracing voxelization from the 3 directions (X, Y and Z). This is left as future work.

Moreover, the algorithm employed here is rather simple, and in particular it does not try to minimize the total volume error between the voxelized surface and the input mesh. I am not aware of any open-source implementation that explicitly tries to minimize this volumetric error, but I would be happy to hear about it.
