# Repulsive Surfaces scenes

Meshes from the scenes of Repulsive Surfaces (Christopher Yu, Caleb
Brakensiek, Henrik Schumacher, Keenan Crane, ACM Transactions on Graphics
2021), copied unchanged from https://github.com/icethrush/repulsive-surfaces
(commit c51b464), under the MIT license in LICENSE (figure8tube.obj and
doubletorusring.obj without their vertex normals, the others without texture
coordinates, bunny_points.obj its vertices only). Read by
`examples/repulsive_paper.c`.

| File | Scene | Settings there |
|---|---|---|
| hopf.obj | HopfLink (a scene of the paper's code, not a figure) | alpha 6, beta 12 (p 6); the area and the barycenter held, 200 iterations |
| borromean.obj | BorromeanRings (a scene of the paper's code, not a figure) | alpha 6, beta 12 (p 6); the area and the barycenter held, 350 iterations |
| doubletorus10k.obj | EntangledDoubleTorus (Fig. 1: linked handcuffs) | alpha 6, beta 12 (p 6); the area and the barycenter held, 700 iterations |
| doubletorusring.obj | DoubleTorusAndRing (DoubleTorusAndRingMesh_14464F.obj; a scene of the paper's code) | alpha 6, beta 12 (p 6); the area and the barycenter held, 200 iterations |
| figure8tube.obj | BrokenFigureEight (BrokenFigureEightTubeMesh_11740F.obj; a scene of the paper's code) | alpha 6, beta 12 (p 6); the area and the barycenter held, 500 iterations |
| walnut_inside.obj | Walnut (inside.obj; Fig. 22, top) | alpha 6, beta 12; the barycenter held, the area 8x over 200 steps, volume_deviation 10, obstacle walnut_shell.obj 0.5, 220 iterations |
| walnut_shell.obj | Walnut (shell2.obj), its obstacle | |
| sphere20k.obj | Reconstruction (sphere20k.obj; Fig. 20, the scene's sphere80k-half at a quarter of the faces, scaled by 0.6 to its radius 1.32 in the example) | alpha 6, beta 12; the volume 0.11x over 200 steps, the barycenter free, point_cloud_obstacle bunny_points.obj 1, 220 iterations |
| bunny_points.obj | Reconstruction (bunny-scaled.obj, its vertices only), the point cloud | |

OBJ meshes: `v x y z` vertices and `f i j k` faces (1-based).
