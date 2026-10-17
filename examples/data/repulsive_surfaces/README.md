# Repulsive Surfaces scenes

Meshes from the scenes of Repulsive Surfaces (Christopher Yu, Caleb
Brakensiek, Henrik Schumacher, Keenan Crane, ACM Transactions on Graphics
2021), copied unchanged from https://github.com/icethrush/repulsive-surfaces
(commit c51b464), under the MIT license in LICENSE (figure8tube.obj and
doubletorusring.obj without their vertex normals). Read by
`examples/repulsive_paper.c`.

| File | Scene | Settings there |
|---|---|---|
| hopf.obj | HopfLink (a scene of the paper's code, not a figure) | alpha 6, beta 12 (p 6); the area and the barycenter held, 200 iterations |
| borromean.obj | BorromeanRings (a scene of the paper's code, not a figure) | alpha 6, beta 12 (p 6); the area and the barycenter held, 350 iterations |
| doubletorus10k.obj | EntangledDoubleTorus (Fig. 1: linked handcuffs) | alpha 6, beta 12 (p 6); the area and the barycenter held, 700 iterations |
| doubletorusring.obj | DoubleTorusAndRing (DoubleTorusAndRingMesh_14464F.obj; a scene of the paper's code) | alpha 6, beta 12 (p 6); the area and the barycenter held, 200 iterations |
| figure8tube.obj | BrokenFigureEight (BrokenFigureEightTubeMesh_11740F.obj; a scene of the paper's code) | alpha 6, beta 12 (p 6); the area and the barycenter held, 500 iterations |

OBJ meshes: `v x y z` vertices and `f i j k` faces (1-based).
