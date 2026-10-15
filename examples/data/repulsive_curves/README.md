# Repulsive Curves scenes

Curves from the scenes of Repulsive Curves (Christopher Yu, Henrik
Schumacher, Keenan Crane, ACM Transactions on Graphics 2021), copied
unchanged from https://github.com/icethrush/repulsive-curves (commit
78a0c15), under the MIT license in LICENSE. Read by
`examples/repulsive_paper.c`.

| File | Scene | Settings there |
|---|---|---|
| freedmancurve600.obj | FreedmanUnknot (the scene uses the 1200-vertex sampling) | alpha 2, beta 4.5; edge lengths, barycenter |
| k5.obj, k33.obj | GraphDrawing | alpha 3, beta 6; edge lengths, barycenter |
| earbuds.obj | CableUntangle | alpha 3, beta 6; edge lengths, barycenter |
| graph-sampled.obj | PlanarGraphDrawing | alpha 3, beta 6; edge lengths, barycenter; in the plane y = 0 |

OBJ curves: `v x y z` vertices and `l i j` edges (1-based).
