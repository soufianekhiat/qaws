# API Reference: qaws_motion.h

Motions over t in [0, 1], the input of sweeps: a body moving and what it
covers. A rigid motion moves points of the body's frame, x(t) = R(t) p +
T(t); a vertex motion moves each vertex of a deforming mesh on its own
trajectory.

```c
qaws_motion_pose s = { position, { w, x, y, z }, velocity, angular_velocity };
qaws_motion_pose e = { ... };
qaws_motion* m;
qaws_motion_create_poses(&m, &s, &e, 0);                    /* cubic; QAWS_MOTION_QUINTIC: accelerations too */
qaws_motion_create_curve(&m, curve, rotation);              /* along a 3D curve, a constant rotation */
qaws_motion_create_vertices(&m, p0, p1, v0, v1, a0, a1, count, flags);   /* per vertex */

qaws_motion_evaluate(&state, m, t);                         /* position, rotation (by rows), velocity, angular velocity */
qaws_motion_transform_points(out, out_velocities, m, t, points, count);
qaws_motion_destroy(m);
```

## How it works

- **Translations** are Hermite polynomials: cubic through the two
  positions and velocities, quintic with the accelerations too
  (`QAWS_MOTION_QUINTIC`). A cubic path given by its ends' values and
  derivatives is reproduced exactly, and each trajectory is a Bezier curve
  in time (the swept triangle a polynomial patch, barycentric in space).
- **Rotations** are R(t) = R0 Cay(c(t)): c(t) a Hermite polynomial of the
  same degree in the Gibbs vector (the axis times tan(angle / 2)) of the
  rotation from pose 0, and Cay(c) = I + 2 ([c]× + [c]×²) / (1 + c·c) the
  Cayley map. The rotation is a rational function of t, its entries
  polynomials over 1 + c·c: exact, without trigonometry, orthonormal for
  any c. c(0) = 0 and c(1) is the relative rotation's (the shorter way
  round). The ends' derivatives come from the angular velocities (world
  frame) through the Gibbs vector's kinematics, ω = R0 (2 / (1 + c·c))
  (c′ + c × c′), inverted as c′ = (ω′ − c × ω′ + (c·ω′) c) / 2 with ω′ in
  pose 0's frame, and its derivative for the angular accelerations.
- **Along a curve:** T(t) is the curve at its parameter mapped linearly
  from t, its velocity the curve's derivative times the range; R(t) is the
  given rotation.
- **Vertex motions:** each vertex its own Hermite polynomial from its
  positions and velocities at both ends (accelerations with the flag).

## Limits

- A relative rotation of a half turn or more has no Gibbs vector:
  INVALID_ARGUMENT (cut the motion in two).
- Affine motions (scaling, shear) are not offered.
- A curve motion borrows its curve: it must outlive the motion.

## Tests

Test 179 (`tests/179_motion.c`, f64):

| Check | Result |
|---|---|
| Two poses, cubic and quintic: both ends' position, rotation, velocity and angular velocity | 3e-16 |
| R(t) orthonormal at 9 times | 4e-16 |
| The angular velocity against R's central differences; the velocity against T's | 5e-10; 8e-10 |
| Quintic: the accelerations and angular accelerations at both ends (differences of h = 1e-4) | 6e-6 |
| A cubic path from its ends' values and derivatives, at 17 times | exact |
| Three vertices' cubics; a Bezier curve's points and derivatives | exact |
| A half turn, a zero quaternion | INVALID_ARGUMENT |
