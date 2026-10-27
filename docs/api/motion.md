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
- **Many turns:** a quaternion says where the body ends, not how many
  turns it took; the angular velocities say it. The total angle is the
  relative rotation's plus the whole turns that bring it nearest the
  angular velocities' mean along its axis (three turns about z: equal end
  rotations, ω = 6π). Past a half turn the motion is cut into equal pieces
  of an eighth of a turn at most, each a pose motion as above: the angle
  along the axis a Hermite polynomial, the velocities' other parts blended
  linearly, the translation's polynomial cut at the joints (exact), so
  positions, orientations and their rates match where pieces meet.
- **Vertex motions:** each vertex its own Hermite polynomial from its
  positions and velocities at both ends (accelerations with the flag).

## Limits

- A half turn the angular velocities do not orient (either way round as
  near): INVALID_ARGUMENT.
- Past a half turn, the off-axis parts of the angular velocities are
  blended linearly between the ends (a tumbling body given by its two end
  poses alone is approximated at the joints).
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
| Three turns about z moving along x (equal end rotations, ω = 6π): against the uniform spin at 601 times | 1.2e-3 |
| A half turn with ω = −π: turns that way (R at t = 1/2) | exact |
| A half turn its velocities do not orient, a zero quaternion | INVALID_ARGUMENT |
