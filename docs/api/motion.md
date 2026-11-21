# API Reference: qaws_motion.h

Motions over t in [0, 1], the input of sweeps: a body moving, and what it covers. A rigid motion moves points of the
body's frame, `x(t) = R(t) p + T(t)`; a vertex motion moves each vertex of a deforming mesh on its own trajectory.
Thread-safe on an immutable motion.

Names are shown unsuffixed; each exists as `_f32` and `_f64` (see [precision.md](precision.md)). Internally the
motion evaluates in double.

```c
qaws_motion_pose s = { position, { w, x, y, z }, velocity, angular_velocity };
qaws_motion_pose e = { ... };
qaws_motion* m;
qaws_motion_create_poses(&m, &s, &e, 0);                    /* cubic; QAWS_MOTION_QUINTIC: accelerations too */
qaws_motion_create_curve(&m, curve, rotation);              /* along a 3D curve, a constant rotation */
qaws_motion_create_curve_frames(&m, curve, rotation, QAWS_MOTION_FRAME_ROTATION_MINIMIZING);   /* turning along it */
qaws_motion_create_vertices(&m, p0, p1, v0, v1, a0, a1, count, flags);   /* per vertex */

qaws_motion_evaluate(&state, m, t);                         /* position, rotation (by rows), velocity, angular velocity */
qaws_motion_transform_points(out, out_velocities, m, t, points, count);
qaws_motion_destroy(m);
```

Related: [curve.md](curve.md), [view.md](view.md) (`qaws_array_view`), [inspect.md](inspect.md)
(rotation-minimizing frames of a curve).

---

## Types and flags

```c
typedef struct qaws_motion qaws_motion;   /* opaque */

#define QAWS_MOTION_QUINTIC 1u      /* the accelerations matched too (quintic Hermite) */

typedef struct qaws_motion_pose
{
    qaws_vec3 position;
    qaws_scalar rotation[4];        /* a quaternion (w, x, y, z), normalized when made */
    qaws_vec3 velocity;
    qaws_vec3 angular_velocity;     /* in the world frame */
    qaws_vec3 acceleration;         /* with QAWS_MOTION_QUINTIC */
    qaws_vec3 angular_acceleration; /* with QAWS_MOTION_QUINTIC, in the world frame */
} qaws_motion_pose;

/* A rigid motion at t: x = rotation p + position (rotation by rows), its velocity and angular velocity (world frame) */
typedef struct qaws_motion_state
{
    qaws_vec3 position;
    qaws_scalar rotation[9];
    qaws_vec3 velocity;
    qaws_vec3 angular_velocity;
} qaws_motion_state;

typedef enum qaws_motion_frame
{
    QAWS_MOTION_FRAME_FIXED = 0,             /* the rotation held (as qaws_motion_create_curve) */
    QAWS_MOTION_FRAME_ROTATION_MINIMIZING    /* the rotation carried along the curve, no turn about its tangent */
} qaws_motion_frame;
```

---

## Creating motions

```c
qaws_status qaws_motion_create_poses(qaws_motion** out_motion, qaws_motion_pose const* start, qaws_motion_pose const* end,
    unsigned int flags);

qaws_status qaws_motion_create_curve(qaws_motion** out_motion, qaws_curve const* curve, qaws_scalar const* rotation);

qaws_status qaws_motion_create_curve_frames(qaws_motion** out_motion, qaws_curve const* curve, qaws_scalar const* rotation,
    qaws_motion_frame frame);

qaws_status qaws_motion_create_vertices(qaws_motion** out_motion, qaws_scalar const* positions0, qaws_scalar const* positions1,
    qaws_scalar const* velocities0, qaws_scalar const* velocities1, qaws_scalar const* accelerations0,
    qaws_scalar const* accelerations1, unsigned int count, unsigned int flags);

qaws_status qaws_motion_create_vertices_borrowed(qaws_motion** out_motion, qaws_array_view const* positions0,
    qaws_array_view const* positions1);

void qaws_motion_destroy(qaws_motion* motion);
```

- `create_poses`: a rigid motion between two poses (see "How it works"). `flags`: 0 (cubic) or
  `QAWS_MOTION_QUINTIC`. `QAWS_STATUS_INVALID_ARGUMENT` for a value not finite, a zero quaternion, or a half turn
  the angular velocities do not orient.
- `create_curve`: a rigid motion along a 3D curve with a constant rotation; `rotation` is a quaternion (w, x, y, z),
  NULL for the identity. `QAWS_STATUS_INVALID_ARGUMENT` for a curve that is not 3D.
- `create_curve_frames`: as `create_curve`, the rotation being the body's at t = 0, turned along the curve by
  `frame` (`QAWS_MOTION_FRAME_FIXED` holds it, `QAWS_MOTION_FRAME_ROTATION_MINIMIZING` carries it).
  `QAWS_STATUS_INVALID_ARGUMENT` as above or for an unknown frame.
- `create_vertices`: `count` vertices (3 scalars each) with positions and velocities at t = 0 and 1, and with
  `QAWS_MOTION_QUINTIC` accelerations (either array may be NULL: zero). The arrays are copied.
  `QAWS_STATUS_INVALID_ARGUMENT` for `count == 0`, a null position or velocity array, or a value not finite.
- `create_vertices_borrowed`: a vertex motion moving each vertex on its chord at a constant rate from `positions0`
  (t = 0) to `positions1` (t = 1). Both views have 3 components and the same non-zero count, and are borrowed for
  the motion's life (read at each evaluation). `QAWS_STATUS_INVALID_ARGUMENT` for invalid views, counts that differ,
  a value not finite.
- The curve of a curve motion is borrowed: keep it until the motion is destroyed.
- `*out_motion` is NULL on failure; `QAWS_STATUS_ALLOCATION_FAILURE` when memory runs out.
  `qaws_motion_destroy` frees a motion; NULL is a no-op.

---

## Queries and evaluation

```c
int qaws_motion_is_per_vertex(qaws_motion const* motion);
unsigned int qaws_motion_get_vertex_count(qaws_motion const* motion);

qaws_status qaws_motion_evaluate(qaws_motion_state* out_state, qaws_motion const* motion, qaws_scalar t);

qaws_status qaws_motion_transform_points(qaws_scalar* out_points, qaws_scalar* out_velocities, qaws_motion const* motion,
    qaws_scalar t, qaws_scalar const* points, unsigned int count);
```

- `is_per_vertex`: 1 for a vertex motion, 0 for a rigid one (and for NULL). `get_vertex_count`: a vertex motion's
  count, 0 otherwise.
- `evaluate`: a rigid motion's state at `t`. `QAWS_STATUS_INVALID_ARGUMENT` for a vertex motion, a null pointer or
  a `t` not finite.
- `transform_points`: points (3 scalars each) moved to `t`. For a rigid motion, `count` points of the body's frame
  become `R(t) p + T(t)`, with velocities `omega x (R p) + T'`. For a vertex motion, pass `points` NULL and `count`
  equal to its vertex count: its own vertices are written. `out_velocities` may be NULL.
  `QAWS_STATUS_INVALID_ARGUMENT` for a null output, a `t` not finite, a vertex count that differs, or missing points.
- `t` outside [0, 1] extends the end pieces of a pose motion; a curve motion clamps the frame lookup to [0, 1].

---

## How it works

- **Translations** are Hermite polynomials: cubic through the two positions and velocities, quintic with the
  accelerations too (`QAWS_MOTION_QUINTIC`). A cubic path given by its ends' values and derivatives is reproduced
  exactly, and each trajectory is a Bezier curve in time.
- **Rotations** are `R(t) = R0 Cay(c(t))`: `c(t)` a Hermite polynomial of the same degree in the Gibbs vector (the
  axis times tan(angle / 2)) of the rotation from pose 0, and `Cay(c) = I + 2 ([c]x + [c]x^2) / (1 + c . c)` the
  Cayley map. The rotation is a rational function of t (entries polynomials over `1 + c . c`): exact, without
  trigonometry, orthonormal for any c. `c(0) = 0` and `c(1)` is the relative rotation's. The ends' derivatives
  come from the angular velocities (world frame) through the Gibbs vector's kinematics,
  `omega = R0 (2 / (1 + c . c)) (c' + c x c')`, and its derivative for the angular accelerations.
- **Many turns:** a quaternion says where the body ends, not how many turns it took; the angular velocities say
  it. The total angle is the relative rotation's (about its axis) plus the whole turns that bring it nearest the
  angular velocities' mean along that axis (a spin of three turns ends where it started, its velocities 6 pi). Past
  a half turn the motion is cut into equal pieces of an eighth of a turn at most, each as above: the angle along
  the axis a Hermite polynomial (its rates the angular velocities' along the axis), the other parts blended
  linearly, the translation's polynomial cut at the joints (exact), so positions, orientations and their rates
  match where the pieces meet.
- **Along a curve:** `T(t)` is the curve at its parameter mapped linearly from t (range start at 0, end at 1), its
  velocity the curve's derivative times the range; `R(t)` the given rotation. Rotation-minimizing:
  `R(t) = F(t) F(0)^T R0`, F the frame of the unit tangent and a normal transported by double reflection (Wang,
  Juttler, Zheng and Liu 2008) from 512 samples, each evaluation one reflection on from the sample before t, so the
  frame is continuous; the angular velocity is `c' x c'' / abs(c')^2` (none along the tangent). A planar curve
  turns the body about the plane's normal only.
- **Vertex motions:** each vertex its own Hermite polynomial from its positions and velocities at both ends
  (accelerations with the flag); borrowed vertex motions move on straight chords.

## Limits

- A half turn the angular velocities do not orient (either way round as near) is refused
  (`QAWS_STATUS_INVALID_ARGUMENT`).
- Past a half turn, the off-axis parts of the angular velocities are blended linearly between the ends (a tumbling
  body given by its two end poses alone is approximated at the joints).
- Affine motions (scaling, shear) are not offered.
