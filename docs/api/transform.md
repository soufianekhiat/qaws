# API Reference: qaws_transform.h

A placement of an object built in its own frame: `world = translation + R (scale * x)`. The scale per axis of the
object's frame comes first, then the rotation, then the translation. The three parts are kept apart, never as a matrix,
so no shear can arise: a box of the object's frame scaled per axis is still a box once the rotation is undone. A
query on a placed object therefore maps into its frame without rebuilding anything built there (a tree, moments).

```c
typedef struct qaws_transform
{
	qaws_scalar rotation[4];      /* a quaternion (w, x, y, z), normalized where used; all zero: none */
	qaws_scalar scale[3];         /* per axis of the object's frame; a component 0 counts as 1; negative mirrors */
	qaws_scalar translation[3];   /* added last */
} qaws_transform;
```

- A zeroed `qaws_transform` is the identity.
- A negative scale component mirrors: the object's mirror image, its inside kept inside (a placed solid stays a solid,
  whatever the parity of its negative components).
- In 2D (x, y) the rotation must turn about z (its x and y components 0); the third components of `scale` and
  `translation` are ignored.

## qaws_transform_points, qaws_transform_points_inverse

```c
qaws_status qaws_transform_points(qaws_scalar* out_points, qaws_transform const* transform, qaws_scalar const* points,
	unsigned int count);
qaws_status qaws_transform_points_inverse(qaws_scalar* out_points, qaws_transform const* transform, qaws_scalar const* points,
	unsigned int count);
```

`count` points, 3 coordinates each, from the object's frame to the world (`_inverse`: back). `out_points` may be
`points`. `QAWS_STATUS_INVALID_ARGUMENT` for a NULL transform, NULL points or output with a count, or a value that is
not finite.
