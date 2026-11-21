# API Reference: qaws_boolean_2d.h

The Boolean operation selector for 2D regions. The header declares only this
enum, so that code built on qaws names the operations with one shared type;
it declares no function. For the region clip types and results see
[clip.md](clip.md); for seams and shatter of regions see
[shatter.md](shatter.md).

```c
typedef enum qaws_boolean_op {
	QAWS_BOOLEAN_UNION = 0,
	QAWS_BOOLEAN_INTERSECTION,
	QAWS_BOOLEAN_DIFFERENCE
} qaws_boolean_op;
```

| Member | Value | Result |
|---|---|---|
| `QAWS_BOOLEAN_UNION` | 0 | in A or in B |
| `QAWS_BOOLEAN_INTERSECTION` | 1 | in A and in B |
| `QAWS_BOOLEAN_DIFFERENCE` | 2 | in A and not in B |

The header has no precision-bearing names.
