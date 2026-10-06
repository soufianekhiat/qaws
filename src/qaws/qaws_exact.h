#ifndef QAWS_EXACT_H
#define QAWS_EXACT_H

/*
 * Exact / certified geometry (C only, built with QAWS_ENABLE_EXACT).
 *
 * Predicates on raw doubles: a finite double is an exact dyadic rational
 * m 2^e, so the sign of a polynomial in double inputs is a well defined
 * exact quantity. Each predicate evaluates in f64 first with a proven error
 * bound; when the bound does not prove the sign, it decodes the inputs
 * exactly and evaluates in multi-limb integers. The f64 filter never
 * decides a zero. Results are certified: never a wrong sign.
 *
 * QAWS_STATUS_EXACT_RANGE_EXCEEDED: the inputs span more binary exponents
 * than the integer budget covers (about 2^900 between the smallest and the
 * largest magnitude); QAWS_STATUS_INVALID_ARGUMENT for non-finite inputs.
 */

#include "qaws_status.h"
#include "qaws_types.h"
#include "qaws_surface_types.h"

typedef enum qaws_exact_sign
{
	QAWS_EXACT_NEGATIVE = -1,
	QAWS_EXACT_ZERO = 0,
	QAWS_EXACT_POSITIVE = 1
} qaws_exact_sign;

/* How a predicate decided (optional output). */
typedef enum qaws_exact_path
{
	QAWS_EXACT_PATH_FILTER = 0,   /* the f64 evaluation and its error bound */
	QAWS_EXACT_PATH_EXACT         /* exact integer evaluation */
} qaws_exact_path;

/* Sign of (a - c) x (b - c): positive when a, b, c turn counterclockwise. */
qaws_status qaws_exact_orient2d(double const a[2], double const b[2], double const c[2], qaws_exact_sign* out_sign,
	qaws_exact_path* out_path);

/* Sign of the volume (a - d) . ((b - d) x (c - d)). */
qaws_status qaws_exact_orient3d(double const a[3], double const b[3], double const c[3], double const d[3],
	qaws_exact_sign* out_sign, qaws_exact_path* out_path);

/* Sign of a / b - c / d (b, d non-zero). */
qaws_status qaws_exact_compare_ratio(double a, double b, double c, double d, qaws_exact_sign* out_sign, qaws_exact_path* out_path);

/* ===================================================================
 * Exact curves
 *
 * A curve is prepared once into integer homogeneous Bezier data:
 * coordinates are quantized onto the lattice x = i 2^space_exp2 (round to
 * nearest, ties to even; |i| below 2^coord_bits), weights are scaled by a
 * common power of two to weight_bits (a common factor cancels in the
 * rational curve; integer weights stay exact), and parameters onto
 * t = T 2^-param_bits. After that boundary nothing rounds: evaluation is
 * division-free De Casteljau on multi-limb integers, derivatives come from
 * integer derivative polygons and the exact quotient rule, and every value
 * is converted to double once, correctly rounded.
 *
 * Supported: Bezier, rational Bezier, B-spline, NURBS, cubic Hermite,
 * uniform Catmull-Rom and polynomial curves, 2D and 3D, degrees up to the
 * integer budget (QAWS_STATUS_EXACT_RANGE_EXCEEDED past it). Other families
 * return QAWS_STATUS_EXACT_UNSUPPORTED for now.
 * =================================================================== */

typedef struct qaws_exact_curve qaws_exact_curve;

typedef enum qaws_numeric_quality
{
	QAWS_NUMERIC_APPROXIMATE = 0,
	QAWS_NUMERIC_CERTIFIED,
	QAWS_NUMERIC_EXACT_RATIONAL,
	QAWS_NUMERIC_UNSUPPORTED
} qaws_numeric_quality;

typedef enum qaws_exact_flag
{
	QAWS_EXACT_FLAG_NONE = 0,
	QAWS_EXACT_FLAG_INPUT_QUANTIZED = 1 << 0,   /* inputs were rounded onto the lattices */
	QAWS_EXACT_FLAG_PREP_QUANTIZED = 1 << 1     /* an approximate preparation was frozen */
} qaws_exact_flag;

typedef struct qaws_exact_desc
{
	int space_exp2;              /* lattice step 2^space_exp2 world units */
	unsigned int coord_bits;     /* usable coordinate magnitude bits (default 26) */
	unsigned int param_bits;     /* parameter lattice bits (default 24) */
	unsigned int weight_bits;    /* weight bits (default 24) */
} qaws_exact_desc;

typedef struct qaws_exact_report
{
	qaws_numeric_quality quality;
	unsigned int flags;
	unsigned int storage_bits;               /* widest homogeneous control value */
	double max_position_quantization_error;  /* world units */
	double max_weight_quantization_error;    /* relative */
	double parameter_quantization_error;     /* last evaluation */
} qaws_exact_report;

/* Defaults: 2^-20 world units, 26 coordinate bits, 24 parameter and weight bits. */
void qaws_exact_desc_default(qaws_exact_desc* desc);

qaws_status qaws_exact_curve_prepare(qaws_exact_desc const* desc, qaws_curve const* curve, qaws_exact_curve** out_curve,
	qaws_exact_report* out_report);
void qaws_exact_curve_destroy(qaws_exact_curve* curve);

/*
 * Position and derivatives at t (quantized onto the parameter lattice),
 * evaluated exactly and rounded once: out receives (order + 1) points of
 * `dimension` doubles, C, C', ..., C^(order), order <= 3. out_report (may be
 * NULL) receives the parameter quantization error.
 */
qaws_status qaws_exact_curve_evaluate(qaws_exact_curve const* curve, double t, unsigned int order, double* out,
	qaws_exact_report* out_report);

/* Number of spans (one exact rational Bezier each). */
unsigned int qaws_exact_curve_span_count(qaws_exact_curve const* curve);

/*
 * Span s as a rational Bezier, rounded once: out_degree, its parameter
 * interval [t0, t1], `dimension` coordinates per control point (world
 * units) and weights relative to the largest one (any of the outputs may
 * be NULL; points and weights hold degree + 1 entries).
 */
qaws_status qaws_exact_curve_span_bezier(qaws_exact_curve const* curve, unsigned int s, unsigned int* out_degree, double* out_t0,
	double* out_t1, double* out_points, double* out_weights);

/*
 * Certified winding number of a closed 2D loop of exact curves (each
 * piece's end equal to the next one's start, exactly) around the double
 * point p: the signed crossings of the ray x > p_x, decided on exact
 * integer Bernstein coefficients by subdivision. The loop must share one
 * exact space (QAWS_STATUS_EXACT_INCOMPATIBLE_SPACE otherwise).
 * QAWS_STATUS_CERTIFICATION_FAILED when p lies on the curve (or within
 * 2^-max_depth of it in parameter, past the subdivision budget): the
 * winding number is then undefined, never guessed.
 */
qaws_status qaws_exact_winding_2d(qaws_exact_curve const* const* pieces, unsigned int count, double const p[2], int* out_winding);

/* ===================================================================
 * Certified intersections
 *
 * An intersection parameter is in general algebraic, not rational, so it
 * is returned as a certified enclosure: each crossing interval holds
 * exactly one intersection, a transversal one, and no intersection lies
 * outside the reported hits. The roots are isolated on exact integer
 * Bernstein coefficients (subdivision and Descartes' rule of signs).
 * =================================================================== */

typedef enum qaws_exact_hit_kind
{
	QAWS_EXACT_HIT_CROSSING = 0,   /* (t_lo, t_hi) holds exactly one transversal intersection */
	QAWS_EXACT_HIT_POINT,          /* an intersection at t_lo == t_hi exactly (enclosed when not a double) */
	QAWS_EXACT_HIT_OVERLAP         /* the curve lies in the line / plane on [t_lo, t_hi] */
} qaws_exact_hit_kind;

typedef struct qaws_exact_hit
{
	qaws_exact_hit_kind kind;
	double t_lo, t_hi;   /* curve parameters; the true value is inside, rounded outward */
} qaws_exact_hit;

/*
 * Intersections of a 2D exact curve with the infinite line through the
 * double points p0 != p1, crossing intervals shrunk to width <= width
 * (0: the narrowest the parameter lattice gives, about 2^-60 of a span).
 * QAWS_STATUS_CERTIFICATION_FAILED at a tangency (a multiple root cannot be
 * told from two close crossings by signs), QAWS_STATUS_BUFFER_TOO_SMALL
 * past capacity (out_count holds the number found so far).
 */
qaws_status qaws_exact_curve_line_hits(qaws_exact_curve const* curve, double const p0[2], double const p1[2], double width,
	qaws_exact_hit* out_hits, unsigned int capacity, unsigned int* out_count);

/* The same for a 3D exact curve and the plane through `point` with normal `normal` (doubles, taken exactly). */
qaws_status qaws_exact_curve_plane_hits(qaws_exact_curve const* curve, double const point[3], double const normal[3], double width,
	qaws_exact_hit* out_hits, unsigned int capacity, unsigned int* out_count);

/* An intersection of two curves: curve a's parameter in [a_lo, a_hi], curve b's in [b_lo, b_hi]. */
typedef struct qaws_exact_pair
{
	qaws_exact_hit_kind kind;   /* POINT: both parameters exact (lo == hi); CROSSING: enclosures */
	double a_lo, a_hi;
	double b_lo, b_hi;
} qaws_exact_pair;

/*
 * Certified intersections of two exact curves (both 2D or both 3D)
 * sharing one exact space. Span against span, the lower-degree span
 * (degree <= 6, product of the degrees <= 64) is implicitized exactly (its
 * Bezout matrix, in a coordinate plane for 3D), the other is substituted:
 * the roots of the resulting integer polynomial are isolated and each is
 * kept only if the implicitized span's own parameter, recovered exactly
 * from the matrix cofactors, is proven inside [0, 1]. In 3D the remaining
 * coordinate must agree: an exact zero test (a modular gcd proved by
 * division). Every intersection is reported once (a knot point by the
 * span ending there), each crossing enclosure holding exactly one.
 * QAWS_STATUS_CERTIFICATION_FAILED for tangencies, a common component
 * (overlapping curves), a crossing through a singular point, or an
 * intersection on a span end at an irrational parameter.
 */
qaws_status qaws_exact_curve_curve_hits(qaws_exact_curve const* a, qaws_exact_curve const* b, qaws_exact_pair* out_pairs, unsigned int capacity,
	unsigned int* out_count);

/*
 * Certified self-intersections of a 2D or 3D exact curve: pairs a < b with
 * C(a) = C(b) (a_lo..a_hi, b_lo..b_hi enclosures). Inside a span, the
 * divided differences (X(s) W(t) - X(t) W(s)) / (s - t) are eliminated the
 * same way (spans proven monotone along a direction are skipped); spans
 * against spans as for two curves. Shared knots and the closing point of a
 * closed curve are not self-intersections; spans on one algebraic curve
 * (a NURBS circle) are told apart exactly (touching ends) or refused
 * (overlap). QAWS_STATUS_CERTIFICATION_FAILED as above, and for an
 * irrational cusp or a span folding back on itself.
 */
qaws_status qaws_exact_curve_self_hits(qaws_exact_curve const* curve, qaws_exact_pair* out_pairs, unsigned int capacity, unsigned int* out_count);

/* ===================================================================
 * Exact tensor-product surfaces
 *
 * Bezier, B-spline and NURBS surfaces are prepared into one integer
 * homogeneous Bezier patch per non-empty knot rectangle, with the same
 * lattices as the curves (coordinates, weights, knots). Partial
 * derivatives up to second order come from integer derivative nets and
 * the exact bivariate quotient rule; the normal Su x Sv is exact too
 * (unnormalized: its unit length needs a square root).
 * =================================================================== */

typedef struct qaws_exact_surface qaws_exact_surface;

qaws_status qaws_exact_surface_prepare(qaws_exact_desc const* desc, qaws_surface const* surface, qaws_exact_surface** out_surface,
	qaws_exact_report* out_report);
void qaws_exact_surface_destroy(qaws_exact_surface* surface);

/*
 * Values at (u, v) (quantized onto the parameter lattices), exact and
 * rounded once, as vec3 entries: order 0 gives S; order 1 adds Su, Sv;
 * order 2 adds Suu, Suv, Svv (out holds 3, 9 or 18 doubles).
 * out_normal (may be NULL) receives Su x Sv, rounded once.
 */
qaws_status qaws_exact_surface_evaluate(qaws_exact_surface const* surface, double u, double v, unsigned int order, double* out,
	double* out_normal);

/* Patches along u and v (one exact rational Bezier patch each). */
void qaws_exact_surface_patch_count(qaws_exact_surface const* surface, unsigned int* out_u_count, unsigned int* out_v_count);

/*
 * Patch (iu, iv) rounded once: degrees, parameter rectangle, control
 * points (row-major, (p + 1) (q + 1) vec3, u index major) and weights
 * relative to the largest (any output may be NULL).
 */
qaws_status qaws_exact_surface_patch_bezier(qaws_exact_surface const* surface, unsigned int iu, unsigned int iv, unsigned int* out_p,
	unsigned int* out_q, double out_rect[4], double* out_points, double* out_weights);

#endif /* QAWS_EXACT_H */
