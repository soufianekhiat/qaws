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
 * (overlapping curves) or a crossing through a singular point; a hit on a
 * span end at an irrational root is decided by an exact gcd test.
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

/*
 * Certified intersections of N exact curves at once (all 2D or all 3D, one
 * exact space). Every span of every curve gets a sound box from its control
 * points; one grid over all of them gives the span pairs that may meet, and
 * only those are certified as by qaws_exact_curve_curve_hits. `families`
 * (optional) skips the pairs of curves with the same id; with
 * QAWS_EXACT_BATCH_SELF each curve's self-intersections are added.
 *
 * Hits are sorted by (curve_a, curve_b, a_lo, b_lo), curve_a < curve_b
 * (curve_a == curve_b for a self-intersection). A curve pair that cannot be
 * certified (tangency, common component, ...) contributes no hits and is
 * counted in stats.uncertified_count; the call then completes every other
 * pair and returns QAWS_STATUS_CERTIFICATION_FAILED. *out_count is the
 * total; past `capacity` the call returns QAWS_STATUS_BUFFER_TOO_SMALL.
 */
#define QAWS_EXACT_BATCH_SELF 1u

typedef struct qaws_exact_batch_desc
{
	qaws_exact_curve const* const* curves;
	unsigned int curve_count;
	unsigned int const* families;   /* optional, one id per curve */
	unsigned int flags;             /* QAWS_EXACT_BATCH_* */
} qaws_exact_batch_desc;

typedef struct qaws_exact_batch_hit
{
	unsigned int curve_a;
	unsigned int curve_b;
	qaws_exact_pair pair;           /* a_* on curve_a, b_* on curve_b */
} qaws_exact_batch_hit;

typedef struct qaws_exact_batch_stats
{
	unsigned int span_count;        /* spans over all curves */
	unsigned int cell_count;        /* grid cells */
	unsigned int candidate_count;   /* span pairs whose boxes overlap */
	unsigned int uncertified_count; /* curve pairs that failed certification */
	unsigned int hit_count;
	unsigned int patch_count;       /* surface patches (curve / surface batch) */
} qaws_exact_batch_stats;

qaws_status qaws_exact_curve_batch_hits(qaws_exact_batch_desc const* desc, qaws_exact_batch_hit* out_hits, unsigned int capacity,
	unsigned int* out_count, qaws_exact_batch_stats* out_stats);

/*
 * Certified closest points: for every query point (exact doubles, world
 * units), the nearest point over the exact curves of desc (families and
 * flags ignored). On a span C = X / W the squared distance to p = P / Q is
 * stationary at the roots of the integer polynomial
 *   g(s) = sum_c (Q X_c - P_c W)(X_c' W - X_c W')
 * (degree 3n - 2); its roots, isolated exactly, and the span ends are the
 * candidates, each with an enclosure of its distance. Spans come from one
 * grid over their sound control boxes, nearest first, while they can still
 * win. The result is certified when the winner's distance enclosure lies
 * below every other candidate's (otherwise the nearest candidate is still
 * reported, certified = 0: an exact tie, e.g. the centre of a circle, or
 * candidates closer than the refinement can separate).
 */
#define QAWS_EXACT_CLOSEST_NONE 0xFFFFFFFFu

typedef struct qaws_exact_closest_point
{
	unsigned int curve;             /* QAWS_EXACT_CLOSEST_NONE without curves */
	double t_lo, t_hi;              /* the nearest point's parameter (lo == hi: exact) */
	double distance_lo, distance_hi;/* its distance, world units */
	int certified;                  /* the nearest point is proven unique */
} qaws_exact_closest_point;

/* points: point_count points of the curves' dimension. stats: span_count,
   candidate_count (spans solved), hit_count, uncertified_count. */
qaws_status qaws_exact_curve_batch_closest(qaws_exact_batch_desc const* desc, double const* points, unsigned int point_count,
	qaws_exact_closest_point* out_points, qaws_exact_batch_stats* out_stats);

/* ===================================================================
 * Certified 2D Boolean operations
 * =================================================================== */

/* A piece of an input boundary: curve a (region 0) or b (region 1) between two parameters. */
typedef struct qaws_exact_piece
{
	unsigned int region;    /* 0: a, 1: b */
	double t0_lo, t0_hi;    /* where the piece starts along its traversal (an enclosure; lo == hi when exact) */
	double t1_lo, t1_hi;    /* where it ends */
	int reversed;           /* traversed against the curve's parameter */
} qaws_exact_piece;

/* A result loop: pieces[first .. first + count - 1], each ending where the next starts. */
typedef struct qaws_exact_loop
{
	unsigned int first, count;
} qaws_exact_loop;

/*
 * Certified Boolean of two regions bounded by closed, simple 2D exact
 * curves (operation: QAWS_BOOLEAN_UNION, _INTERSECTION or _DIFFERENCE,
 * a minus b). The result boundary is made of pieces of the inputs
 * themselves (nothing is refitted): the boundaries are cut at their
 * certified crossings, each piece is classified by the certified winding
 * number of the other boundary around an exact rational point inside it
 * (non-zero: inside), and the kept pieces are linked into loops (two kept
 * pieces meet at a crossing; at a touching point the same curve goes on).
 * Loops are traced in the direction of their first piece.
 * QAWS_STATUS_CERTIFICATION_FAILED when a crossing cannot be certified
 * (a tangency, an overlap) or a boundary is not simple.
 */
qaws_status qaws_exact_boolean_2d(qaws_exact_curve const* a, qaws_exact_curve const* b, unsigned int operation, qaws_exact_piece* out_pieces,
	unsigned int piece_capacity, unsigned int* out_piece_count, qaws_exact_loop* out_loops, unsigned int loop_capacity, unsigned int* out_loop_count);

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

/* A line / surface intersection: patch parameters and the line parameter, enclosed. */
typedef struct qaws_exact_surface_hit
{
	qaws_exact_hit_kind kind;   /* POINT: u, v and t exact (lo == hi); CROSSING: enclosures */
	double u_lo, u_hi;
	double v_lo, v_hi;
	double t_lo, t_hi;          /* x = p0 + t (p1 - p0) */
} qaws_exact_surface_hit;

/*
 * Certified intersections of the line through the double points p0 != p1
 * with an exact surface (exact ray casting: keep t >= 0 for a ray). Per
 * patch, the line is the meet of two planes; their equations in (u, v)
 * are eliminated exactly (a Bezout matrix in u with polynomial entries in
 * v, product of the degrees bounded by the root budget), the roots in v
 * isolated, u recovered from the cofactors and proven in [0, 1]; t is
 * enclosed over the exact sub-patch. Every intersection is reported once,
 * sorted by t. QAWS_STATUS_CERTIFICATION_FAILED for a tangent line, a line
 * lying on the surface, or a hit on a patch edge at an irrational
 * parameter.
 */
qaws_status qaws_exact_surface_line_hits(qaws_exact_surface const* surface, double const p0[3], double const p1[3], qaws_exact_surface_hit* out_hits,
	unsigned int capacity, unsigned int* out_count);

/* A curve / surface intersection: the curve's parameter and the surface's, enclosed. */
typedef struct qaws_exact_curve_surface_hit
{
	double t_lo, t_hi;   /* curve parameter */
	double u_lo, u_hi;
	double v_lo, v_hi;
} qaws_exact_curve_surface_hit;

/*
 * Certified intersections of a 3D exact curve with an exact surface (one
 * exact space). Per span and patch, the three equations
 * X_c(u, v) W_C(t) - C_c(t) W(u, v) = 0 have exact integer Bernstein
 * coefficients on the (u, v, t) box; boxes are subdivided exactly and
 *   - excluded when one equation's coefficients all have one strict sign;
 *   - certified to hold exactly one intersection when, after an
 *     approximate-inverse preconditioning G = Y F (any Y is valid: the
 *     test is exact), G passes the Poincare-Miranda sign test on the box
 *     faces (existence) and its interval Jacobian, bounded by Bernstein
 *     coefficients, is strictly diagonally dominant (G injective: unique).
 * Certified boxes are then shrunk by bisection (about 2^-40). Every
 * intersection is reported once. QAWS_STATUS_CERTIFICATION_FAILED when the
 * subdivision budget runs out (a tangency, the curve lying on the
 * surface, or crossings closer than the budget separates).
 */
qaws_status qaws_exact_curve_surface_hits(qaws_exact_curve const* curve, qaws_exact_surface const* surface, qaws_exact_curve_surface_hit* out_hits,
	unsigned int capacity, unsigned int* out_count);

/*
 * Certified intersections of N 3D exact curves with M exact surfaces at
 * once (one exact space). Every curve span and surface patch gets a sound
 * box from its control points; one grid over all of them gives the
 * span / patch pairs that may meet, and only those are certified as by
 * qaws_exact_curve_surface_hits. Hits are sorted by (curve, surface, t_lo),
 * a root on a span or patch edge reported once. A curve / surface pair that
 * cannot be certified contributes no hits and is counted in
 * stats.uncertified_count (the call then returns
 * QAWS_STATUS_CERTIFICATION_FAILED); past `capacity`,
 * QAWS_STATUS_BUFFER_TOO_SMALL. stats.span_count / patch_count count the pieces.
 */
typedef struct qaws_exact_surface_batch_desc
{
	qaws_exact_curve const* const* curves;
	unsigned int curve_count;
	qaws_exact_surface const* const* surfaces;
	unsigned int surface_count;
} qaws_exact_surface_batch_desc;

typedef struct qaws_exact_curve_surface_batch_hit
{
	unsigned int curve;
	unsigned int surface;
	qaws_exact_curve_surface_hit hit;
} qaws_exact_curve_surface_batch_hit;

qaws_status qaws_exact_curve_surface_batch_hits(qaws_exact_surface_batch_desc const* desc, qaws_exact_curve_surface_batch_hit* out_hits,
	unsigned int capacity, unsigned int* out_count, qaws_exact_batch_stats* out_stats);

/* A point of a surface / surface intersection curve: (u1, v1) on a, (u2, v2) on b, enclosed. */
typedef struct qaws_exact_ssi_point
{
	double u1_lo, u1_hi, v1_lo, v1_hi;
	double u2_lo, u2_hi, v2_lo, v2_hi;
} qaws_exact_ssi_point;

/* A branch: points[first .. first + count - 1]; closed when the last joins the first. */
typedef struct qaws_exact_ssi_branch
{
	unsigned int first, count;
	int closed;
} qaws_exact_ssi_branch;

/*
 * Certified surface / surface intersection curves (one exact space). Per
 * patch pair, the three equations X^a_c W^b - X^b_c W^a = 0 in
 * (u1, v1, u2, v2) have exact integer Bernstein coefficients; the 4D box
 * is subdivided exactly until each box is excluded or regular: a 3 x 3
 * minor of its interval Jacobian keeps one strict sign, so every solution
 * arc in it is a monotone graph over the remaining variable and ends on
 * the box boundary. The arcs' end points on the box faces are certified
 * points (three equations in three unknowns); a regular box with none
 * holds no curve, one with two holds exactly one arc joining them.
 * min_depth (0..10) subdivides regular boxes at least that many times per
 * direction: more points along the curves, still certified.
 * Branches chain those points: between consecutive points lies exactly
 * one smooth intersection arc, and every intersection curve is covered.
 * QAWS_STATUS_CERTIFICATION_FAILED for tangential contact, overlapping
 * surfaces or a singular intersection point (the budget runs out).
 */
qaws_status qaws_exact_surface_surface_hits(qaws_exact_surface const* a, qaws_exact_surface const* b, unsigned int min_depth, qaws_exact_ssi_point* out_points,
	unsigned int point_capacity, unsigned int* out_point_count, qaws_exact_ssi_branch* out_branches, unsigned int branch_capacity,
	unsigned int* out_branch_count);

/*
 * Certified intersection curves of N exact surfaces at once (one exact
 * space). Every patch gets a sound box from its control net; one grid over
 * all patches gives the patch pairs of different surfaces (and families,
 * when given) whose boxes overlap, and each surface pair with any is
 * solved as by qaws_exact_surface_surface_hits on those patch pairs only.
 * Branches are grouped by (surface_a, surface_b), surface_a < surface_b;
 * their points index the shared point buffer. A surface pair that cannot be
 * certified contributes nothing and is counted in stats.uncertified_count
 * (the call then returns QAWS_STATUS_CERTIFICATION_FAILED); when the
 * buffers run out the call stops with QAWS_STATUS_BUFFER_TOO_SMALL.
 * stats.patch_count counts the patches, candidate_count the patch pairs.
 */
typedef struct qaws_exact_ssi_batch_desc
{
	qaws_exact_surface const* const* surfaces;
	unsigned int surface_count;
	unsigned int const* families;   /* optional, one id per surface */
	unsigned int min_depth;         /* as for qaws_exact_surface_surface_hits */
} qaws_exact_ssi_batch_desc;

typedef struct qaws_exact_ssi_batch_branch
{
	unsigned int surface_a;
	unsigned int surface_b;
	qaws_exact_ssi_branch branch;   /* first indexes the shared point buffer */
} qaws_exact_ssi_batch_branch;

qaws_status qaws_exact_surface_batch_hits(qaws_exact_ssi_batch_desc const* desc, qaws_exact_ssi_point* out_points, unsigned int point_capacity,
	unsigned int* out_point_count, qaws_exact_ssi_batch_branch* out_branches, unsigned int branch_capacity, unsigned int* out_branch_count,
	qaws_exact_batch_stats* out_stats);

/*
 * Certified closest points on surfaces: for every query point (exact
 * doubles, world units), the nearest point over the exact surfaces of desc
 * (families and min_depth ignored). On a patch S = X / W the squared
 * distance to p = P / Q is stationary where
 *   G_u = sum_c (Q X_c - P_c W)(X_c,u W - X_c W_u) = 0,
 *   G_v = sum_c (Q X_c - P_c W)(X_c,v W - X_c W_v) = 0;
 * with 3w - 1 = 0 the certified 3D solver isolates its interior roots,
 * and the patch edges, exact curves, add theirs and the corners (each edge
 * and corner once over the surface). Every candidate gets an enclosure of
 * its distance; the nearest is certified when its enclosure lies below
 * every other one (a tie, e.g. the centre of a sphere, stays uncertified).
 * Patches up to bicubic.
 */
typedef struct qaws_exact_surface_closest_point
{
	unsigned int surface;           /* QAWS_EXACT_CLOSEST_NONE without surfaces */
	double u_lo, u_hi, v_lo, v_hi;  /* the nearest point's parameters */
	double distance_lo, distance_hi;/* its distance, world units */
	int certified;
} qaws_exact_surface_closest_point;

/* stats: patch_count, candidate_count (patches solved), hit_count, uncertified_count */
qaws_status qaws_exact_surface_batch_closest(qaws_exact_ssi_batch_desc const* desc, double const* points, unsigned int point_count,
	qaws_exact_surface_closest_point* out_points, qaws_exact_batch_stats* out_stats);

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
