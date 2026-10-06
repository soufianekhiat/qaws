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
 * Supported: Bezier and rational Bezier, 2D and 3D, degrees up to the
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

#endif /* QAWS_EXACT_H */
