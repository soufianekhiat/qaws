#ifndef QAWS_DIFF_TYPES_H
#define QAWS_DIFF_TYPES_H

#include "qaws_types.h"
#include "qaws_status.h"

/*
 * Differentiation model.
 *
 * Every differentiable qaws object exposes typed parameter fields
 * (control points, weights, knots, radius, ...). Tangents and adjoints
 * are stored in views that mirror those fields, so callers never have to
 * pack a flat parameter vector. Flattening, bounds and optimizer variables
 * belong to the caller.
 *
 *   Tangent : forward propagation  y' = J x'
 *   Adjoint : reverse propagation  x_bar += J^T y_bar
 *   Tangent2: second directional derivative along the same tangent
 *
 * Derivatives with respect to the evaluation coordinate (t, or u/v) and
 * derivatives with respect to parameters are two independent axes. Results
 * are returned as spatial jets (see qaws_curve_jet_3d / qaws_surface_jet):
 * the tangent jet holds d/de of every spatial derivative.
 *
 * No numerical differentiation is ever returned through this contract.
 */

/* ===================================================================
 * Parameter fields
 *
 * Numeric values are stable: they are part of qaws_param_key.
 * Append only.
 * =================================================================== */

typedef enum qaws_diff_field
{
	QAWS_FIELD_NONE = 0,
	QAWS_FIELD_CONTROL_POINTS = 1,
	QAWS_FIELD_WEIGHTS = 2,
	QAWS_FIELD_KNOTS = 3,
	QAWS_FIELD_U_KNOTS = 4,
	QAWS_FIELD_V_KNOTS = 5,
	QAWS_FIELD_POINTS = 6,           /* interpolated points (Hermite, Catmull-Rom, ...) */
	QAWS_FIELD_DERIVATIVES = 7,      /* Hermite tangents */
	QAWS_FIELD_COEFFICIENTS = 8,     /* monomial coefficients */
	QAWS_FIELD_KEY_TIMES = 9,
	QAWS_FIELD_CENTER = 10,
	QAWS_FIELD_RADIUS = 11,
	QAWS_FIELD_RADIUS_B = 12,
	QAWS_FIELD_ANGLE_START = 13,
	QAWS_FIELD_ANGLE_END = 14,
	QAWS_FIELD_AXIS_U = 15,
	QAWS_FIELD_AXIS_V = 16,
	QAWS_FIELD_OFFSET_DISTANCE = 17,
	QAWS_FIELD_SCALE = 18,
	QAWS_FIELD_CURVATURE = 19,
	QAWS_FIELD_CURVATURE_RATE = 20,
	QAWS_FIELD_DIRECTION = 21,
	QAWS_FIELD_COUNT
} qaws_diff_field;

typedef enum qaws_value_type
{
	QAWS_VALUE_SCALAR = 1,
	QAWS_VALUE_VEC2 = 2,
	QAWS_VALUE_VEC3 = 3
} qaws_value_type;

typedef enum qaws_value_domain
{
	QAWS_DOMAIN_GENERIC = 0,
	QAWS_DOMAIN_POSITION,
	QAWS_DOMAIN_DIRECTION,      /* tangent vector / free vector */
	QAWS_DOMAIN_LENGTH,
	QAWS_DOMAIN_ANGLE,
	QAWS_DOMAIN_WEIGHT,
	QAWS_DOMAIN_TIME,
	QAWS_DOMAIN_PARAMETRIC,     /* curve/surface parameter coordinate */
	QAWS_DOMAIN_ARC_LENGTH,
	QAWS_DOMAIN_CURVATURE
} qaws_value_domain;

typedef enum qaws_value_constraint
{
	QAWS_CONSTRAINT_NONE = 0,
	QAWS_CONSTRAINT_POSITIVE,       /* > 0 */
	QAWS_CONSTRAINT_NON_NEGATIVE,   /* >= 0 */
	QAWS_CONSTRAINT_MONOTONIC,      /* elements non-decreasing (knot vectors) */
	QAWS_CONSTRAINT_UNIT_INTERVAL   /* in [0, 1] */
} qaws_value_constraint;

/* ===================================================================
 * Differentiability classes, capabilities, runtime validity
 * =================================================================== */

typedef enum qaws_diff_class
{
	QAWS_DIFF_SMOOTH = 0,
	QAWS_DIFF_PIECEWISE_SMOOTH,  /* smooth inside pieces (knot spans, segments) */
	QAWS_DIFF_ACTIVE_SET,        /* smooth while a discrete choice is unchanged */
	QAWS_DIFF_SUBGRADIENT,       /* non-smooth (max/min); a subgradient is returned */
	QAWS_DIFF_DISCRETE,          /* only derivatives under frozen topology exist */
	QAWS_DIFF_UNSUPPORTED
} qaws_diff_class;

typedef enum qaws_diff_capability
{
	QAWS_CAP_TANGENT    = 1 << 0,
	QAWS_CAP_ADJOINT    = 1 << 1,
	QAWS_CAP_TANGENT2   = 1 << 2,
	QAWS_CAP_DIRECT_HVP = 1 << 3,
	QAWS_CAP_LOCAL_SUPPORT = 1 << 4,  /* qaws_*_local_support is available */
	QAWS_CAP_LINEAR     = 1 << 5      /* output is linear in all fields (weights in support are exact) */
} qaws_diff_capability;

/* Discrete states a derivative may rely on. */
typedef enum qaws_diff_freeze
{
	QAWS_FREEZE_NONE       = 0,
	QAWS_FREEZE_SPAN       = 1 << 0,   /* knot span / segment selection */
	QAWS_FREEZE_ACTIVE_SET = 1 << 1,   /* closest primitive, branch, extremum */
	QAWS_FREEZE_TOPOLOGY   = 1 << 2,   /* trim loops, boolean topology, mesh connectivity */
	QAWS_FREEZE_SAMPLE_SET = 1 << 3    /* sample counts / adaptive sampling decisions */
} qaws_diff_freeze;

/* Runtime validity, ordered from best to worst. */
typedef enum qaws_diff_validity
{
	QAWS_DIFF_VALID = 0,
	QAWS_DIFF_VALID_LOCALLY,     /* valid while the frozen discrete state is unchanged */
	QAWS_DIFF_AT_BOUNDARY,       /* evaluated exactly on a piece boundary: one-sided */
	QAWS_DIFF_ILL_CONDITIONED,   /* implicit solve near singular */
	QAWS_DIFF_AMBIGUOUS,         /* several competing solutions */
	QAWS_DIFF_INVALID
} qaws_diff_validity;

typedef enum qaws_diff_accumulation
{
	QAWS_ACCUMULATE_SCATTER = 0, /* each sample adds into the parameter adjoints */
	QAWS_ACCUMULATE_GATHER,      /* each parameter gathers its samples (deterministic order) */
	QAWS_ACCUMULATE_TILED        /* tile-local accumulation, one flush per tile */
} qaws_diff_accumulation;

/* Coordinate kind of an evaluation parameter. */
typedef enum qaws_coordinate_kind
{
	QAWS_COORDINATE_PARAMETRIC = 0,
	QAWS_COORDINATE_NORMALIZED,
	QAWS_COORDINATE_ARC_LENGTH,
	QAWS_COORDINATE_TIME,
	QAWS_COORDINATE_WORLD_DISTANCE
} qaws_coordinate_kind;

/* ===================================================================
 * Field description (schema)
 * =================================================================== */

typedef struct qaws_field_desc
{
	qaws_diff_field field;
	qaws_value_type value_type;
	unsigned int count;               /* number of elements */
	qaws_value_domain domain;
	qaws_value_constraint constraint;
	qaws_scalar default_lower;        /* -inf/+inf encoded as -/+QAWS_DIFF_UNBOUNDED */
	qaws_scalar default_upper;
	qaws_diff_class diff_class;
	unsigned int capabilities;        /* qaws_diff_capability bits for this field */
} qaws_field_desc;

#define QAWS_DIFF_UNBOUNDED ((qaws_scalar)1e30)

/* ===================================================================
 * Parameter keys: stable structural identity of one scalar parameter.
 *
 *   child path / field / element / component
 *   e.g. "child[0]/control_points/4/z"
 * =================================================================== */

#define QAWS_PARAM_KEY_MAX_DEPTH 4

typedef struct qaws_param_key
{
	unsigned short child[QAWS_PARAM_KEY_MAX_DEPTH];
	unsigned char depth;      /* number of valid child indices */
	unsigned char field;      /* qaws_diff_field */
	unsigned char component;  /* 0..2 */
	unsigned char components; /* components of the field; 1 = scalar (no component in text) */
	unsigned int element;
} qaws_param_key;

/* ===================================================================
 * Differential storage views
 *
 * A view addresses one field of tangent or adjoint storage:
 *   element i, component c  ->  data[i * stride + c]
 * Missing fields are treated as zero tangents / inactive adjoints.
 * =================================================================== */

typedef struct qaws_field_view
{
	qaws_diff_field field;
	qaws_scalar* data;
	unsigned int count;               /* elements */
	unsigned int components;          /* 1, 2 or 3 */
	unsigned int stride;              /* scalars between elements; 0 = components */
	unsigned char const* active;      /* optional per-element activity, NULL = all active */
	unsigned int component_mask;      /* bit c set = component active; 0 = all */
} qaws_field_view;

typedef struct qaws_diff_views
{
	qaws_field_view* fields;
	unsigned int field_count;
	struct qaws_diff_views* children;  /* for derived objects, indexed like child keys */
	unsigned int child_count;
} qaws_diff_views;

/* ===================================================================
 * Spatial jets
 *
 * Curve jet: d[k] = d^k C / dt^k, k = 0..3.
 * channels uses QAWS_EVAL_FLAG_POSITION/D1/D2/D3 (bit k = order k).
 * =================================================================== */

#define QAWS_CURVE_JET_ORDER 3

typedef struct qaws_curve_jet_2d
{
	qaws_vec2 d[QAWS_CURVE_JET_ORDER + 1];
	unsigned int channels;
} qaws_curve_jet_2d;

typedef struct qaws_curve_jet_3d
{
	qaws_vec3 d[QAWS_CURVE_JET_ORDER + 1];
	unsigned int channels;
} qaws_curve_jet_3d;

/* Surface jet channels: all partial derivatives up to total order 3. */
typedef enum qaws_surface_jet_channel
{
	QAWS_SJET_P   = 1 << 0,
	QAWS_SJET_U   = 1 << 1,
	QAWS_SJET_V   = 1 << 2,
	QAWS_SJET_UU  = 1 << 3,
	QAWS_SJET_UV  = 1 << 4,
	QAWS_SJET_VV  = 1 << 5,
	QAWS_SJET_UUU = 1 << 6,
	QAWS_SJET_UUV = 1 << 7,
	QAWS_SJET_UVV = 1 << 8,
	QAWS_SJET_VVV = 1 << 9,
	QAWS_SJET_ORDER1 = QAWS_SJET_P | QAWS_SJET_U | QAWS_SJET_V,
	QAWS_SJET_ORDER2 = QAWS_SJET_ORDER1 | QAWS_SJET_UU | QAWS_SJET_UV | QAWS_SJET_VV,
	QAWS_SJET_ORDER3 = QAWS_SJET_ORDER2 | QAWS_SJET_UUU | QAWS_SJET_UUV | QAWS_SJET_UVV | QAWS_SJET_VVV
} qaws_surface_jet_channel;

#define QAWS_SURFACE_JET_COUNT 10

/* d[i] is the channel with bit i, i.e. d[0]=S, d[1]=Su, d[2]=Sv, d[3]=Suu,
   d[4]=Suv, d[5]=Svv, d[6]=Suuu, d[7]=Suuv, d[8]=Suvv, d[9]=Svvv. */
typedef struct qaws_surface_jet
{
	qaws_vec3 d[QAWS_SURFACE_JET_COUNT];
	unsigned int channels;
} qaws_surface_jet;

/* ===================================================================
 * Local support: which parameters influence one evaluation.
 *
 * For linear families (QAWS_CAP_LINEAR) the weights are exact:
 *   d^k C / dt^k = sum_j weights[r][k][j] * field_r[first_r + j]
 * Weights are given for derivative orders 0..order.
 * =================================================================== */

typedef enum qaws_support_kind
{
	QAWS_SUPPORT_LOCAL = 0,   /* a few contiguous elements per field */
	QAWS_SUPPORT_WIDE,        /* many elements, still bounded */
	QAWS_SUPPORT_GLOBAL,      /* every element may contribute */
	QAWS_SUPPORT_IMPLICIT     /* defined through a solve */
} qaws_support_kind;

#define QAWS_DIFF_MAX_ORDER 5
#define QAWS_DIFF_MAX_SUPPORT 32
#define QAWS_DIFF_MAX_RANGES 4

typedef struct qaws_support_range
{
	qaws_diff_field field;
	unsigned int first;
	unsigned int count;
} qaws_support_range;

typedef struct qaws_local_support
{
	qaws_support_kind kind;
	unsigned int range_count;
	qaws_support_range ranges[QAWS_DIFF_MAX_RANGES];
	int has_weights;
	unsigned int order;
	qaws_scalar weights[QAWS_DIFF_MAX_RANGES][QAWS_DIFF_MAX_ORDER + 1][QAWS_DIFF_MAX_SUPPORT];
} qaws_local_support;

/* Tensor-product surface support: weight(i,j) = u_weights[a][i] * v_weights[b][j]
   for the partial d^(a+b) S / du^a dv^b. Elements are addressed as
   (u_first + i) * v_stride + (v_first + j). */
typedef struct qaws_surface_local_support
{
	qaws_support_kind kind;
	qaws_diff_field field;
	unsigned int u_first, u_count;
	unsigned int v_first, v_count;
	unsigned int v_stride;
	int has_weights;
	unsigned int order;
	qaws_scalar u_weights[QAWS_DIFF_MAX_ORDER + 1][QAWS_DIFF_MAX_SUPPORT];
	qaws_scalar v_weights[QAWS_DIFF_MAX_ORDER + 1][QAWS_DIFF_MAX_SUPPORT];
} qaws_surface_local_support;

/* ===================================================================
 * Context and report
 * =================================================================== */

typedef struct qaws_diff_report
{
	qaws_diff_class diff_class;       /* worst class touched */
	qaws_diff_validity validity;      /* worst runtime validity */
	unsigned int frozen_used;         /* qaws_diff_freeze bits the result relies on */
	qaws_scalar condition_number;     /* largest seen, 1 when not applicable */
	qaws_scalar residual;             /* largest implicit residual */
	qaws_scalar branch_gap;           /* smallest gap to a competing solution */
	unsigned int worst_index;         /* sample index of the worst validity */
	unsigned int evaluation_count;
} qaws_diff_report;

typedef struct qaws_diff_context
{
	unsigned int order;                     /* 1 = first order, 2 = also second */
	unsigned int frozen;                    /* qaws_diff_freeze bits the caller accepts */
	qaws_diff_accumulation accumulation;
	unsigned int tile_size;                 /* samples per tile, 0 = default */
	qaws_allocator const* allocator;        /* scratch allocations, NULL = malloc */
	qaws_diff_report* report;               /* optional output */
} qaws_diff_context;

#endif /* QAWS_DIFF_TYPES_H */
