/*
 * Test 181: attributes through Booleans and sweeps (qaws_mesh_3d_boolean,
 * qaws_mesh_3d_sweep, qaws_mesh_sweep_boolean), and the colour rules
 *
 *   - colour spaces on an edge split: black to white in Oklab (the midpoint
 *     Oklab's L = 0.5: linear light 0.125, not sRGB 0.5); transparent black to red in sRGB with
 *     alpha (premultiplied: the midpoint is red at alpha 0.5, not dark);
 *   - a unit cube less a cube shifted by half (volume 0.875): on A's faces
 *     an affine vertex field, a corner UV, an sRGB colour linear in light
 *     (so encode(x, y, z) exactly at every corner, the cut ones included),
 *     on B's their own field; a face material on both, a face NORMAL on A
 *     alone (B's faces take theirs, turned with the carved pieces): every
 *     corner exact, every normal along its face;
 *   - a cube swept along x carrying a field constant along the motion (each
 *     corner's value its own y + 2 z) and face normals (walls' their own);
 *     a ball's radial vertex normals along a segment, against the faces as
 *     well as on the ball itself (the walls' normals its silhouette's);
 *   - carving: the slot's faces carry the tool's material, their area the
 *     slot's (floor and walls, 4), all normals along their faces;
 *   - edge attributes through edits: a split edge's halves keep its value,
 *     the new edges zeros; a flipped edge's new diagonal zeros, the others
 *     kept;
 *   - edge attributes and features through A - B: a crease weight on A's
 *     box edges stays on their remaining pieces (length 10.5); the result's
 *     features are exactly its creases (A's and B's kept, the cut curves
 *     added: length 15);
 *   - tangents with handedness: unit directions, the sign +1 on A's pieces,
 *     -1 on B's turned ones;
 *   - vectors and orientations: an edge split mixes VECTORs linearly (the
 *     midpoint of (2, 0, 0) and (0, 2, 0) is (1, 1, 0), not unit) and
 *     QUATERNIONs on one side (the identity and minus a quarter turn about z
 *     give the eighth turn, not a near zero); a Boolean (A's vectors per
 *     vertex, B's per corner: matched) leaves B's turned
 *     pieces' vectors as they were (a NORMAL would be negated); a box swept
 *     through a quarter turn about z: at the end pose's outer corners the vector
 *     turned to (0, 1, 0) and the quaternion the quarter turn's, at the
 *     start's both as given;
 *   - through subdivision: an ANGLE of 3 and -3 radians meets at pi (the
 *     short way round), not 0; QUATERNIONs of 170 and -170 degrees about z
 *     meet at the half turn (sign-free, no flip through the identity); edge
 *     values on a quad's edges kept by their pieces after two Catmull-Clark
 *     levels (4 each, the inner edges 0), on two triangles' edges after one
 *     Loop level (2 each);
 *   - features along a sweep's walls: a cube (its creases features) swept
 *     along (1, 0.5, 0): the result's features are exactly its creases, the
 *     long edges its corners trace included; along (1, 0, 0) (its x edges
 *     sliding along themselves): 16 of each; a 16-sided prism swept along
 *     (1, 0, 0): its features exactly its creases past 0.5 radians, the
 *     lines its rims trace where they run along the motion included;
 *     swept along a curve while
 *     turning a quarter (16 slabs): closed, one piece (no sliver left between
 *     two slabs' twisted walls), its features exactly its edges turning by
 *     more than 25 degrees (the slabs' chords bend less), and the cube at
 *     101 times along the motion inside it (points 0.01 within its sides, at
 *     mid height and by the top: nothing a turning face sweeps left out);
 *   - a sweep simplified (qaws_sweep_desc.simplify): a featured cube with
 *     a LINEAR field, per-face UVs (seams) and edge values, along a curve
 *     while turning a quarter: every corner's values the plain sweep's at
 *     its place (a point just inside its triangle), the edge values'
 *     lengths and the features as the plain sweep's (the collapses keep
 *     each wedge's values, the lines carrying values and the features), on
 *     far fewer faces;
 *   - a face sliding along itself: a cube with its face's number on each face's corners, along a curve while turning a
 *     quarter about z: every upward face of the sweep (its top, which the cube's top face sweeps in its own plane) carries
 *     the top's number, none a side's (the side faces' hulls' top walls lie in that plane: they take the top face's values
 *     at the edge that swept them);
 *   - edge values through the connector: two tubes joined, each tube's edges
 *     their values (a's 1 + i, b's 100 + i), the bridge's 0;
 *   - refused: one name of two elements, a NORMAL of 2 components.
 */

#include "test_common.h"
#include "qaws_mesh_sweep.h"
#include "qaws_mesh_subdivide.h"
#include "qaws_mesh_connect.h"
#include "qaws_mesh_bvh.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static double ac_encode(double l)
{
	return l <= 0.0031308 ? 12.92 * l : 1.055 * pow(l, 1 / 2.4) - 0.055;
}

/* the box [x0, x1] x [y0, y1] x [z0, z1], outward */
static qaws_mesh_3d* ac_box(double x0, double x1, double y0, double y1, double z0, double z1)
{
	static unsigned int const f[36] = { 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3 };
	qaws_scalar x[24];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	unsigned int i;
	for (i = 0; i < 8; i++)
	{
		x[3 * i] = (qaws_scalar)((i == 1 || i == 2 || i == 5 || i == 6) ? x1 : x0);
		x[3 * i + 1] = (qaws_scalar)((i == 2 || i == 3 || i == 6 || i == 7) ? y1 : y0);
		x[3 * i + 2] = (qaws_scalar)(i >= 4 ? z1 : z0);
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 8; d.face_vertices = f; d.face_count = 12;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

static unsigned int ac_add(qaws_mesh_3d* m, char const* name, qaws_mesh_element e, unsigned int n, qaws_attribute_rule rule, qaws_scalar const* v)
{
	qaws_attribute_desc d;
	unsigned int id = QAWS_MESH_NONE;
	d.name = name;
	d.element = e;
	d.components = n;
	d.rule = rule;
	qaws_mesh_3d_add_attribute(&id, m, v, &d);
	return id;
}

/* per vertex: f(x, y, z) = a x + b y + c z + k */
static void ac_vertex_field(qaws_mesh_3d* m, char const* name, double a, double b, double c, double k)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), i;
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nv);
	for (i = 0; i < nv; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
		v[i] = (qaws_scalar)(a * p.x + b * p.y + c * p.z + k);
	}
	ac_add(m, name, QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_LINEAR, v);
	free(v);
}

/* a triangle's unit normal */
static void ac_tri_normal(qaws_mesh_3d const* m, unsigned int f, double* n)
{
	unsigned int c[3], k;
	qaws_vec3 a, b, e;
	double u[3], w[3], l;
	qaws_mesh_3d_get_triangle(c, m, f);
	a = qaws_mesh_3d_get_vertex(m, c[0]);
	b = qaws_mesh_3d_get_vertex(m, c[1]);
	e = qaws_mesh_3d_get_vertex(m, c[2]);
	u[0] = (double)b.x - a.x; u[1] = (double)b.y - a.y; u[2] = (double)b.z - a.z;
	w[0] = (double)e.x - a.x; w[1] = (double)e.y - a.y; w[2] = (double)e.z - a.z;
	n[0] = u[1] * w[2] - u[2] * w[1];
	n[1] = u[2] * w[0] - u[0] * w[2];
	n[2] = u[0] * w[1] - u[1] * w[0];
	l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
	for (k = 0; k < 3; k++) n[k] = l > 0 ? n[k] / l : 0;
}

/* a face NORMAL attribute: each triangle's own normal */
static void ac_face_normals(qaws_mesh_3d* m)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), i;
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nf);
	for (i = 0; i < nf; i++)
	{
		double n[3];
		ac_tri_normal(m, i, n);
		v[3 * i] = (qaws_scalar)n[0]; v[3 * i + 1] = (qaws_scalar)n[1]; v[3 * i + 2] = (qaws_scalar)n[2];
	}
	ac_add(m, "nrm", QAWS_MESH_FACE, 3, QAWS_ATTRIBUTE_NORMAL, v);
	free(v);
}

static void ac_face_constant(qaws_mesh_3d* m, char const* name, double value)
{
	unsigned int nf = qaws_mesh_3d_get_face_count(m), i;
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * nf);
	for (i = 0; i < nf; i++) v[i] = (qaws_scalar)value;
	ac_add(m, name, QAWS_MESH_FACE, 1, QAWS_ATTRIBUTE_NEAREST, v);
	free(v);
}

/* an attribute's values (NULL if absent) and element */
static qaws_scalar* ac_get(qaws_mesh_3d const* m, char const* name, qaws_mesh_element* e, unsigned int* n)
{
	unsigned int id;
	qaws_attribute_desc d;
	qaws_scalar* v;
	unsigned int count;
	if (qaws_mesh_3d_find_attribute(&id, m, name) != QAWS_STATUS_OK)
		return NULL;
	qaws_mesh_3d_get_attribute_desc(&d, m, id);
	count = d.element == QAWS_MESH_VERTEX ? qaws_mesh_3d_get_vertex_count(m) : d.element == QAWS_MESH_CORNER ? qaws_mesh_3d_get_corner_count(m) : qaws_mesh_3d_get_face_count(m);
	v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * d.components * count + 1);
	qaws_mesh_3d_get_attribute(v, m, id);
	*e = d.element;
	*n = d.components;
	return v;
}

/* the largest misalignment of a face NORMAL attribute with the faces: 1 - min dot */
static double ac_normal_error(qaws_mesh_3d const* m)
{
	qaws_mesh_element e;
	unsigned int n, nf = qaws_mesh_3d_get_face_count(m), i;
	qaws_scalar* v = ac_get(m, "nrm", &e, &n);
	double worst = 0;
	if (!v || e != QAWS_MESH_FACE)
	{
		free(v);
		return 2;
	}
	for (i = 0; i < nf; i++)
	{
		double g[3];
		ac_tri_normal(m, i, g);
		worst = fmax(worst, 1 - (g[0] * v[3 * i] + g[1] * v[3 * i + 1] + g[2] * v[3 * i + 2]));
	}
	free(v);
	return worst;
}

static void test_colour_spaces(void)
{
	static unsigned int const t[3] = { 0, 1, 2 };
	qaws_scalar x[9] = { 0, 0, 0, 1, 0, 0, 0, 1, 0 };
	qaws_scalar ok[9] = { 0, 0, 0, 1, 1, 1, 0, 0, 0 }, sr[12] = { 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1 }, o[4];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	unsigned int v = 0, a0, a1;
	double l = ac_encode(0.125), e0, e1;
	char msg[300];
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 3; d.face_vertices = t; d.face_count = 1;
	qaws_mesh_3d_create(&m, &d);
	a0 = ac_add(m, "oklab", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_OKLAB, ok);
	a1 = ac_add(m, "srgba", QAWS_MESH_VERTEX, 4, QAWS_ATTRIBUTE_SRGB, sr);
	qaws_mesh_3d_split_edge(&v, m, NULL, NULL, 0, 1);
	{
		qaws_scalar* all = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 4 * qaws_mesh_3d_get_vertex_count(m));
		qaws_mesh_3d_get_attribute(all, m, a0);
		e0 = fabs((double)all[3 * v] - l) + fabs((double)all[3 * v + 1] - l) + fabs((double)all[3 * v + 2] - l);
		qaws_mesh_3d_get_attribute(all, m, a1);
		memcpy(o, all + 4 * v, sizeof o);
		e1 = fabs((double)o[0] - 1) + fabs((double)o[1]) + fabs((double)o[2]) + fabs((double)o[3] - 0.5);
		free(all);
	}
	sprintf(msg, "an edge split between black and white in Oklab: %.6f (%.6f, L = 0.5), error %.1e; transparent to red in sRGB: (%.3f %.3f %.3f %.3f), error %.1e",
		l, l, e0, (double)o[0], (double)o[1], (double)o[2], (double)o[3], e1);
	printf("  %s\n", msg);
	TEST_ASSERT(e0 < 1e-9 && e1 < 1e-12, msg);
	qaws_mesh_3d_destroy(m);
}

static void test_boolean(void)
{
	qaws_mesh_3d *a = ac_box(0, 1, 0, 1, 0, 1), *b = ac_box(0.5, 1.5, 0.5, 1.5, 0.5, 1.5), *o = NULL;
	qaws_boolean_report r;
	qaws_status st;
	unsigned int i, j, nf = 0, from_a = 0, from_b = 0, nfield = 0, nuv = 0, ncol = 0, nmat = 0, nsrc = 0;
	double ef = 0, eu = 0, ec = 0, en, vol = 0;
	qaws_mesh_element ee[5] = { QAWS_MESH_VERTEX, QAWS_MESH_VERTEX, QAWS_MESH_VERTEX, QAWS_MESH_VERTEX, QAWS_MESH_VERTEX };
	qaws_scalar *field = NULL, *uv = NULL, *col = NULL, *mat = NULL, *src = NULL;
	char msg[600];
	/* A: field x + 2y + 3z, uv (x + y, z) per corner, colour encode(x, y, z), material 1, its normals; B: field 10 + x, material 2 */
	ac_vertex_field(a, "field", 1, 2, 3, 0);
	ac_vertex_field(b, "field", 1, 0, 0, 10);
	{
		unsigned int nc = qaws_mesh_3d_get_corner_count(a), nv = qaws_mesh_3d_get_vertex_count(a);
		unsigned int* fv = (unsigned int*)malloc(sizeof(unsigned int) * nc);
		qaws_scalar* u = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nc > nv ? nc : nv));
		qaws_mesh_3d_get_faces(fv, NULL, a);
		for (i = 0; i < nc; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(a, fv[i]);
			u[2 * i] = p.x + p.y;
			u[2 * i + 1] = p.z;
		}
		ac_add(a, "uv", QAWS_MESH_CORNER, 2, QAWS_ATTRIBUTE_LINEAR, u);
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(a, i);
			u[3 * i] = (qaws_scalar)ac_encode(p.x);
			u[3 * i + 1] = (qaws_scalar)ac_encode(p.y);
			u[3 * i + 2] = (qaws_scalar)ac_encode(p.z);
		}
		ac_add(a, "color", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_SRGB, u);
		free(fv);
		free(u);
	}
	ac_face_constant(a, "material", 1);
	ac_face_constant(b, "material", 2);
	ac_face_normals(a);
	st = qaws_mesh_3d_boolean(&o, &r, a, b, QAWS_BOOLEAN_DIFFERENCE, NULL);
	if (st == QAWS_STATUS_OK && o)
	{
		field = ac_get(o, "field", &ee[0], &nfield);
		uv = ac_get(o, "uv", &ee[1], &nuv);
		col = ac_get(o, "color", &ee[2], &ncol);
		mat = ac_get(o, "material", &ee[3], &nmat);
		src = ac_get(o, "source", &ee[4], &nsrc);
		nf = qaws_mesh_3d_get_face_count(o);
	}
	for (i = 0; field && uv && col && mat && src && i < nf; i++)
	{
		unsigned int c[3];
		int in_a = (unsigned int)src[i] < 12;
		qaws_vec3 p[3];
		qaws_mesh_3d_get_triangle(c, o, i);
		for (j = 0; j < 3; j++) p[j] = qaws_mesh_3d_get_vertex(o, c[j]);
		vol += ((double)p[0].x * ((double)p[1].y * p[2].z - (double)p[1].z * p[2].y) + (double)p[0].y * ((double)p[1].z * p[2].x - (double)p[1].x * p[2].z) +
			(double)p[0].z * ((double)p[1].x * p[2].y - (double)p[1].y * p[2].x)) / 6;
		from_a += (unsigned int)in_a;
		from_b += (unsigned int)!in_a;
		ef = fmax(ef, fabs((double)mat[i] - (in_a ? 1 : 2)));
		for (j = 0; j < 3; j++)
		{
			double x = p[j].x, y = p[j].y, z = p[j].z;
			unsigned int k = 3 * i + j;
			ef = fmax(ef, fabs((double)field[k] - (in_a ? x + 2 * y + 3 * z : 10 + x)));
			eu = fmax(eu, fabs((double)uv[2 * k] - (in_a ? x + y : 0)) + fabs((double)uv[2 * k + 1] - (in_a ? z : 0)));
			ec = fmax(ec, fabs((double)col[3 * k] - (in_a ? ac_encode(x) : 0)) + fabs((double)col[3 * k + 1] - (in_a ? ac_encode(y) : 0)) +
				fabs((double)col[3 * k + 2] - (in_a ? ac_encode(z) : 0)));
		}
	}
	en = o ? ac_normal_error(o) : 2;
	sprintf(msg, "A - B with attributes: volume %.12f (0.875), %u faces (%u of A, %u of B); errors: field and material %.1e, uv %.1e, sRGB colour %.1e, normals 1 - dot %.1e; field per %s, material per %s",
		vol, nf, from_a, from_b, ef, eu, ec, en, ee[0] == QAWS_MESH_CORNER ? "corner" : "?", ee[3] == QAWS_MESH_FACE ? "face" : "?");
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && fabs(vol - 0.875) < 1e-12 && from_a && from_b && ef < 1e-9 && eu < 1e-9 && ec < 1e-9 && en < 1e-9 &&
		ee[0] == QAWS_MESH_CORNER && ee[1] == QAWS_MESH_CORNER && ee[2] == QAWS_MESH_CORNER && ee[3] == QAWS_MESH_FACE && r.closed, msg);
	free(field); free(uv); free(col); free(mat); free(src);
	qaws_mesh_3d_destroy(a); qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o);
}

static qaws_motion* ac_move(double x0, double y0, double z0, double x1, double y1, double z1)
{
	qaws_motion_pose s, e;
	qaws_motion* m = NULL;
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = e.rotation[0] = 1;
	s.position.x = (qaws_scalar)x0; s.position.y = (qaws_scalar)y0; s.position.z = (qaws_scalar)z0;
	e.position.x = (qaws_scalar)x1; e.position.y = (qaws_scalar)y1; e.position.z = (qaws_scalar)z1;
	s.velocity.x = e.velocity.x = (qaws_scalar)(x1 - x0);
	s.velocity.y = e.velocity.y = (qaws_scalar)(y1 - y0);
	s.velocity.z = e.velocity.z = (qaws_scalar)(z1 - z0);
	qaws_motion_create_poses(&m, &s, &e, 0);
	return m;
}

/* an icosahedron split once, on the sphere of radius r: 80 faces */
static qaws_mesh_3d* ac_ball(double r)
{
	static int const ic[60] = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8, 3, 9, 4, 3, 4, 2,
		3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
	double p = (1 + sqrt(5.0)) / 2, x[3 * 42];
	double const v0[36] = { -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, 0, 0, -1, p, 0, 1, p, 0, -1, -p, 0, 1, -p, p, 0, -1, p, 0, 1, -p, 0, -1, -p, 0, 1 };
	unsigned int t[3 * 80], nv = 12, nt = 0, i, j, mid[30][3], nm = 0;
	qaws_scalar xs[3 * 42];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < 36; i++) x[i] = v0[i];
	for (i = 0; i < 20; i++)
	{
		unsigned int c[3], e;
		for (e = 0; e < 3; e++)
		{
			unsigned int a = (unsigned int)ic[3 * i + e], b = (unsigned int)ic[3 * i + (e + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a, k;
			for (k = 0; k < nm && !(mid[k][0] == lo && mid[k][1] == hi); k++) {}
			if (k == nm)
			{
				mid[nm][0] = lo; mid[nm][1] = hi; mid[nm][2] = nv;
				for (j = 0; j < 3; j++) x[3 * nv + j] = 0.5 * (x[3 * a + j] + x[3 * b + j]);
				nv++; nm++;
			}
			c[e] = mid[k][2];
		}
		t[3 * nt] = (unsigned int)ic[3 * i]; t[3 * nt + 1] = c[0]; t[3 * nt + 2] = c[2]; nt++;
		t[3 * nt] = (unsigned int)ic[3 * i + 1]; t[3 * nt + 1] = c[1]; t[3 * nt + 2] = c[0]; nt++;
		t[3 * nt] = (unsigned int)ic[3 * i + 2]; t[3 * nt + 1] = c[2]; t[3 * nt + 2] = c[1]; nt++;
		t[3 * nt] = c[0]; t[3 * nt + 1] = c[1]; t[3 * nt + 2] = c[2]; nt++;
	}
	for (i = 0; i < nv; i++)
	{
		double l = sqrt(x[3 * i] * x[3 * i] + x[3 * i + 1] * x[3 * i + 1] + x[3 * i + 2] * x[3 * i + 2]);
		for (j = 0; j < 3; j++) xs[3 * i + j] = (qaws_scalar)(r * x[3 * i + j] / l);
	}
	memset(&d, 0, sizeof d);
	d.positions = xs; d.vertex_count = nv; d.face_vertices = t; d.face_count = nt;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

static void test_sweep(void)
{
	qaws_mesh_3d *c = ac_box(0, 1, 0, 1, 0, 1), *o = NULL, *ball = ac_ball(0.3), *ob = NULL;
	qaws_motion *m = ac_move(0, 0, 0, 1, 0, 0), *mb = ac_move(0, 0, 0, 1, 0.5, 0);
	qaws_sweep_report r, rb;
	qaws_status st, sb;
	qaws_mesh_element e = QAWS_MESH_VERTEX, eb = QAWS_MESH_VERTEX;
	unsigned int n, nf = 0, i, j;
	double ef = 0, en, worst = 1, mean = 0, worst0 = 1, mean0 = 0;
	qaws_scalar *field = NULL, *bn = NULL;
	char msg[500];
	ac_vertex_field(c, "field", 0, 1, 2, 0);
	ac_face_normals(c);
	{
		/* the ball's radial normals per vertex */
		unsigned int nv = qaws_mesh_3d_get_vertex_count(ball);
		qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * nv);
		for (i = 0; i < nv; i++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(ball, i);
			v[3 * i] = p.x / (qaws_scalar)0.3; v[3 * i + 1] = p.y / (qaws_scalar)0.3; v[3 * i + 2] = p.z / (qaws_scalar)0.3;
		}
		ac_add(ball, "n", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_NORMAL, v);
		/* (on the ball itself: its vertex normals against its faces) */
		for (i = 0; i < qaws_mesh_3d_get_face_count(ball); i++)
		{
			unsigned int t[3];
			double g[3];
			qaws_mesh_3d_get_triangle(t, ball, i);
			ac_tri_normal(ball, i, g);
			for (j = 0; j < 3; j++)
			{
				double dt = g[0] * v[3 * t[j]] + g[1] * v[3 * t[j] + 1] + g[2] * v[3 * t[j] + 2];
				worst0 = fmin(worst0, dt);
				mean0 += dt / (3.0 * qaws_mesh_3d_get_face_count(ball));
			}
		}
		free(v);
	}
	st = qaws_mesh_3d_sweep(&o, &r, c, m, NULL);
	sb = qaws_mesh_3d_sweep(&ob, &rb, ball, mb, NULL);
	if (st == QAWS_STATUS_OK && o)
	{
		field = ac_get(o, "field", &e, &n);
		nf = qaws_mesh_3d_get_face_count(o);
	}
	for (i = 0; field && i < nf; i++)
	{
		unsigned int t[3];
		qaws_mesh_3d_get_triangle(t, o, i);
		for (j = 0; j < 3; j++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(o, t[j]);
			ef = fmax(ef, fabs((double)field[3 * i + j] - ((double)p.y + 2 * (double)p.z)));
		}
	}
	en = o ? ac_normal_error(o) : 2;
	if (sb == QAWS_STATUS_OK && ob && (bn = ac_get(ob, "n", &eb, &n)) != NULL)
	{
		unsigned int nfb = qaws_mesh_3d_get_face_count(ob);
		for (i = 0; i < nfb; i++)
		{
			double g[3];
			ac_tri_normal(ob, i, g);
			for (j = 0; j < 3; j++)
			{
				double dt = g[0] * bn[9 * i + 3 * j] + g[1] * bn[9 * i + 3 * j + 1] + g[2] * bn[9 * i + 3 * j + 2];
				worst = fmin(worst, dt);
				mean += dt / (3.0 * nfb);
			}
		}
	}
	sprintf(msg, "a cube swept along x: field error %.1e (per %s), face normals 1 - dot %.1e, closed %d; a ball's vertex normals along a segment, per %s: dot with the faces mean %.4f, least %.4f (on the ball %.4f, %.4f)",
		ef, e == QAWS_MESH_CORNER ? "corner" : "?", en, r.closed, eb == QAWS_MESH_CORNER ? "corner" : "?", mean, worst, mean0, worst0);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && sb == QAWS_STATUS_OK && field && bn && e == QAWS_MESH_CORNER && eb == QAWS_MESH_CORNER && ef < 1e-9 && en < 1e-9 &&
		mean > mean0 - 0.01 && worst > worst0 - 0.05 && r.closed && rb.closed, msg);
	free(field); free(bn);
	qaws_mesh_3d_destroy(c); qaws_mesh_3d_destroy(o); qaws_mesh_3d_destroy(ball); qaws_mesh_3d_destroy(ob);
	qaws_motion_destroy(m); qaws_motion_destroy(mb);
}

static void test_carving(void)
{
	qaws_mesh_3d *block = ac_box(0, 4, 0, 2, 0, 1), *tool = ac_box(-0.25, 0.25, -0.25, 0.25, -0.25, 0.25), *o = NULL;
	qaws_motion* slot = ac_move(-1, 1, 1, 5, 1, 1);
	qaws_motion const* one[1];
	qaws_sweep_boolean_report r;
	qaws_status st;
	qaws_mesh_element e = QAWS_MESH_VERTEX;
	unsigned int n, nf = 0, i;
	double area2 = 0, en;
	qaws_scalar* mat = NULL;
	char msg[400];
	one[0] = slot;
	ac_face_constant(block, "material", 1);
	ac_face_constant(tool, "material", 2);
	ac_face_normals(block);
	ac_face_normals(tool);
	st = qaws_mesh_sweep_boolean(&o, &r, tool, one, block, 1, QAWS_BOOLEAN_DIFFERENCE, NULL);
	if (st == QAWS_STATUS_OK && o)
	{
		mat = ac_get(o, "material", &e, &n);
		nf = qaws_mesh_3d_get_face_count(o);
	}
	for (i = 0; mat && i < nf; i++)
		if (mat[i] == 2)
		{
			unsigned int t[3];
			qaws_vec3 a, b, c;
			double u[3], w[3];
			qaws_mesh_3d_get_triangle(t, o, i);
			a = qaws_mesh_3d_get_vertex(o, t[0]);
			b = qaws_mesh_3d_get_vertex(o, t[1]);
			c = qaws_mesh_3d_get_vertex(o, t[2]);
			u[0] = (double)b.x - a.x; u[1] = (double)b.y - a.y; u[2] = (double)b.z - a.z;
			w[0] = (double)c.x - a.x; w[1] = (double)c.y - a.y; w[2] = (double)c.z - a.z;
			area2 += 0.5 * sqrt(pow(u[1] * w[2] - u[2] * w[1], 2) + pow(u[2] * w[0] - u[0] * w[2], 2) + pow(u[0] * w[1] - u[1] * w[0], 2));
		}
	en = o ? ac_normal_error(o) : 2;
	sprintf(msg, "a slot carved: the tool's material on %.12f of the surface (4: the floor and two walls), normals 1 - dot %.1e, removed %.12f, closed %d",
		area2, en, (double)r.changed_volume, r.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && mat && e == QAWS_MESH_FACE && fabs(area2 - 4) < 1e-9 && en < 1e-9 && fabs((double)r.changed_volume - 0.5) < 1e-12 && r.closed, msg);
	free(mat);
	qaws_mesh_3d_destroy(block); qaws_mesh_3d_destroy(tool); qaws_mesh_3d_destroy(o);
	qaws_motion_destroy(slot);
}

/* an edge attribute's value on edge (a, b), or -1 */
static double ac_edge_value(qaws_mesh_3d const* m, unsigned int id, unsigned int a, unsigned int b)
{
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), i;
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * ne);
	double r = -1;
	qaws_mesh_3d_get_edges(e, m);
	qaws_mesh_3d_get_attribute(v, m, id);
	for (i = 0; i < ne; i++)
		if ((e[2 * i] == a && e[2 * i + 1] == b) || (e[2 * i] == b && e[2 * i + 1] == a))
			r = (double)v[i];
	free(e);
	free(v);
	return r;
}

static void test_edge_edits(void)
{
	/* the square 0 1 2 3 as two triangles over the diagonal 0 2; each edge's value 1 + its index in this list */
	static unsigned int const t[6] = { 0, 1, 2, 0, 2, 3 };
	static unsigned int const pairs[10] = { 0, 1, 1, 2, 2, 3, 3, 0, 0, 2 };
	qaws_scalar x[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_attribute_desc ad;
	unsigned int id = 0, v = 0, ne, i, j;
	unsigned int e[16];
	qaws_scalar vals[8];
	double h0, h1, cross, flipped, kept;
	char msg[300];
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 4; d.face_vertices = t; d.face_count = 2;
	qaws_mesh_3d_create(&m, &d);
	ne = qaws_mesh_3d_get_edge_count(m);
	qaws_mesh_3d_get_edges(e, m);
	for (i = 0; i < ne; i++)
		for (j = 0; j < 5; j++)
			if ((e[2 * i] == pairs[2 * j] && e[2 * i + 1] == pairs[2 * j + 1]) || (e[2 * i] == pairs[2 * j + 1] && e[2 * i + 1] == pairs[2 * j]))
				vals[i] = (qaws_scalar)(1 + j);
	ad.name = "crease"; ad.element = QAWS_MESH_EDGE; ad.components = 1; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&id, m, vals, &ad);
	/* the edge 0 1 split: its halves 0 v and v 1 keep 1; the new edge v 2 zero */
	qaws_mesh_3d_split_edge(&v, m, NULL, NULL, 0, 1);
	h0 = ac_edge_value(m, id, 0, v);
	h1 = ac_edge_value(m, id, v, 1);
	cross = ac_edge_value(m, id, v, 2);
	/* the diagonal 0 2 flipped: the new one v 3 zero, the border kept */
	qaws_mesh_3d_flip_edge(m, 0, 2);
	flipped = ac_edge_value(m, id, v, 3);
	kept = ac_edge_value(m, id, 2, 3) + ac_edge_value(m, id, 3, 0) + ac_edge_value(m, id, 1, 2);
	sprintf(msg, "edge values through edits: the split edge's halves %.0f %.0f (1 1), the new edge %.0f (0); after a flip the new diagonal %.0f (0), the border kept %.0f (3 + 4 + 2 = 9)",
		h0, h1, cross, flipped, kept);
	printf("  %s\n", msg);
	TEST_ASSERT(h0 == 1 && h1 == 1 && cross == 0 && flipped == 0 && kept == 9, msg);
	qaws_mesh_3d_destroy(m);
}

static void test_edge_boolean(void)
{
	qaws_mesh_3d *a = ac_box(0, 1, 0, 1, 0, 1), *b = ac_box(0.5, 1.5, 0.5, 1.5, 0.5, 1.5), *o = NULL;
	qaws_status st;
	unsigned int ne = 0, nfe = 0, id = 0, mism = 0;
	double crease = 0, feat = 0, tan_err = 0;
	int signs_ok = 1;
	char msg[400];
	/* A: crease 1 on its box edges, 0 on the diagonals; both: features from their 90-degree edges; tangents (x, +1) on A, (y, +1) on B */
	qaws_mesh_3d_detect_features(a, (qaws_scalar)0.5, 0, 0);
	qaws_mesh_3d_detect_features(b, (qaws_scalar)0.5, 0, 0);
	{
		unsigned int n = qaws_mesh_3d_get_edge_count(a), k;
		unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * n);
		qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * n);
		qaws_scalar ta[32], tb[32];
		qaws_mesh_3d_get_edges(e, a);
		for (k = 0; k < n; k++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(a, e[2 * k]), q = qaws_mesh_3d_get_vertex(a, e[2 * k + 1]);
			int same = (p.x == q.x) + (p.y == q.y) + (p.z == q.z);
			v[k] = (qaws_scalar)(same == 2 ? 1 : 0);
		}
		ac_add(a, "crease", QAWS_MESH_EDGE, 1, QAWS_ATTRIBUTE_LINEAR, v);
		for (k = 0; k < 8; k++)
		{
			ta[4 * k] = 1; ta[4 * k + 1] = 0; ta[4 * k + 2] = 0; ta[4 * k + 3] = 1;
			tb[4 * k] = 0; tb[4 * k + 1] = 1; tb[4 * k + 2] = 0; tb[4 * k + 3] = 1;
		}
		ac_add(a, "tangent", QAWS_MESH_VERTEX, 4, QAWS_ATTRIBUTE_TANGENT, ta);
		ac_add(b, "tangent", QAWS_MESH_VERTEX, 4, QAWS_ATTRIBUTE_TANGENT, tb);
		free(e);
		free(v);
	}
	st = qaws_mesh_3d_boolean(&o, NULL, a, b, QAWS_BOOLEAN_DIFFERENCE, NULL);
	if (st == QAWS_STATUS_OK && o && qaws_mesh_3d_find_attribute(&id, o, "crease") == QAWS_STATUS_OK)
	{
		unsigned int nf = qaws_mesh_3d_get_face_count(o), nfeat = qaws_mesh_3d_get_feature_edge_count(o), f, j, k;
		unsigned int *e, *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfeat + 1)), *t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
		qaws_scalar *v, *tg = NULL, *src = NULL;
		qaws_mesh_element el;
		unsigned int nc;
		double* fn = (double*)malloc(sizeof(double) * 3 * nf);
		ne = qaws_mesh_3d_get_edge_count(o);
		e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
		v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * ne);
		qaws_mesh_3d_get_edges(e, o);
		qaws_mesh_3d_get_attribute(v, o, id);
		qaws_mesh_3d_get_feature_edges(fe, o);
		qaws_mesh_3d_get_triangles(t, o);
		for (f = 0; f < nf; f++) ac_tri_normal(o, f, fn + 3 * f);
		nfe = nfeat;
		for (k = 0; k < ne; k++)
		{
			qaws_vec3 p = qaws_mesh_3d_get_vertex(o, e[2 * k]), q = qaws_mesh_3d_get_vertex(o, e[2 * k + 1]);
			double len = sqrt(pow((double)p.x - q.x, 2) + pow((double)p.y - q.y, 2) + pow((double)p.z - q.z, 2)), turn = 0;
			int is_feat = 0, side = 0;
			unsigned int g[2] = { 0, 0 };
			for (j = 0; j < nfeat; j++)
				if ((fe[2 * j] == e[2 * k] && fe[2 * j + 1] == e[2 * k + 1]) || (fe[2 * j] == e[2 * k + 1] && fe[2 * j + 1] == e[2 * k]))
					is_feat = 1;
			for (f = 0; f < nf && side < 2; f++)
				for (j = 0; j < 3; j++)
					if ((t[3 * f + j] == e[2 * k] && t[3 * f + (j + 1) % 3] == e[2 * k + 1]) || (t[3 * f + j] == e[2 * k + 1] && t[3 * f + (j + 1) % 3] == e[2 * k]))
						g[side++] = f;
			if (side == 2)
				turn = 1 - (fn[3 * g[0]] * fn[3 * g[1]] + fn[3 * g[0] + 1] * fn[3 * g[1] + 1] + fn[3 * g[0] + 2] * fn[3 * g[1] + 2]);
			crease += (double)v[k] * len;
			feat += is_feat * len;
			mism += (unsigned int)((turn > 1e-9) != is_feat);
		}
		/* tangents: A's pieces (x, +1), B's turned ones (y, -1) */
		tg = ac_get(o, "tangent", &el, &nc);
		src = ac_get(o, "source", &el, &nc);
		for (f = 0; tg && src && f < nf; f++)
			for (j = 0; j < 3; j++)
			{
				qaws_scalar const* q = tg + 4 * (3 * f + j);
				int in_a = (unsigned int)src[f] < 12;
				tan_err = fmax(tan_err, fabs((double)q[0] - (in_a ? 1 : 0)) + fabs((double)q[1] - (in_a ? 0 : 1)) + fabs((double)q[2]));
				signs_ok &= (double)q[3] == (in_a ? 1 : -1);
			}
		free(e); free(v); free(fe); free(t); free(fn); free(tg); free(src);
	}
	sprintf(msg, "A - B: the crease weight along %.12f (10.5: A's box edges left), %u features along %.12f (15: every crease, the cuts added), %u edges a feature without a crease or the reverse; tangents error %.1e, handedness %s",
		crease, nfe, feat, mism, tan_err, signs_ok ? "+1 on A, -1 on B's turned pieces" : "wrong");
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && fabs(crease - 10.5) < 1e-9 && fabs(feat - 15) < 1e-9 && mism == 0 && tan_err < 1e-12 && signs_ok, msg);
	qaws_mesh_3d_destroy(a); qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o);
}

/* 1 - |q . r| for unit quaternions: 0 for one rotation */
static double ac_quat_error(qaws_scalar const* q, double const* r)
{
	return 1 - fabs((double)q[0] * r[0] + (double)q[1] * r[1] + (double)q[2] * r[2] + (double)q[3] * r[3]);
}

static void test_vector_quaternion(void)
{
	static unsigned int const t[6] = { 0, 1, 2, 0, 2, 3 };
	qaws_scalar x[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };
	qaws_scalar vec[12] = { 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0 }, quat[16] = { 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0 };
	double const h = sqrt(0.5), eighth[4] = { cos(3.14159265358979323846 / 8), 0, 0, sin(3.14159265358979323846 / 8) };
	double const quarter[4] = { h, 0, 0, h }, ident[4] = { 1, 0, 0, 0 };
	double vmid = 1, qmid = 1, keep = 0, vend = 0, qend = 0, vstart = 0, qstart = 0;
	unsigned int v = 0, iv, iq, nend = 0, nstart = 0, i, j;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d *m = NULL, *a = ac_box(0, 1, 0, 1, 0, 1), *b = ac_box(0.5, 1.5, 0.5, 1.5, 0.5, 1.5), *o = NULL, *box = ac_box(1, 2, -0.1, 0.1, 0, 0.6), *s = NULL;
	qaws_motion_pose p0, p1;
	qaws_motion* mo = NULL;
	qaws_mesh_element e;
	unsigned int n;
	char msg[500];
	/* an edge split: vertex 1 holds (0, 2, 0) and minus a quarter turn about z */
	quat[4] = (qaws_scalar)-h; quat[7] = (qaws_scalar)-h;
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 4; d.face_vertices = t; d.face_count = 2;
	qaws_mesh_3d_create(&m, &d);
	vec[0] = 2; vec[3] = 0; vec[4] = 2;
	iv = ac_add(m, "v", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_VECTOR, vec);
	iq = ac_add(m, "q", QAWS_MESH_VERTEX, 4, QAWS_ATTRIBUTE_QUATERNION, quat);
	if (qaws_mesh_3d_split_edge(&v, m, NULL, NULL, 0, 1) == QAWS_STATUS_OK)
	{
		qaws_scalar* av = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * qaws_mesh_3d_get_vertex_count(m));
		qaws_scalar* aq = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 4 * qaws_mesh_3d_get_vertex_count(m));
		qaws_mesh_3d_get_attribute(av, m, iv);
		qaws_mesh_3d_get_attribute(aq, m, iq);
		vmid = fabs((double)av[3 * v] - 1) + fabs((double)av[3 * v + 1] - 1) + fabs((double)av[3 * v + 2]);
		qmid = ac_quat_error(aq + 4 * v, eighth);
		free(av); free(aq);
	}
	/* A - B, vectors (1, 0, 0) on A, (0, 1, 0) on B */
	{
		qaws_scalar va[24], vb[108];
		qaws_scalar *ov, *src;
		for (i = 0; i < 8; i++)
		{
			va[3 * i] = 1; va[3 * i + 1] = 0; va[3 * i + 2] = 0;
		}
		for (i = 0; i < 36; i++)
		{
			vb[3 * i] = 0; vb[3 * i + 1] = 1; vb[3 * i + 2] = 0;
		}
		ac_add(a, "v", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_VECTOR, va);
		ac_add(b, "v", QAWS_MESH_CORNER, 3, QAWS_ATTRIBUTE_VECTOR, vb);   /* (per corner, A's per vertex: they match) */
		keep = 1;
		if (qaws_mesh_3d_boolean(&o, NULL, a, b, QAWS_BOOLEAN_DIFFERENCE, NULL) == QAWS_STATUS_OK && o && (ov = ac_get(o, "v", &e, &n)) != NULL)
		{
			src = ac_get(o, "source", &e, &n);
			keep = 0;
			for (i = 0; src && i < qaws_mesh_3d_get_face_count(o); i++)
				for (j = 0; j < 3; j++)
				{
					qaws_scalar const* q = ov + 9 * i + 3 * j;
					int in_a = (unsigned int)src[i] < 12;
					keep = fmax(keep, fabs((double)q[0] - in_a) + fabs((double)q[1] - !in_a) + fabs((double)q[2]));
				}
			free(ov); free(src);
		}
	}
	/* a box swept through a quarter turn about z */
	{
		qaws_scalar bv[24], bq[32];
		for (i = 0; i < 8; i++)
		{
			bv[3 * i] = 1; bv[3 * i + 1] = 0; bv[3 * i + 2] = 0;
			bq[4 * i] = 1; bq[4 * i + 1] = 0; bq[4 * i + 2] = 0; bq[4 * i + 3] = 0;
		}
		ac_add(box, "v", QAWS_MESH_VERTEX, 3, QAWS_ATTRIBUTE_VECTOR, bv);
		ac_add(box, "q", QAWS_MESH_VERTEX, 4, QAWS_ATTRIBUTE_QUATERNION, bq);
		memset(&p0, 0, sizeof p0);
		memset(&p1, 0, sizeof p1);
		p0.rotation[0] = 1;
		p1.rotation[0] = (qaws_scalar)h; p1.rotation[3] = (qaws_scalar)h;
		p0.angular_velocity.z = p1.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
		qaws_motion_create_poses(&mo, &p0, &p1, 0);
		if (mo && qaws_mesh_3d_sweep(&s, NULL, box, mo, NULL) == QAWS_STATUS_OK && s)
		{
			qaws_scalar *sv = ac_get(s, "v", &e, &n), *sq = ac_get(s, "q", &e, &n);
			for (i = 0; sv && sq && i < qaws_mesh_3d_get_face_count(s); i++)
			{
				unsigned int tv[3];
				qaws_mesh_3d_get_triangle(tv, s, i);
				for (j = 0; j < 3; j++)
				{
					qaws_vec3 p = qaws_mesh_3d_get_vertex(s, tv[j]);
					double px = (double)p.x, py = (double)p.y, pz = (double)p.z;
					int at_z = fabs(pz) < 1e-12 || fabs(pz - 0.6) < 1e-12;
					qaws_scalar const *cv = sv + 9 * i + 3 * j, *cq = sq + 12 * i + 4 * j;
					/* the end pose's outer corners (x in {-0.1, 0.1}, y = 2) and the start's (x = 2): reached at those poses alone (the inner ones are
					   crossed mid-motion by the other inner edge) */
					if (at_z && (fabs(fabs(px) - 0.1) < 1e-9) && fabs(py - 2) < 1e-9)
					{
						vend = fmax(vend, fabs((double)cv[0]) + fabs((double)cv[1] - 1) + fabs((double)cv[2]));
						qend = fmax(qend, ac_quat_error(cq, quarter));
						nend++;
					}
					if (at_z && (fabs(fabs(py) - 0.1) < 1e-9) && fabs(px - 2) < 1e-9)
					{
						vstart = fmax(vstart, fabs((double)cv[0] - 1) + fabs((double)cv[1]) + fabs((double)cv[2]));
						qstart = fmax(qstart, ac_quat_error(cq, ident));
						nstart++;
					}
				}
			}
			free(sv); free(sq);
		}
	}
	sprintf(msg, "a split's VECTOR midpoint off (1, 1, 0) by %.1e, QUATERNION off the eighth turn by %.1e; A - B: vectors off their operand's by %.1e (B's turned pieces kept); a quarter-turn sweep: %u end corners (vector %.1e, quaternion %.1e), %u start corners (%.1e, %.1e)",
		vmid, qmid, keep, nend, vend, qend, nstart, vstart, qstart);
	printf("  %s\n", msg);
	TEST_ASSERT(vmid < 1e-12 && qmid < 1e-12 && keep < 1e-12 && nend > 0 && nstart > 0 && vend < 1e-9 && qend < 1e-9 && vstart < 1e-9 && qstart < 1e-9, msg);
	qaws_mesh_3d_destroy(m); qaws_mesh_3d_destroy(a); qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(o); qaws_mesh_3d_destroy(box); qaws_mesh_3d_destroy(s);
	qaws_motion_destroy(mo);
}

/* the vertex of m nearest p */
static unsigned int ac_nearest(qaws_mesh_3d const* m, double x, double y)
{
	unsigned int n = qaws_mesh_3d_get_vertex_count(m), i, best = 0;
	double d = 1e300;
	for (i = 0; i < n; i++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, i);
		double e = ((double)p.x - x) * ((double)p.x - x) + ((double)p.y - y) * ((double)p.y - y);
		if (e < d)
		{
			d = e;
			best = i;
		}
	}
	return best;
}

/* the count of m's edges whose edge attribute "w" is v */
static unsigned int ac_edges_valued(qaws_mesh_3d const* m, double v)
{
	unsigned int id, ne = qaws_mesh_3d_get_edge_count(m), i, c = 0;
	qaws_scalar* w = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (ne + 1));
	if (qaws_mesh_3d_find_attribute(&id, m, "w") == QAWS_STATUS_OK)
	{
		qaws_mesh_3d_get_attribute(w, m, id);
		for (i = 0; i < ne; i++) c += (double)w[i] == v;
	}
	free(w);
	return c;
}

static void test_subdivision_rules(void)
{
	static unsigned int const f4[4] = { 0, 1, 2, 3 }, s4[1] = { 4 }, t6[6] = { 0, 1, 2, 0, 2, 3 };
	qaws_scalar x[12] = { 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0 };
	double const c85 = cos(85 * 3.14159265358979323846 / 180), s85 = sin(85 * 3.14159265358979323846 / 180);
	qaws_scalar ang[4] = { 3, -3, -3, 3 }, q[16];
	qaws_scalar wv[5];
	qaws_mesh_3d_desc d;
	qaws_mesh_3d *quad = NULL, *tri = NULL, *o = NULL, *ot = NULL;
	qaws_mesh_subdivide_desc sd;
	double ea = 1, eq = 1, unit = 0;
	unsigned int cnt[6] = { 0, 0, 0, 0, 0, 0 }, tcnt[6] = { 0, 0, 0, 0, 0, 0 }, ne = 0, nte = 0, i, j;
	char msg[500];
	unsigned int e[10];
	/* the left side 170 degrees about z, the right side -170 */
	for (i = 0; i < 4; i++)
	{
		int left = i == 0 || i == 3;
		q[4 * i] = (qaws_scalar)c85; q[4 * i + 1] = q[4 * i + 2] = 0; q[4 * i + 3] = (qaws_scalar)(left ? s85 : -s85);
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 4; d.face_vertices = f4; d.face_sizes = s4; d.face_count = 1;
	qaws_mesh_3d_create(&quad, &d);
	ac_add(quad, "a", QAWS_MESH_VERTEX, 1, QAWS_ATTRIBUTE_ANGLE, ang);
	ac_add(quad, "q", QAWS_MESH_VERTEX, 4, QAWS_ATTRIBUTE_QUATERNION, q);
	ne = qaws_mesh_3d_get_edge_count(quad);
	for (i = 0; i < ne; i++) wv[i] = (qaws_scalar)(1 + i);
	ac_add(quad, "w", QAWS_MESH_EDGE, 1, QAWS_ATTRIBUTE_LINEAR, wv);
	memset(&sd, 0, sizeof sd);
	sd.scheme = QAWS_MESH_SUBDIVIDE_CATMULL_CLARK;
	sd.levels = 2;
	if (qaws_mesh_3d_subdivide(&o, quad, &sd) == QAWS_STATUS_OK && o)
	{
		qaws_mesh_element el;
		unsigned int n, c = ac_nearest(o, 0.5, 0.5), nv = qaws_mesh_3d_get_vertex_count(o);
		qaws_scalar *va = ac_get(o, "a", &el, &n), *vq = ac_get(o, "q", &el, &n);
		double half[4] = { 0, 0, 0, 1 };
		if (va && vq)
		{
			ea = fabs(fabs((double)va[c]) - 3.14159265358979323846);
			eq = ac_quat_error(vq + 4 * c, half);
			for (i = 0; i < nv; i++)
				unit = fmax(unit, fabs(sqrt((double)vq[4 * i] * vq[4 * i] + (double)vq[4 * i + 1] * vq[4 * i + 1] + (double)vq[4 * i + 2] * vq[4 * i + 2] + (double)vq[4 * i + 3] * vq[4 * i + 3]) - 1));
		}
		for (j = 0; j < 6; j++) cnt[j] = ac_edges_valued(o, j);
		free(va); free(vq);
	}
	/* two triangles, Loop once: each of the 5 edges in 2 pieces */
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 4; d.face_vertices = t6; d.face_count = 2;
	qaws_mesh_3d_create(&tri, &d);
	nte = qaws_mesh_3d_get_edge_count(tri);
	qaws_mesh_3d_get_edges(e, tri);
	for (i = 0; i < nte; i++) wv[i] = (qaws_scalar)(1 + i);
	ac_add(tri, "w", QAWS_MESH_EDGE, 1, QAWS_ATTRIBUTE_LINEAR, wv);
	sd.scheme = QAWS_MESH_SUBDIVIDE_LOOP;
	sd.levels = 1;
	if (qaws_mesh_3d_subdivide(&ot, tri, &sd) == QAWS_STATUS_OK && ot)
		for (j = 0; j < 6; j++) tcnt[j] = ac_edges_valued(ot, j);
	sprintf(msg, "subdivided: the ANGLE at the centre off +-pi by %.1e, the QUATERNION off the half turn by %.1e (unit within %.1e); Catmull-Clark x2 edge values 1-4 on %u %u %u %u edges (4 each), 0 on %u (inner); Loop x1 1-5 on %u %u %u %u %u (2 each), 0 on %u (inner: 3 per triangle)",
		ea, eq, unit, cnt[1], cnt[2], cnt[3], cnt[4], cnt[0], tcnt[1], tcnt[2], tcnt[3], tcnt[4], tcnt[5], tcnt[0]);
	printf("  %s\n", msg);
	TEST_ASSERT(ea < 1e-12 && eq < 1e-12 && unit < 1e-12 && cnt[1] == 4 && cnt[2] == 4 && cnt[3] == 4 && cnt[4] == 4 && ne == 4 && nte == 5 &&
		tcnt[1] == 2 && tcnt[2] == 2 && tcnt[3] == 2 && tcnt[4] == 2 && tcnt[5] == 2 && tcnt[0] == 6, msg);
	qaws_mesh_3d_destroy(quad); qaws_mesh_3d_destroy(tri); qaws_mesh_3d_destroy(o); qaws_mesh_3d_destroy(ot);
}

/* m's feature length, crease length (edges whose faces' normals turn by more than turn: 1 - cos) and the edges a feature
   without a crease or the reverse */
static void ac_crease_check(qaws_mesh_3d const* m, double turn_min, double* feat, double* crease, unsigned int* mism)
{
	unsigned int ne = qaws_mesh_3d_get_edge_count(m), nf = qaws_mesh_3d_get_triangle_count(m), nfeat = qaws_mesh_3d_get_feature_edge_count(m), k, j, f;
	unsigned int *e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne), *fe = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (nfeat + 1));
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * nf);
	double* fn = (double*)malloc(sizeof(double) * 3 * nf);
	*feat = *crease = 0;
	*mism = 0;
	qaws_mesh_3d_get_edges(e, m);
	qaws_mesh_3d_get_feature_edges(fe, m);
	qaws_mesh_3d_get_triangles(t, m);
	for (f = 0; f < nf; f++) ac_tri_normal(m, f, fn + 3 * f);
	for (k = 0; k < ne; k++)
	{
		qaws_vec3 p = qaws_mesh_3d_get_vertex(m, e[2 * k]), q = qaws_mesh_3d_get_vertex(m, e[2 * k + 1]);
		double len = sqrt(pow((double)p.x - q.x, 2) + pow((double)p.y - q.y, 2) + pow((double)p.z - q.z, 2)), turn = 0;
		int is_feat = 0, side = 0;
		unsigned int g[2] = { 0, 0 };
		for (j = 0; j < nfeat; j++)
			if ((fe[2 * j] == e[2 * k] && fe[2 * j + 1] == e[2 * k + 1]) || (fe[2 * j] == e[2 * k + 1] && fe[2 * j + 1] == e[2 * k]))
				is_feat = 1;
		for (f = 0; f < nf && side < 2; f++)
			for (j = 0; j < 3; j++)
				if ((t[3 * f + j] == e[2 * k] && t[3 * f + (j + 1) % 3] == e[2 * k + 1]) || (t[3 * f + j] == e[2 * k + 1] && t[3 * f + (j + 1) % 3] == e[2 * k]))
					g[side++] = f;
		if (side == 2)
			turn = 1 - (fn[3 * g[0]] * fn[3 * g[1]] + fn[3 * g[0] + 1] * fn[3 * g[1] + 1] + fn[3 * g[0] + 2] * fn[3 * g[1] + 2]);
		*feat += is_feat * len;
		*crease += (turn > turn_min) * len;
		*mism += (unsigned int)((turn > turn_min) != is_feat);
	}
	free(e); free(fe); free(t); free(fn);
}

/* a closed 16-sided prism of radius 0.5 round the z axis, z in [0, 1], its caps fanned */
static qaws_mesh_3d* ac_prism(void)
{
	qaws_scalar x[96];
	unsigned int t[3 * 60], i, k = 0;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < 32; i++)
	{
		double a = 2 * 3.14159265358979323846 * (i % 16) / 16;
		x[3 * i] = (qaws_scalar)(0.5 * cos(a)); x[3 * i + 1] = (qaws_scalar)(0.5 * sin(a)); x[3 * i + 2] = (qaws_scalar)(i < 16 ? 0 : 1);
	}
	for (i = 0; i < 16; i++)
	{
		unsigned int a = i, b = (i + 1) % 16;
		t[k++] = a; t[k++] = b; t[k++] = b + 16;
		t[k++] = a; t[k++] = b + 16; t[k++] = a + 16;
	}
	for (i = 1; i + 1 < 16; i++)
	{
		t[k++] = 0; t[k++] = i + 1; t[k++] = i;
		t[k++] = 16; t[k++] = 16 + i; t[k++] = 17 + i;
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 32; d.face_vertices = t; d.face_count = k / 3;
	qaws_mesh_3d_create(&m, &d);
	return m;
}

static void test_sweep_features(void)
{
	qaws_mesh_3d *c = ac_box(0, 1, 0, 1, 0, 1), *o = NULL;
	qaws_motion* m = ac_move(0, 0, 0, 1, 0.5, 0);
	double feat = 0, crease = 0;
	unsigned int mism = 1000;
	qaws_status st;
	char msg[300];
	qaws_mesh_3d_detect_features(c, (qaws_scalar)0.5, 0, 0);
	st = qaws_mesh_3d_sweep(&o, NULL, c, m, NULL);
	if (st == QAWS_STATUS_OK && o)
		ac_crease_check(o, 1e-9, &feat, &crease, &mism);
	sprintf(msg, "a featured cube swept along (1, 0.5, 0): features along %.12f, creases along %.12f, %u edges a feature without a crease or the reverse",
		feat, crease, mism);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && mism == 0 && fabs(feat - crease) < 1e-12 && feat > 12, msg);
	qaws_mesh_3d_destroy(o);
	qaws_motion_destroy(m);
	/* along (1, 0, 0): its x edges slide along themselves */
	o = NULL;
	mism = 1000;
	feat = crease = 0;
	m = ac_move(0, 0, 0, 1, 0, 0);
	st = qaws_mesh_3d_sweep(&o, NULL, c, m, NULL);
	if (st == QAWS_STATUS_OK && o)
		ac_crease_check(o, 1e-9, &feat, &crease, &mism);
	sprintf(msg, "a featured cube swept along (1, 0, 0): features along %.12f, creases along %.12f, %u edges a feature without a crease or the reverse",
		feat, crease, mism);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && o && mism == 0 && fabs(feat - 16) < 1e-12 && fabs(crease - 16) < 1e-12, msg);
	qaws_mesh_3d_destroy(c); qaws_mesh_3d_destroy(o);
	qaws_motion_destroy(m);
	{
		/* a prism's rims swept sideways: the creases where the rims run along the motion (their envelope) traced by
		   vertices that are no feature corners */
		qaws_mesh_3d* pr = ac_prism();
		o = NULL;
		mism = 1000;
		feat = crease = 0;
		m = ac_move(0, 0, 0, 1, 0, 0);
		qaws_mesh_3d_detect_features(pr, (qaws_scalar)0.5, 0, 0);
		st = qaws_mesh_3d_sweep(&o, NULL, pr, m, NULL);
		if (st == QAWS_STATUS_OK && o)
			ac_crease_check(o, 1 - cos(0.5), &feat, &crease, &mism);
		sprintf(msg, "a featured 16-sided prism swept along (1, 0, 0): features along %.12f, creases along %.12f, %u edges a feature without a crease or the reverse",
			feat, crease, mism);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && o && mism == 0 && fabs(feat - crease) < 1e-12 && crease > 4, msg);
		qaws_mesh_3d_destroy(pr); qaws_mesh_3d_destroy(o);
		qaws_motion_destroy(m);
	}
	{
		/* along a curve, turning a quarter */
		qaws_mesh_3d* cube = ac_box(-0.3, 0.3, -0.3, 0.3, -0.3, 0.3);
		qaws_motion_pose s, e;
		qaws_motion* mt = NULL;
		qaws_sweep_desc sd;
		qaws_sweep_report r;
		qaws_mesh_topology tp;
		unsigned int uncovered;
		o = NULL;
		memset(&s, 0, sizeof s);
		memset(&e, 0, sizeof e);
		memset(&tp, 0, sizeof tp);
		s.rotation[0] = 1;
		e.rotation[0] = (qaws_scalar)sqrt(0.5); e.rotation[3] = (qaws_scalar)sqrt(0.5);
		s.position.x = -0.5; s.position.y = 1; s.position.z = 1.05;
		e.position.x = 4.5; e.position.y = 1; e.position.z = 1.05;
		s.velocity.x = e.velocity.x = 5; s.velocity.y = (qaws_scalar)2.4; e.velocity.y = (qaws_scalar)-2.4;
		s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
		qaws_motion_create_poses(&mt, &s, &e, 0);
		qaws_mesh_3d_detect_features(cube, (qaws_scalar)0.5, 0, 0);
		memset(&sd, 0, sizeof sd);
		sd.tolerance = (qaws_scalar)0.004;
		mism = 1000;
		st = qaws_mesh_3d_sweep(&o, &r, cube, mt, &sd);
		if (st == QAWS_STATUS_OK && o)
		{
			qaws_mesh_3d_get_topology(&tp, o);
			ac_crease_check(o, 1 - cos(25 * 3.14159265358979323846 / 180), &feat, &crease, &mism);
		}
		{
			/* the cube at 101 times: points 0.01 within its four sides, at mid height and 0.01 below its top */
			qaws_mesh_bvh* bv = NULL;
			qaws_scalar pts[3 * 160], wn[160];
			unsigned int ti, pi2;
			uncovered = 1000;
			if (o && qaws_mesh_bvh_create_mesh_3d(&bv, o, NULL) == QAWS_STATUS_OK)
			{
				uncovered = 0;
				for (ti = 0; ti <= 100; ti++)
				{
					qaws_motion_state ms;
					qaws_motion_evaluate(&ms, mt, (qaws_scalar)(ti / 100.0));
					for (pi2 = 0; pi2 < 160; pi2++)
					{
						unsigned int sd2 = pi2 / 40 % 4;
						double a = -0.29 + 0.58 * (pi2 % 20) / 19.0, lx = sd2 == 0 ? 0.29 : sd2 == 1 ? -0.29 : a, ly = sd2 == 2 ? 0.29 : sd2 == 3 ? -0.29 : a, lz = pi2 % 40 < 20 ? 0 : 0.29;
						pts[3 * pi2] = (qaws_scalar)(ms.rotation[0] * lx + ms.rotation[1] * ly + ms.rotation[2] * lz + ms.position.x);
						pts[3 * pi2 + 1] = (qaws_scalar)(ms.rotation[3] * lx + ms.rotation[4] * ly + ms.rotation[5] * lz + ms.position.y);
						pts[3 * pi2 + 2] = (qaws_scalar)(ms.rotation[6] * lx + ms.rotation[7] * ly + ms.rotation[8] * lz + ms.position.z);
					}
					qaws_mesh_bvh_winding_number(wn, bv, pts, 160, NULL);
					for (pi2 = 0; pi2 < 160; pi2++)
						uncovered += wn[pi2] < 0.5;
				}
			}
			qaws_mesh_bvh_destroy(bv);
		}
		sprintf(msg, "a featured cube along a curve, turning a quarter (%u slabs): closed %d, %u piece(s); features along %.6f, edges turning past 25 degrees along %.6f, %u edges off; %u of the cube's points at 101 times outside",
			r.slabs, tp.closed, tp.components, feat, crease, mism, uncovered);
		printf("  %s\n", msg);
		TEST_ASSERT(st == QAWS_STATUS_OK && o && tp.closed && tp.components == 1 && mism == 0 && uncovered == 0, msg);
		qaws_mesh_3d_destroy(cube); qaws_mesh_3d_destroy(o);
		qaws_motion_destroy(mt);
	}
}

/* an open tube of 8 sides and one row about z, from z0 to z1, its edges valued base + index */
static qaws_mesh_3d* ac_tube(double z0, double z1, double base)
{
	qaws_scalar x[48], w[64];
	unsigned int t[48], i, k = 0, ne;
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	for (i = 0; i < 16; i++)
	{
		double a = 2 * 3.14159265358979323846 * (i % 8) / 8;
		x[3 * i] = (qaws_scalar)(0.5 * cos(a)); x[3 * i + 1] = (qaws_scalar)(0.5 * sin(a)); x[3 * i + 2] = (qaws_scalar)(i < 8 ? z0 : z1);
	}
	for (i = 0; i < 8; i++)
	{
		unsigned int a = i, b = (i + 1) % 8;
		t[k++] = a; t[k++] = b; t[k++] = b + 8;
		t[k++] = a; t[k++] = b + 8; t[k++] = a + 8;
	}
	memset(&d, 0, sizeof d);
	d.positions = x; d.vertex_count = 16; d.face_vertices = t; d.face_count = 16;
	qaws_mesh_3d_create(&m, &d);
	ne = qaws_mesh_3d_get_edge_count(m);
	for (i = 0; i < ne; i++) w[i] = (qaws_scalar)(base + i);
	ac_add(m, "w", QAWS_MESH_EDGE, 1, QAWS_ATTRIBUTE_LINEAR, w);
	return m;
}

/* the largest difference between t's corner values of the corner attribute name (LINEAR) and r's at the same places (a point a
   millionth of the way into t's triangle, on r's triangle of the same plane holding it, interpolated); *unmatched: places on
   no such triangle */
static double ac_corner_match(qaws_mesh_3d const* r, qaws_mesh_3d const* t, char const* name, unsigned int* unmatched)
{
	qaws_mesh_element er, et;
	unsigned int nr = 0, nt = 0, ntr = qaws_mesh_3d_get_triangle_count(r), ntt = qaws_mesh_3d_get_triangle_count(t), f, g, j, k;
	qaws_scalar *vr = ac_get(r, name, &er, &nr), *vt = ac_get(t, name, &et, &nt);
	double worst = 0;
	*unmatched = 0;
	for (f = 0; vr && vt && f < ntt; f++)
	{
		unsigned int c[3];
		double P[3][3], n[3];
		qaws_mesh_3d_get_triangle(c, t, f);
		ac_tri_normal(t, f, n);
		for (j = 0; j < 3; j++)
		{
			qaws_vec3 x = qaws_mesh_3d_get_vertex(t, c[j]);
			P[j][0] = x.x; P[j][1] = x.y; P[j][2] = x.z;
		}
		for (j = 0; j < 3; j++)
		{
			double x[3], best = 1e30, val[4] = { 0, 0, 0, 0 };
			for (k = 0; k < 3; k++)
				x[k] = P[j][k] + 1e-6 * ((P[0][k] + P[1][k] + P[2][k]) / 3 - P[j][k]);
			for (g = 0; g < ntr; g++)
			{
				unsigned int d[3], i;
				double Q[3][3], m[3], u[3], v[3], q[3], d00 = 0, d01 = 0, d11 = 0, d20 = 0, d21 = 0, den, w1, w2, hgt = 0, out;
				ac_tri_normal(r, g, m);
				if (m[0] * n[0] + m[1] * n[1] + m[2] * n[2] < 0.999)
					continue;
				qaws_mesh_3d_get_triangle(d, r, g);
				for (i = 0; i < 3; i++)
				{
					qaws_vec3 y = qaws_mesh_3d_get_vertex(r, d[i]);
					Q[i][0] = y.x; Q[i][1] = y.y; Q[i][2] = y.z;
				}
				for (k = 0; k < 3; k++)
				{
					u[k] = Q[1][k] - Q[0][k];
					v[k] = Q[2][k] - Q[0][k];
					q[k] = x[k] - Q[0][k];
					hgt += q[k] * m[k];
					d00 += u[k] * u[k]; d01 += u[k] * v[k]; d11 += v[k] * v[k]; d20 += q[k] * u[k]; d21 += q[k] * v[k];
				}
				den = d00 * d11 - d01 * d01;
				w1 = (d11 * d20 - d01 * d21) / den;
				w2 = (d00 * d21 - d01 * d20) / den;
				out = fabs(hgt) + fmax(0, -fmin(fmin(w1, w2), 1 - w1 - w2));
				if (out < best)
				{
					best = out;
					for (k = 0; k < nr && k < 4; k++)
						val[k] = (1 - w1 - w2) * (double)vr[(3 * g) * nr + k] + w1 * (double)vr[(3 * g + 1) * nr + k] + w2 * (double)vr[(3 * g + 2) * nr + k];
				}
			}
			if (best > 1e-9)
			{
				(*unmatched)++;
				continue;
			}
			for (k = 0; k < nt && k < 4; k++)
				worst = fmax(worst, fabs(val[k] - (double)vt[(3 * f + j) * nt + k]));
		}
	}
	free(vr); free(vt);
	return worst;
}

/* the length of m's edges whose value (edge attribute name) rounds to each of 1 .. 5, into l[0 .. 4] */
static void ac_edge_lengths(qaws_mesh_3d const* m, char const* name, double* l)
{
	unsigned int id, ne = qaws_mesh_3d_get_edge_count(m), i, *ev = (unsigned int*)malloc(sizeof(unsigned int) * (2 * ne + 1));
	qaws_scalar* v = (qaws_scalar*)malloc(sizeof(qaws_scalar) * (ne + 1));
	for (i = 0; i < 5; i++)
		l[i] = 0;
	if (qaws_mesh_3d_find_attribute(&id, m, name) == QAWS_STATUS_OK && ev && v)
	{
		qaws_mesh_3d_get_edges(ev, m);
		qaws_mesh_3d_get_attribute(v, m, id);
		for (i = 0; i < ne; i++)
		{
			long k = lround((double)v[i]);
			qaws_vec3 a = qaws_mesh_3d_get_vertex(m, ev[2 * i]), b = qaws_mesh_3d_get_vertex(m, ev[2 * i + 1]);
			if (k >= 1 && k <= 5)
				l[k - 1] += sqrt(((double)a.x - b.x) * ((double)a.x - b.x) + ((double)a.y - b.y) * ((double)a.y - b.y) + ((double)a.z - b.z) * ((double)a.z - b.z));
		}
	}
	free(ev); free(v);
}

static void test_sweep_simplify(void)
{
	qaws_mesh_3d *c = ac_box(0, 1, 0, 1, 0, 0.6), *o = NULL, *s = NULL;
	qaws_scalar uv[72], *ev;
	unsigned int f[36], i, k, ne, fo = 0, fs = 0;
	qaws_motion_pose p0, p1;
	qaws_motion* mo = NULL;
	qaws_sweep_desc d;
	qaws_sweep_report ro, rs;
	double const h = sqrt(0.5);
	double lo[5], ls[5], err = 0, vals = 1;
	unsigned int um = 1, um2 = 1;
	char msg[400];
	qaws_mesh_3d_detect_features(c, (qaws_scalar)0.5, 0, 0);
	ac_vertex_field(c, "lin", 1, 2, 3, 0);
	qaws_mesh_3d_get_triangles(f, c);
	for (i = 0; i < 12; i++)
		for (k = 0; k < 3; k++)
		{
			qaws_vec3 x = qaws_mesh_3d_get_vertex(c, f[3 * i + k]);
			uv[6 * i + 2 * k] = (qaws_scalar)(x.x + x.z + (double)(i / 2));   /* (each side its own chart: seams at the creases) */
			uv[6 * i + 2 * k + 1] = (qaws_scalar)(x.y - x.z);
		}
	ac_add(c, "uv", QAWS_MESH_CORNER, 2, QAWS_ATTRIBUTE_LINEAR, uv);
	ne = qaws_mesh_3d_get_edge_count(c);
	ev = (qaws_scalar*)malloc(sizeof(qaws_scalar) * ne);
	for (i = 0; i < ne; i++)
		ev[i] = (qaws_scalar)(1 + i % 5);
	ac_add(c, "e", QAWS_MESH_EDGE, 1, QAWS_ATTRIBUTE_LINEAR, ev);
	free(ev);
	memset(&p0, 0, sizeof p0);
	memset(&p1, 0, sizeof p1);
	p0.rotation[0] = 1;
	p1.rotation[0] = (qaws_scalar)h; p1.rotation[3] = (qaws_scalar)h;
	p0.angular_velocity.z = p1.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
	p1.position.x = 3; p1.position.y = 1;
	p0.velocity.x = 4; p1.velocity.x = 2; p1.velocity.y = 2;
	qaws_motion_create_poses(&mo, &p0, &p1, 0);
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)0.004;
	memset(&ro, 0, sizeof ro);
	memset(&rs, 0, sizeof rs);
	if (mo && qaws_mesh_3d_sweep(&o, &ro, c, mo, &d) == QAWS_STATUS_OK && o)
	{
		d.simplify = 1;
		if (qaws_mesh_3d_sweep(&s, &rs, c, mo, &d) == QAWS_STATUS_OK && s)
		{
			vals = fmax(ac_corner_match(o, s, "lin", &um), ac_corner_match(o, s, "uv", &um2));
			ac_edge_lengths(o, "e", lo);
			ac_edge_lengths(s, "e", ls);
			for (i = 0; i < 5; i++)
				err = fmax(err, fabs(lo[i] - ls[i]) / fmax(lo[i], 1e-300));
			fo = qaws_mesh_3d_get_face_count(o);
			fs = qaws_mesh_3d_get_face_count(s);
		}
	}
	sprintf(msg, "a sweep simplified: %u faces from %u, its corners' values off the plain sweep's there by %.1e (%u places off it), edge values' lengths by %.1e (relative), features %u (%u), closed %d",
		fs, fo, vals, um + um2, err, s ? qaws_mesh_3d_get_feature_edge_count(s) : 0, o ? qaws_mesh_3d_get_feature_edge_count(o) : 0, rs.closed);
	printf("  %s\n", msg);
	TEST_ASSERT(o && s && fs > 0 && 2 * fs < fo && vals < 1e-5 && um + um2 == 0 && err < 1e-9 && rs.closed && qaws_mesh_3d_get_feature_edge_count(s) == qaws_mesh_3d_get_feature_edge_count(o), msg);
	qaws_mesh_3d_destroy(c); qaws_mesh_3d_destroy(o); qaws_mesh_3d_destroy(s);
	qaws_motion_destroy(mo);
}

static void test_connect_edges(void)
{
	qaws_mesh_3d *a = ac_tube(0, 1, 1), *b = ac_tube(1.6, 2.6, 100), *m = NULL;
	qaws_connect_desc cd;
	unsigned int ia = 0, ib = 0, im = 0, ne = 0, bad = 0, kept = 0, bridge = 0, i;
	qaws_status st;
	char msg[300];
	memset(&cd, 0, sizeof cd);
	cd.quads = 1;
	qaws_mesh_3d_find_attribute(&ia, a, "w");
	qaws_mesh_3d_find_attribute(&ib, b, "w");
	st = qaws_mesh_3d_connect(&m, NULL, a, 8, b, 0, &cd);
	if (st == QAWS_STATUS_OK && m && qaws_mesh_3d_find_attribute(&im, m, "w") == QAWS_STATUS_OK)
	{
		unsigned int* e;
		qaws_scalar* w;
		ne = qaws_mesh_3d_get_edge_count(m);
		e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * ne);
		w = (qaws_scalar*)malloc(sizeof(qaws_scalar) * ne);
		qaws_mesh_3d_get_edges(e, m);
		qaws_mesh_3d_get_attribute(w, m, im);
		for (i = 0; i < ne; i++)
		{
			unsigned int x = e[2 * i], y = e[2 * i + 1];
			double want = x < 16 && y < 16 ? ac_edge_value(a, ia, x, y) : x >= 16 && y >= 16 && x < 32 && y < 32 ? ac_edge_value(b, ib, x - 16, y - 16) : 0;
			if (want < 0)
				want = 0;   /* (two of a's vertices on no edge of a: a bridge edge) */
			bad += (double)w[i] != want;
			kept += want != 0;
			bridge += want == 0;
		}
		free(e); free(w);
	}
	sprintf(msg, "two tubes joined: %u edges, %u keeping their tube's value, %u of the bridge at 0, %u off", ne, kept, bridge, bad);
	printf("  %s\n", msg);
	TEST_ASSERT(st == QAWS_STATUS_OK && m && bad == 0 && kept == 64 && bridge > 0, msg);   /* (32 edges per tube) */
	qaws_mesh_3d_destroy(a); qaws_mesh_3d_destroy(b); qaws_mesh_3d_destroy(m);
}

static void test_refusals(void)
{
	qaws_mesh_3d *a = ac_box(0, 1, 0, 1, 0, 1), *b = ac_box(0.5, 1.5, 0.5, 1.5, 0.5, 1.5), *o = NULL;
	qaws_attribute_desc d;
	unsigned int id;
	qaws_status s1, s2;
	char msg[200];
	ac_vertex_field(a, "field", 1, 0, 0, 0);
	ac_face_constant(b, "field", 1);
	s1 = qaws_mesh_3d_boolean(&o, NULL, a, b, QAWS_BOOLEAN_UNION, NULL);
	d.name = "n2";
	d.element = QAWS_MESH_VERTEX;
	d.components = 2;
	d.rule = QAWS_ATTRIBUTE_NORMAL;
	s2 = qaws_mesh_3d_add_attribute(&id, a, NULL, &d);
	sprintf(msg, "refused: \"field\" per vertex on A and per face on B (%d), a NORMAL of 2 components (%d)", (int)s1, (int)s2);
	printf("  %s\n", msg);
	TEST_ASSERT(s1 == QAWS_STATUS_INVALID_ARGUMENT && !o && s2 == QAWS_STATUS_INVALID_ARGUMENT, msg);
	qaws_mesh_3d_destroy(a); qaws_mesh_3d_destroy(b);
}

static void test_sweep_top(void)
{
	qaws_mesh_3d *c = ac_box(-0.2, 0.2, -0.2, 0.2, -0.2, 0.2), *o = NULL;
	qaws_motion* m = NULL;
	qaws_motion_pose s, e;
	qaws_sweep_desc d;
	qaws_attribute_desc ad;
	qaws_scalar side[36], *v = NULL;
	qaws_mesh_element el;
	unsigned int i, k, id, n, nf, nc, *fv = NULL, *fs = NULL, up = 0, top = 0;
	char msg[300];
	for (i = 0; i < 36; i++)
		side[i] = (qaws_scalar)(i / 6);
	memset(&ad, 0, sizeof ad);
	ad.name = "side"; ad.element = QAWS_MESH_CORNER; ad.components = 1; ad.rule = QAWS_ATTRIBUTE_LINEAR;
	qaws_mesh_3d_add_attribute(&id, c, side, &ad);
	memset(&s, 0, sizeof s);
	memset(&e, 0, sizeof e);
	s.rotation[0] = 1;
	e.rotation[0] = (qaws_scalar)cos(3.14159265358979323846 / 4);
	e.rotation[3] = (qaws_scalar)sin(3.14159265358979323846 / 4);
	s.angular_velocity.z = e.angular_velocity.z = (qaws_scalar)(3.14159265358979323846 / 2);
	e.position.x = 2;
	s.velocity.x = 1; s.velocity.y = 2;
	e.velocity.x = 1; e.velocity.y = -2;
	qaws_motion_create_poses(&m, &s, &e, 0);
	memset(&d, 0, sizeof d);
	d.tolerance = (qaws_scalar)0.01;
	qaws_mesh_3d_sweep(&o, NULL, c, m, &d);
	if (o && (v = ac_get(o, "side", &el, &n)) != NULL && el == QAWS_MESH_CORNER)
	{
		nf = qaws_mesh_3d_get_face_count(o);
		nc = qaws_mesh_3d_get_corner_count(o);
		fv = (unsigned int*)malloc(sizeof(unsigned int) * (nc + 1));
		fs = (unsigned int*)malloc(sizeof(unsigned int) * (nf + 1));
		qaws_mesh_3d_get_faces(fv, fs, o);
		for (i = 0, k = 0; i < nf; k += fs[i++])
		{
			qaws_vec3 a = qaws_mesh_3d_get_vertex(o, fv[k]), b = qaws_mesh_3d_get_vertex(o, fv[k + 1]), q = qaws_mesh_3d_get_vertex(o, fv[k + 2]);
			double u[3] = { b.x - a.x, b.y - a.y, b.z - a.z }, w[3] = { q.x - a.x, q.y - a.y, q.z - a.z }, nz, l;
			unsigned int j, all = 1;
			nz = u[0] * w[1] - u[1] * w[0];
			l = sqrt(pow(u[1] * w[2] - u[2] * w[1], 2) + pow(u[2] * w[0] - u[0] * w[2], 2) + nz * nz);
			if (!(l > 0) || nz < 0.999 * l)
				continue;
			up++;
			for (j = 0; j < fs[i]; j++)
				all &= fabs((double)v[k + j] - 1) < 1e-9;
			top += all;
		}
	}
	sprintf(msg, "a face sliding along itself: %u upward faces of the sweep, %u of them the top's number at every corner", up, top);
	printf("  %s\n", msg);
	TEST_ASSERT(o && up > 0 && top == up, msg);
	free(v); free(fv); free(fs);
	qaws_motion_destroy(m);
	qaws_mesh_3d_destroy(c); qaws_mesh_3d_destroy(o);
}

int test_181_mesh_attribute_carry_main(void)
{
	g_pass = 0;
	g_fail = 0;
	printf("Test 181: attributes through Booleans and sweeps\n");
	test_colour_spaces();
	test_boolean();
	test_sweep();
	test_carving();
	test_edge_edits();
	test_edge_boolean();
	test_vector_quaternion();
	test_subdivision_rules();
	test_sweep_features();
	test_sweep_simplify();
	test_sweep_top();
	test_connect_edges();
	test_refusals();
	printf("  Results: %d passed, %d failed\n", g_pass, g_fail);
	return g_fail;
}
