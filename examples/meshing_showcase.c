/*
 * meshing_showcase.c - the meshing module: meshes, references, remeshing,
 * quads, tetrahedra and hexahedra
 *
 * Writes PNG figures into ./showcase/:
 *   mesh1_core.png     the corner store: a torus and its dual (the stored
 *                      adjacency), a cube's creases kept through random edits,
 *                      a linear attribute kept exact through edits
 *   mesh2_reference.png  the surface between the vertices: a coarse blob flat and
 *                      as its PN reference, a coarse can rounded without features
 *                      and sharp with its creases, curvature from the reference
 *   mesh3_remesh.png   remeshing on the reference: uniform and adaptive by
 *                      curvature, creases kept, attributes and their seams transferred
 *   mesh4_quads.png    quads by Instant Meshes on the reference: a cross field and
 *                      its singularities, quad-dominant and pure, along creases, a torus
 *   mesh5_generate.png meshes of exact geometry: a B-spline patch in triangles,
 *                      quads and adaptive, a NURBS torus, a shell of six faces
 *   mesh6_tets.png     tetrahedra filling closed surfaces, cut open: a torus, a box,
 *                      a blob
 *   mesh7_hexes.png    hexahedra on a grid padded to the surface (a box, a torus) and
 *                      the torus's frame field with its singular edges
 *   mesh8_frame_hexes.png  hexahedra following frame fields (integer-grid maps): a turned
 *                      box against the grid, a bent bar, a twisted bar
 *   mesh9..mesh14      quality measures, angle-bounded remeshing, CVT remeshing,
 *                      smoothest cross fields, quads by min-cost flow, quad finishing
 *   mesh16_guidance.png  quads guided: density paint, a stroke of directions, an
 *                      edge loop, pinned vertices; placed cones, an attractor
 *   mesh18_repair.png  repair for farq: a broken torus soup (holes, a gap, faces
 *                      turned, duplicates) made closed and outward, the fills tinted
 *   mesh19_layers.png  boundary layers: prisms grown into a blob, hexahedra into a
 *                      box in quads, the blob's layers with tetrahedra inside
 *   mesh20_curved.png  curved elements: a coarse torus in quads straight and at
 *                      order 2 through its nodes on the torus, a patch in cubic triangles
 *   mesh21_metric.png  metric remeshing: a tanh layer remeshed to its Hessian, a
 *                      cylinder stretched along its axis
 *   mesh22_bvh.png     spatial queries on a torus soup with holes: ray traced
 *                      through the BVH, its winding number and signed distance
 *   mesh23_quality.png quality reports: a 2D mesh by shape, a surface by smallest
 *                      angle, tetrahedra by smallest dihedral, with histograms
 *   mesh24_repair_more.png  repair completed: a box's hole filled along its
 *                      creases (the corner back), crossing faces found exactly
 *   mesh25_watertight.png  watertight surfaces: a holed sphere closed, two
 *                      crossing cubes' union, an open sheet wrapped at an offset
 *   mesh26_intersect.png  self-intersections resolved: two spheres cut along
 *                      their seam (union and lens by picking faces), three cubes
 *   mesh27_boolean.png mesh Booleans: two spheres united, intersected, one taken
 *                      from the other; the CSG emblem from five Booleans
 *   mesh28_geodesic.png  geodesics by the heat method: a torus's distance and
 *                      paths, a bumpy sphere's nearest of three sources, paths
 *                      round a hill
 *   mesh29_smooth.png  smoothing and fairing: a noisy torus by Laplacian, Taubin
 *                      and bi-Laplacian steps, a roof's crease kept as a feature,
 *                      a sphere's cap faired as a membrane and a thin plate
 *   mesh30_deform.png  deformation and parameterization: a bar bent by ARAP and
 *                      by Laplacian editing, a hill flattened with its (u, v) grid
 *   mesh31_deform_2.png  the ARAP paper's tests: a plate with 4 x 4 bumps lifted
 *                      on one side, a square bar twisted half a turn; ARAP against
 *                      Laplacian and biharmonic editing
 *   mesh32_decimate.png  decimation by quadric error: a rippled torus at 20%, 5% and
 *                      1.5% with its certified deviation, a fine cube with its creases
 *   mesh33_decimate_views.png  what decimation keeps: a rippled sphere's outline from one
 *                      view, its top from views above, a can's rims as creases, importance
 *   mesh34_remesh_attributes.png  remeshing from attributes: isolines kept as lines,
 *                      two attributes' Hessian metrics intersected, graded
 *   mesh35_connect.png the mesh connector: two tubes at an angle joined (faired too),
 *                      a funnel between tubes of 32 and 12 sides
 *   mesh36_connect_more.png  the connector finished: a planar strip, an attribute
 *                      carried, two hexahedral blocks joined by hexahedra
 *   mesh37_repulsive_mesh_curves.png  repulsive curves on a mesh: loops across a fold
 *                      by chords and along the surface, a star on a sphere smoothed
 *   mesh38_repulsive_mesh_curves_2.png  repulsive curves on a terrain of Gaussian
 *                      hills: a branching tree spreading, a flower loop opening
 *   mesh39_offset.png  offsets of meshes: contoured, mitred (a slot closed), a shell,
 *                      a sheet thickened
 *   mesh40_subdivide.png  subdivision: a cage, Catmull-Clark, semi-sharp creases, Loop
 *                      with a crease and a colour carried
 *   mesh41_volume_improve.png  tetrahedra improved: slivers smoothed and flipped away,
 *                      the border kept; a blob's cells; dihedral histograms
 *   mesh42_volume_curve.png  curved volume elements: a sphere's tetrahedra and a ball's
 *                      hexahedra, straight and quadratic on a sphere
 *   mesh43_tetrahedra_border.png  tetrahedra keeping the border: Schonhardt's polyhedron,
 *                      its recovery's border points moved inside; a sphere's rings
 *   mesh44_polycube.png  hexahedra of polycube maps: a sphere, a blob, a ring of genus
 *                      one, a bent bar
 *   mesh45_quads_pairing.png  quads from the input's faces: triangles paired, paired and
 *                      split (pure), split
 *   mesh46_quad_layout.png  quad layouts: GLOBAL's quads parted by their motorcycle graph
 *                      into grids (a blob, a can, a torus, a wavy patch)
 *   mesh47_denoise.png  feature-preserving denoising: a noisy box by Taubin and
 *                      bilateral, a noisy blob
 *   mesh48_points.png  watertight surfaces from oriented point clouds: a ring's points,
 *                      noisy points, a cube's points with sharp edges
 *   mesh49_geodesic_exact.png  exact geodesics: a cube from a corner, a hill by the
 *                      heat method and exact
 *   mesh50_butterfly.png  the modified butterfly beside Loop: through the cage's
 *                      vertices, or inside the cage
 *   mesh51_decimate_attributes.png  decimation with the colour in the quadric: a torus
 *                      painted in bands, by the shape alone, with the colour
 *   mesh52_offset_round.png  rounded offsets: the U moved (mitred), rounded out, rounded in
 *   mesh53_volume_cubic.png  cubic volume elements: a ball of hexahedra straight, Q2, Q3
 *   mesh54_subdivide_adaptive.png  adaptive subdivision: a cube by uniform and adaptive Loop
 *
 * Each figure lives in meshing_showcase/NN_*.c, compiled as one translation
 * unit with the rasterizer in example_raster.h.
 */

#include "qaws.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <direct.h>
#define MAKE_DIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MAKE_DIR(p) mkdir(p, 0755)
#endif

#include "example_raster.h"

#define MS_PI 3.14159265358979323846

/* ---- shared helpers ---------------------------------------------------- */

static unsigned int ms_rand(unsigned int* s)
{
	*s = *s * 1664525u + 1013904223u;
	return *s >> 8;
}

static double ms_uniform(unsigned int* s)
{
	return (ms_rand(s) & 0xFFFFFF) / (double)0x1000000;
}

/* blue to red through white */
static void ms_colormap(double t, double* rgb)
{
	static double const k[5][3] = { { 59, 76, 192 }, { 141, 176, 254 }, { 221, 221, 221 }, { 244, 154, 123 }, { 180, 4, 38 } };
	double s;
	int i, j;
	t = t < 0 ? 0 : t > 1 ? 1 : t;
	s = t * 4;
	i = (int)s;
	if (i > 3)
		i = 3;
	s -= i;
	for (j = 0; j < 3; j++)
		rgb[j] = k[i][j] + s * (k[i + 1][j] - k[i][j]);
}

static qaws_mesh_3d* ms_make(double const* x, unsigned int nv, unsigned int const* f, unsigned int const* fs, unsigned int nf)
{
	qaws_mesh_3d_desc d;
	qaws_mesh_3d* m = NULL;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv ? nv : 1));
	unsigned int i;
	for (i = 0; i < 3 * nv; i++)
		p[i] = (qaws_scalar)x[i];
	memset(&d, 0, sizeof d);
	d.positions = p;
	d.vertex_count = nv;
	d.face_vertices = f;
	d.face_sizes = fs;
	d.face_count = nf;
	qaws_mesh_3d_create(&m, &d);
	free(p);
	return m;
}

/* the mesh's positions as doubles (caller frees) */
static double* ms_positions(qaws_mesh_3d const* m)
{
	unsigned int n = qaws_mesh_3d_get_vertex_count(m), i;
	qaws_scalar* p = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (n ? n : 1));
	double* x = (double*)malloc(sizeof(double) * 3 * (n ? n : 1));
	qaws_mesh_3d_get_positions(p, m);
	for (i = 0; i < 3 * n; i++)
		x[i] = (double)p[i];
	free(p);
	return x;
}

/* A mesh drawn: smooth (vertex normals) or flat, a colour per vertex
   (vrgb, 3 each) or one colour, a wireframe of its edges (wire may be NULL). */
static void ms_draw(er_canvas* c, qaws_mesh_3d const* m, double const* rgb, double const* vrgb, int smooth, double const* wire, double wire_width)
{
	unsigned int nv = qaws_mesh_3d_get_vertex_count(m), nt = qaws_mesh_3d_get_triangle_count(m), ne = qaws_mesh_3d_get_edge_count(m), i;
	double* x = ms_positions(m);
	double* n = NULL;
	unsigned int* t = (unsigned int*)malloc(sizeof(unsigned int) * 3 * (nt ? nt : 1));
	unsigned int* e = (unsigned int*)malloc(sizeof(unsigned int) * 2 * (ne ? ne : 1));
	qaws_mesh_3d_get_triangles(t, m);
	if (smooth)
	{
		qaws_scalar* ns = (qaws_scalar*)malloc(sizeof(qaws_scalar) * 3 * (nv ? nv : 1));
		n = (double*)malloc(sizeof(double) * 3 * (nv ? nv : 1));
		qaws_mesh_3d_get_vertex_normals(ns, m);
		for (i = 0; i < 3 * nv; i++)
			n[i] = (double)ns[i];
		free(ns);
	}
	for (i = 0; (rgb || vrgb) && i < nt; i++)   /* (no colour: the edges only) */
	{
		unsigned int a = t[3 * i], b = t[3 * i + 1], d = t[3 * i + 2];
		er_triangle(c, x + 3 * a, x + 3 * b, x + 3 * d, n ? n + 3 * a : NULL, n ? n + 3 * b : NULL, n ? n + 3 * d : NULL,
			vrgb ? vrgb + 3 * a : rgb, vrgb ? vrgb + 3 * b : rgb, vrgb ? vrgb + 3 * d : rgb);
	}
	if (wire)
	{
		qaws_mesh_3d_get_edges(e, m);
		for (i = 0; i < ne; i++)
			er_line(c, x + 3 * e[2 * i], x + 3 * e[2 * i + 1], wire, wire_width);
	}
	free(x); free(n); free(t); free(e);
}

/* a panel fitted to a mesh, with its title and a line of numbers */
static void ms_panel(er_canvas* c, int x, int y, int w, int h, qaws_mesh_3d const* m, double yaw, double pitch, char const* title,
	char const* sub)
{
	double* p = ms_positions(m);
	er_panel(c, x, y + 72, w, h - 72, p, qaws_mesh_3d_get_vertex_count(m), yaw, pitch, 0.06);
	free(p);
	er_text(c, x + 14, y + 12, title, 3, 30, 30, 30);
	if (sub)
		er_text(c, x + 14, y + 46, sub, 2, 95, 95, 95);
}

#include "meshing_showcase/01_core.c"
#include "meshing_showcase/02_reference.c"
#include "meshing_showcase/03_remesh.c"
#include "meshing_showcase/04_quads.c"
#include "meshing_showcase/05_generate.c"
#include "meshing_showcase/06_tets.c"
#include "meshing_showcase/07_hexes.c"
#include "meshing_showcase/08_frame_hexes.c"
#include "meshing_showcase/09_measures.c"
#include "meshing_showcase/10_remesh_angles.c"
#include "meshing_showcase/11_remesh_cvt.c"
#include "meshing_showcase/12_cross_field.c"
#include "meshing_showcase/13_quad_flow.c"
#include "meshing_showcase/14_quad_improve.c"
#include "meshing_showcase/15_quad_global.c"
#include "meshing_showcase/16_guidance.c"
#include "meshing_showcase/17_feature_size.c"
#include "meshing_showcase/18_repair.c"
#include "meshing_showcase/19_layers.c"
#include "meshing_showcase/20_curved.c"
#include "meshing_showcase/21_metric.c"
#include "meshing_showcase/22_bvh.c"
#include "meshing_showcase/23_quality.c"
#include "meshing_showcase/24_repair_more.c"
#include "meshing_showcase/25_watertight.c"
#include "meshing_showcase/26_intersect.c"
#include "meshing_showcase/27_boolean.c"
#include "meshing_showcase/28_geodesic.c"
#include "meshing_showcase/29_smooth.c"
#include "meshing_showcase/30_deform.c"
#include "meshing_showcase/31_deform_2.c"
#include "meshing_showcase/32_decimate.c"
#include "meshing_showcase/33_decimate_views.c"
#include "meshing_showcase/34_remesh_attributes.c"
#include "meshing_showcase/35_connect.c"
#include "meshing_showcase/36_connect_more.c"
#include "meshing_showcase/37_repulsive_mesh_curves.c"
#include "meshing_showcase/38_repulsive_mesh_curves_2.c"
#include "meshing_showcase/39_offset.c"
#include "meshing_showcase/40_subdivide.c"
#include "meshing_showcase/41_volume_improve.c"
#include "meshing_showcase/42_volume_curve.c"
#include "meshing_showcase/43_tetrahedra_border.c"
#include "meshing_showcase/44_polycube.c"
#include "meshing_showcase/45_quads_pairing.c"
#include "meshing_showcase/46_quad_layout.c"
#include "meshing_showcase/47_denoise.c"
#include "meshing_showcase/48_points.c"
#include "meshing_showcase/49_geodesic_exact.c"
#include "meshing_showcase/50_butterfly.c"
#include "meshing_showcase/51_decimate_attributes.c"
#include "meshing_showcase/52_offset_round.c"
#include "meshing_showcase/53_volume_cubic.c"
#include "meshing_showcase/54_subdivide_adaptive.c"

/* the figures to draw: all, or those whose keys the first argument lists ("58"; figures 10 to 35 are a to z, 36 to 54 are A to S) */
int main(int argc, char** argv)
{
	static void (*const demo[])(void) = { demo_core, demo_reference, demo_remesh, demo_quads, demo_generate, demo_tets, demo_hexes, demo_frame_hexes,
		demo_measures, demo_remesh_angles, demo_remesh_cvt, demo_cross_field, demo_quad_flow, demo_quad_improve,
		demo_quad_global, demo_guidance, demo_feature_size, demo_repair, demo_layers, demo_curved, demo_metric, demo_bvh, demo_quality, demo_repair_more, demo_watertight,
		demo_intersect, demo_boolean, demo_geodesic, demo_smooth, demo_deform, demo_deform_2, demo_decimate, demo_decimate_views, demo_remesh_attributes, demo_connect, demo_connect_more,
		demo_repulsive_mesh_curves, demo_repulsive_mesh_curves_2, demo_offset, demo_subdivide, demo_volume_improve, demo_volume_curve, demo_tetrahedra_border, demo_polycube, demo_quads_pairing, demo_quad_layout, demo_denoise, demo_points, demo_geodesic_exact, demo_butterfly, demo_decimate_attributes, demo_offset_round, demo_volume_cubic, demo_subdivide_adaptive };
	static char const key[] = "123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRS";
	unsigned int i;
	setvbuf(stdout, NULL, _IONBF, 0);
	MAKE_DIR("showcase");
	for (i = 0; i < sizeof demo / sizeof *demo; i++)
		if (argc < 2 || strchr(argv[1], key[i]))
			demo[i]();
	return 0;
}
