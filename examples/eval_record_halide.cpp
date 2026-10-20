/*
 * Curve and surface records evaluated by Halide pipelines, the knot spans found at run time (an Expr per element: the
 * parameters are read from a buffer), against the same record cores compiled for C (eval_record_halide_ref.c), f64.
 *
 * Halide is not part of the build (docs/capabilities.md). By hand, from a Visual Studio x64 prompt, with HALIDE the
 * root of a Halide install:
 *   cl /nologo /O2 /MD /std:c++17 /EHsc /I %HALIDE%\include examples\eval_record_halide.cpp examples\eval_record_halide_ref.c
 *      /Fe:working_dir\eval_record_halide.exe /link %HALIDE%\lib\Halide.lib
 *   copy %HALIDE%\bin\Halide.dll working_dir & working_dir\eval_record_halide.exe
 * It prints the largest difference of each record against C and returns non-zero past 1e-12 (relative).
 */

#define QAWS_SCALAR_IS_FLOAT 0
#include "../src/qaws/qaws_halide.h"

static Halide::Func g_knots, g_points, g_weights;
#define QAWS_BUFFER_KNOT(i) g_knots(i)
#define QAWS_BUFFER_POINT(i) g_points(i)
#define QAWS_BUFFER_WEIGHT(i) g_weights(i)
#include "../src/qaws/core/qaws_curve_record_core.h"
#include "../src/qaws/core/qaws_surface_record_core.h"

#include <cmath>
#include <cstdio>
#include <vector>

extern "C" void eval_record_ref_curve(double* out, int const* rec, double const* knots, double const* points, double const* weights,
	double const* t, int n);
extern "C" void eval_record_ref_surface(double* out, int const* rec, double const* knots, double const* points, double const* weights,
	double const* uv, int n);

/* the buffers as Funcs, clamped to their extent (the indices the data decide have no bounds Halide can infer; a valid
   index is left as it is) */
static void er_bind(Halide::Buffer<double>& k, Halide::Buffer<double>& p, Halide::Buffer<double>& w)
{
	Halide::Var x("x");
	g_knots = Halide::Func("knots");
	g_points = Halide::Func("points");
	g_weights = Halide::Func("weights");
	g_knots(x) = k(Halide::clamp(x, 0, k.width() - 1));
	g_points(x) = p(Halide::clamp(x, 0, p.width() - 1));
	g_weights(x) = w(Halide::clamp(x, 0, w.width() - 1));
}

static double er_compare(char const* name, std::vector<double> const& got, std::vector<double> const& ref)
{
	double worst = 0, scale = 0;
	for (size_t i = 0; i < ref.size(); i++)
		scale = std::fmax(scale, std::fabs(ref[i]));
	for (size_t i = 0; i < ref.size(); i++)
		worst = std::fmax(worst, std::fabs(got[i] - ref[i]));
	std::printf("  %-44s %6zu values, largest difference %.3g (%.3g of the largest value)\n", name, ref.size(), worst, worst / scale);
	return worst / scale;
}

int main()
{
	using namespace Halide;
	/* two curves in shared buffers: a clamped cubic B-spline with a double interior knot (3D), and a rational quadratic
	   (2D) whose knot repeated at the end leaves empty spans before it */
	std::vector<double> knots = { 0, 0, 0, 0, 0.1, 0.25, 0.25, 0.4, 0.55, 0.7, 0.9, 1, 1, 1, 1,
		0, 0, 0, 0.2, 0.35, 0.5, 0.75, 1, 1, 1, 1.2, 1.4 };
	std::vector<double> points, weights = { 1, 0.7, 1.3, 0.9, 1, 2, 0.5, 1, 1.1 };
	for (int i = 0; i < 11; i++)
	{
		points.push_back(std::cos(0.7 * i) * (1 + 0.1 * i));
		points.push_back(std::sin(0.9 * i) - 0.2 * i);
		points.push_back(0.3 * i - 0.05 * i * i);
	}
	for (int i = 0; i < 9; i++)
	{
		points.push_back(std::cos(0.8 * i) + 0.15 * i);
		points.push_back(std::sin(1.1 * i) * 0.8);
	}
	int const recs[2][7] = { { 3, 11, 3, 0, 0, 0, 0 }, { 2, 9, 2, 1, 15, 33, 0 } };
	double worst = 0;
	std::printf("Halide records against C (f64)\n");
	for (int c = 0; c < 2; c++)
	{
		int const* rc = recs[c];
		double t0 = knots[rc[4] + rc[0]], t1 = knots[rc[4] + rc[1]];
		std::vector<double> t;
		/* across the domain, every knot exactly, and a little outside (clamped) */
		for (int i = 0; i <= 256; i++) t.push_back(t0 - 0.01 * (t1 - t0) + 1.02 * (t1 - t0) * i / 256.0);
		for (int i = rc[0]; i <= rc[1]; i++) t.push_back(knots[rc[4] + i]);
		int n = (int)t.size();
		Buffer<double> kb(knots.data(), (int)knots.size()), pb(points.data(), (int)points.size()), wb(weights.data(), (int)weights.size());
		Buffer<double> tb(t.data(), n);
		er_bind(kb, pb, wb);
		qaws_curve_record r;
		r.degree = rc[0]; r.count = rc[1]; r.dimension = rc[2]; r.rational = rc[3];
		r.knot_offset = rc[4]; r.point_offset = rc[5]; r.weight_offset = rc[6]; r.reserved = 0;
		Var i("i");
		qaws_eval_3d e = qaws_record_eval_3d(r, qaws_scalar(Expr(tb(i))), 7);
		Func f("curve");
		f(i) = Tuple({ e.position.x, e.position.y, e.position.z, e.d1.x, e.d1.y, e.d1.z, e.d2.x, e.d2.y, e.d2.z });
		Realization R = f.realize({ n });
		std::vector<double> got(9 * n), ref(9 * n);
		for (int k = 0; k < 9; k++)
		{
			Buffer<double> b = R[k];
			for (int j = 0; j < n; j++) got[9 * j + k] = b(j);
		}
		eval_record_ref_curve(ref.data(), rc, knots.data(), points.data(), weights.data(), t.data(), n);
		worst = std::fmax(worst, er_compare(c ? "NURBS curve, empty spans at the end" : "B-spline curve, double interior knot", got, ref));
	}
	/* a surface, cubic in u and quadratic in v, non-uniform knots, rational and not */
	{
		std::vector<double> sk = { 0, 0, 0, 0, 0.3, 0.45, 1, 1, 1, 1, 0, 0, 0, 0.4, 0.6, 1, 1, 1 };
		std::vector<double> sp, sw;
		for (int a = 0; a < 6; a++)
			for (int b = 0; b < 5; b++)
			{
				sp.push_back(a * 0.4);
				sp.push_back(b * 0.5 + 0.05 * a);
				sp.push_back(0.3 * std::sin(1.3 * a + 0.7 * b));
				sw.push_back(1 + 0.4 * std::sin(a + 2.0 * b));
			}
		std::vector<double> uv;
		for (int a = 0; a <= 16; a++)
			for (int b = 0; b <= 16; b++)
			{
				uv.push_back(a / 16.0);
				uv.push_back(b / 16.0);
			}
		uv.push_back(0.3); uv.push_back(0.4);      /* on interior knots */
		uv.push_back(0.45); uv.push_back(0.6);
		int n = (int)uv.size() / 2;
		for (int rational = 0; rational < 2; rational++)
		{
			int const rc[9] = { 3, 2, 6, 5, rational, 0, 10, 0, 0 };
			Buffer<double> kb(sk.data(), (int)sk.size()), pb(sp.data(), (int)sp.size()), wb(sw.data(), (int)sw.size());
			Buffer<double> ub(uv.data(), 2, n);
			er_bind(kb, pb, wb);
			qaws_surface_record r;
			r.u_degree = rc[0]; r.v_degree = rc[1]; r.u_count = rc[2]; r.v_count = rc[3]; r.rational = rc[4];
			r.u_knot_offset = rc[5]; r.v_knot_offset = rc[6]; r.point_offset = rc[7]; r.weight_offset = rc[8];
			r.reserved0 = r.reserved1 = r.reserved2 = 0;
			Var i("i");
			qaws_surface_eval_3d e = qaws_surface_record_eval(r, qaws_scalar(Expr(ub(0, i))), qaws_scalar(Expr(ub(1, i))), 127);
			qaws_vec3 const* v[7] = { &e.position, &e.du, &e.dv, &e.duu, &e.duv, &e.dvv, &e.normal };
			std::vector<Expr> outs;
			for (int k = 0; k < 7; k++) { outs.push_back(v[k]->x); outs.push_back(v[k]->y); outs.push_back(v[k]->z); }
			Func f("surface");
			f(i) = Tuple(outs);
			Realization R = f.realize({ n });
			std::vector<double> got(21 * n), ref(21 * n);
			for (int k = 0; k < 21; k++)
			{
				Buffer<double> b = R[k];
				for (int j = 0; j < n; j++) got[21 * j + k] = b(j);
			}
			eval_record_ref_surface(ref.data(), rc, sk.data(), sp.data(), sw.data(), uv.data(), n);
			worst = std::fmax(worst, er_compare(rational ? "NURBS surface (3 x 2), with normals" : "B-spline surface (3 x 2), with normals", got, ref));
		}
	}
	std::printf("%s: largest relative difference %.3g\n", worst <= 1e-12 ? "OK" : "FAILED", worst);
	return worst <= 1e-12 ? 0 : 1;
}
