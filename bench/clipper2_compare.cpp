// Clipper2 against qaws_clip64 (exact int64) and qaws_clip (curves, on
// polylines): the same inputs, timings, path counts and areas.
//
// Built only with -DQAWS_BUILD_CLIPPER2_BENCH=ON -DQAWS_CLIPPER2_DIR=<Clipper2
// checkout>; Clipper2 is compiled into this program, never into qaws.
//
//   1. Clipper2's "complex polygons" benchmark: two random polygons of n
//      vertices in 800 x 600 (each crossing itself all over), intersection,
//      non-zero
//   2. the union of many small random squares, all subjects
//   3. Clipper2's test corpus (Polygons.txt), every record

#include "clipper2/clipper.h"

extern "C" {
#include "qaws.h"
}

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

using namespace Clipper2Lib;

static double now()
{
	return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static unsigned long long g_rng = 0x9E3779B97F4A7C15ull;
static int rnd(int n)
{
	g_rng ^= g_rng << 13; g_rng ^= g_rng >> 7; g_rng ^= g_rng << 17;
	return (int)(g_rng % (unsigned long long)n);
}

struct Timing { double t; size_t paths; double area; };

template <class F> static Timing best_of(int reps, F f)
{
	Timing best = { 1e30, 0, 0 };
	for (int i = 0; i < reps; i++)
	{
		double t0 = now();
		Timing r = f();
		r.t = now() - t0;
		if (r.t < best.t) best = r;
	}
	return best;
}

static Timing run_clipper(ClipType ct, FillRule fr, Paths64 const& s, Paths64 const& c)
{
	Paths64 sol = BooleanOp(ct, fr, s, c);
	Timing r = { 0, sol.size(), Area(sol) };
	return r;
}

static Timing run_clip64(qaws_clip_type ct, qaws_fill_rule fr, Paths64 const& s, Paths64 const& c)
{
	std::vector<std::vector<int64_t>> store;
	std::vector<qaws_path64> sp, cp;
	for (auto const* set : { &s, &c })
		for (auto const& p : *set)
		{
			std::vector<int64_t> v;
			for (auto const& q : p) { v.push_back(q.x); v.push_back(q.y); }
			store.push_back(v);
		}
	size_t k = 0;
	for (auto const& p : s) { qaws_path64 q = { store[k].data(), (unsigned)p.size() }; sp.push_back(q); k++; }
	for (auto const& p : c) { qaws_path64 q = { store[k].data(), (unsigned)p.size() }; cp.push_back(q); k++; }
	qaws_clip64_desc d;
	memset(&d, 0, sizeof(d));
	d.subjects = sp.data(); d.subject_count = (unsigned)sp.size();
	d.clips = cp.data(); d.clip_count = (unsigned)cp.size();
	d.clip_type = ct; d.fill_rule = fr;
	qaws_clip64_result* r = nullptr;
	Timing t = { 0, 0, 0 };
	if (qaws_clip64_execute(&d, &r) == QAWS_STATUS_OK)
	{
		unsigned n = qaws_clip64_result_get_path_count(r);
		t.paths = n;
		for (unsigned i = 0; i < n; i++)
		{
			int64_t const* p; unsigned m;
			qaws_clip64_result_get_path(r, i, &p, &m);
			t.area += 0.5 * qaws_path64_area2(p, m);
		}
	}
	qaws_clip64_result_destroy(r);
	return t;
}

static Timing run_clip(qaws_clip_type ct, qaws_fill_rule fr, Paths64 const& s, Paths64 const& c)
{
	std::vector<qaws_curve*> curves;
	std::vector<qaws_curve const*> views;
	std::vector<qaws_path_2d> sp, cp;
	for (auto const* set : { &s, &c })
		for (auto const& p : *set)
		{
			std::vector<qaws_scalar> v;
			for (auto const& q : p) { v.push_back((qaws_scalar)q.x); v.push_back((qaws_scalar)q.y); }
			qaws_curve* cv = nullptr;
			qaws_curve_create_polyline_2d(v.data(), (unsigned)p.size(), 1, &cv);
			curves.push_back(cv);
		}
	views.assign(curves.begin(), curves.end());
	for (size_t i = 0; i < views.size(); i++)
	{
		qaws_path_2d q = { &views[i], 1, 1 };
		(i < s.size() ? sp : cp).push_back(q);
	}
	qaws_clip_desc d;
	memset(&d, 0, sizeof(d));
	d.subjects = sp.data(); d.subject_count = (unsigned)sp.size();
	d.clips = cp.data(); d.clip_count = (unsigned)cp.size();
	d.clip_type = ct; d.fill_rule = fr;
	qaws_clip_result* r = nullptr;
	Timing t = { 0, 0, 0 };
	if (qaws_clip_execute(&d, &r) == QAWS_STATUS_OK)
	{
		unsigned n = qaws_clip_result_get_path_count(r);
		t.paths = n;
		for (unsigned i = 0; i < n; i++)
		{
			qaws_path_2d p; qaws_scalar a = 0;
			qaws_clip_result_get_path(r, i, &p);
			qaws_path_compute_area_2d(&p, &a);
			t.area += a;
		}
	}
	qaws_clip_result_destroy(r);
	for (auto* c2 : curves) qaws_curve_destroy(c2);
	return t;
}

static Path64 random_poly(int w, int h, unsigned n)
{
	Path64 p;
	for (unsigned i = 0; i < n; i++) p.push_back(Point64(rnd(w), rnd(h)));
	return p;
}

static void row(char const* name, Timing const& c2, Timing const& e, Timing const& f)
{
	std::printf("| %s | %.4f s, %zu, %.0f | %.4f s (x%.1f), %zu, %.0f | %.4f s (x%.1f), %zu, %.0f |\n", name,
		c2.t, c2.paths, c2.area, e.t, e.t / c2.t, e.paths, e.area, f.t, f.t / c2.t, f.paths, f.area);
	std::fflush(stdout);
}

int main(int argc, char** argv)
{
	int reps = 3;
	std::printf("| case | Clipper2 Clipper64 | qaws_clip64 (exact) | qaws_clip (polylines) |\n|---|---|---|---|\n");
	// 1. Clipper2's complex polygons benchmark
	for (unsigned n : { 50u, 100u, 200u, 400u })
	{
		Paths64 s = { random_poly(800, 600, n) }, c = { random_poly(800, 600, n) };
		char name[64];
		std::snprintf(name, sizeof(name), "complex polygons, %u vertices each", n);
		Timing a = best_of(reps, [&] { return run_clipper(ClipType::Intersection, FillRule::NonZero, s, c); });
		Timing b = best_of(reps, [&] { return run_clip64(QAWS_CLIP_INTERSECTION, QAWS_FILL_NON_ZERO, s, c); });
		Timing f = best_of(reps, [&] { return run_clip(QAWS_CLIP_INTERSECTION, QAWS_FILL_NON_ZERO, s, c); });
		row(name, a, b, f);
	}
	// 2. union of many small squares
	for (unsigned n : { 500u, 2000u, 8000u })
	{
		Paths64 s, c;
		for (unsigned i = 0; i < n; i++)
		{
			int x = rnd(4000), y = rnd(4000), w = 20 + rnd(60);
			s.push_back(MakePath({ x, y, x + w, y, x + w, y + w, x, y + w }));
		}
		char name[64];
		std::snprintf(name, sizeof(name), "union of %u squares", n);
		Timing a = best_of(reps, [&] { return run_clipper(ClipType::Union, FillRule::NonZero, s, c); });
		Timing b = best_of(reps, [&] { return run_clip64(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, s, c); });
		Timing f = best_of(reps, [&] { return run_clip(QAWS_CLIP_UNION, QAWS_FILL_NON_ZERO, s, c); });
		row(name, a, b, f);
	}
	// 3. the corpus
	if (argc > 1)
	{
		std::ifstream in(argv[1]);
		std::string line;
		std::vector<std::pair<std::pair<ClipType, FillRule>, std::pair<Paths64, Paths64>>> recs;
		ClipType ct = ClipType::Union; FillRule fr = FillRule::NonZero;
		Paths64 s, c; int section = 0;
		auto flush = [&] { if (!s.empty() || !c.empty()) recs.push_back({ { ct, fr }, { s, c } }); s.clear(); c.clear(); };
		while (std::getline(in, line))
		{
			if (line.rfind("CAPTION:", 0) == 0) { flush(); section = 0; }
			else if (line.rfind("CLIPTYPE:", 0) == 0)
				ct = line.find("INTERSECTION") != std::string::npos ? ClipType::Intersection : line.find("UNION") != std::string::npos ? ClipType::Union :
					line.find("DIFFERENCE") != std::string::npos ? ClipType::Difference : ClipType::Xor;
			else if (line.rfind("FILLRULE:", 0) == 0)
				fr = line.find("EVENODD") != std::string::npos ? FillRule::EvenOdd : line.find("NONZERO") != std::string::npos ? FillRule::NonZero :
					line.find("POSITIVE") != std::string::npos ? FillRule::Positive : FillRule::Negative;
			else if (line.rfind("SUBJECTS_OPEN", 0) == 0) section = 0;
			else if (line.rfind("SUBJECTS", 0) == 0) section = 1;
			else if (line.rfind("CLIPS", 0) == 0) section = 2;
			else if (section && line.find_first_of("0123456789") != std::string::npos)
			{
				Path64 p; const char* q = line.c_str(); char* e; std::vector<int64_t> v;
				while (*q) { long long x = std::strtoll(q, &e, 10); if (e == q) { q++; continue; } v.push_back(x); q = e; }
				for (size_t i = 0; i + 1 < v.size(); i += 2) p.push_back(Point64(v[i], v[i + 1]));
				(section == 1 ? s : c).push_back(p);
			}
		}
		flush();
		auto qct = [](ClipType t) { return t == ClipType::Intersection ? QAWS_CLIP_INTERSECTION : t == ClipType::Union ? QAWS_CLIP_UNION :
			t == ClipType::Difference ? QAWS_CLIP_DIFFERENCE : QAWS_CLIP_XOR; };
		auto qfr = [](FillRule f) { return f == FillRule::EvenOdd ? QAWS_FILL_EVEN_ODD : f == FillRule::NonZero ? QAWS_FILL_NON_ZERO :
			f == FillRule::Positive ? QAWS_FILL_POSITIVE : QAWS_FILL_NEGATIVE; };
		Timing a = best_of(reps, [&] { Timing t = { 0, 0, 0 }; for (auto& r : recs) { Timing u = run_clipper(r.first.first, r.first.second, r.second.first, r.second.second); t.paths += u.paths; t.area += u.area; } return t; });
		Timing b = best_of(reps, [&] { Timing t = { 0, 0, 0 }; for (auto& r : recs) { Timing u = run_clip64(qct(r.first.first), qfr(r.first.second), r.second.first, r.second.second); t.paths += u.paths; t.area += u.area; } return t; });
		Timing f = best_of(reps, [&] { Timing t = { 0, 0, 0 }; for (auto& r : recs) { Timing u = run_clip(qct(r.first.first), qfr(r.first.second), r.second.first, r.second.second); t.paths += u.paths; t.area += u.area; } return t; });
		char name[64];
		std::snprintf(name, sizeof(name), "Clipper2's Polygons.txt, %zu records", recs.size());
		row(name, a, b, f);
	}
	return 0;
}
