/*
 * GPU parity of the tangent-point kernels (examples/repulsive.hlsl, R6):
 * Direct3D 11 compute, Windows only, outside the library (no runtime
 * dependency). Each case runs the kernels on the GPU (float) and compares
 * them with the CPU library (double): the energy, its gradient, and the
 * metric's product (small cases: against the dense metric). Then the time
 * of one energy and gradient against the CPU direct evaluator and
 * Barnes-Hut.
 *
 * Build (f64 library):
 *   cl /O2 /MD /DQAWS_SCALAR_IS_FLOAT=0 /Isrc/qaws examples/repulsive_gpu_parity.c qaws.lib d3d11.lib d3dcompiler.lib
 * Run from the repository root (it compiles examples/repulsive.hlsl).
 */

#define COBJMACROS
#include <windows.h>
#include <d3d11.h>
#include <d3dcompiler.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_repulsive.h"
#include "internal/qaws_internal_repulsive.h"

typedef struct gp_ctx
{
	ID3D11Device* dev;
	ID3D11DeviceContext* ctx;
	ID3D11ComputeShader *element, *gather;
} gp_ctx;

typedef struct gp_params
{
	unsigned int element_count, vertex_count, k, pad0;
	float alpha, beta, ex, pad1;
} gp_params;

static double gp_now(void)
{
	LARGE_INTEGER f, c;
	QueryPerformanceFrequency(&f);
	QueryPerformanceCounter(&c);
	return (double)c.QuadPart / (double)f.QuadPart;
}

/* includes resolved against the directories the sources live in (the
   standard handler resolves nested relative includes from the root file) */
static HRESULT __stdcall gp_open(ID3DInclude* self, D3D_INCLUDE_TYPE type, LPCSTR name, LPCVOID parent, LPCVOID* data, UINT* bytes)
{
	static char const* const dirs[3] = { "examples/", "src/qaws/", "src/qaws/core/" };
	unsigned int d;
	(void)self; (void)type; (void)parent;
	for (d = 0; d < 3; d++)
	{
		char path[512];
		FILE* fp;
		long n;
		char* buf;
		snprintf(path, sizeof(path), "%s%s", dirs[d], name);
		if (!(fp = fopen(path, "rb")))
			continue;
		fseek(fp, 0, SEEK_END); n = ftell(fp); fseek(fp, 0, SEEK_SET);
		buf = (char*)malloc(n + 1);
		fread(buf, 1, n, fp);
		fclose(fp);
		*data = buf; *bytes = (UINT)n;
		return S_OK;
	}
	return E_FAIL;
}

static HRESULT __stdcall gp_close(ID3DInclude* self, LPCVOID data)
{
	(void)self;
	free((void*)data);
	return S_OK;
}

static ID3DIncludeVtbl gp_include_vtbl = { gp_open, gp_close };
static ID3DInclude gp_include = { &gp_include_vtbl };

static ID3D11ComputeShader* gp_shader(gp_ctx* g, char const* entry)
{
	ID3DBlob *blob = NULL, *err = NULL;
	ID3D11ComputeShader* cs = NULL;
	HRESULT hr = D3DCompileFromFile(L"examples/repulsive.hlsl", NULL, &gp_include, entry, "cs_5_0", D3DCOMPILE_OPTIMIZATION_LEVEL3, 0, &blob, &err);
	if (FAILED(hr))
	{
		printf("compile %s failed: %s\n", entry, err ? (char const*)ID3D10Blob_GetBufferPointer(err) : "?");
		return NULL;
	}
	ID3D11Device_CreateComputeShader(g->dev, ID3D10Blob_GetBufferPointer(blob), ID3D10Blob_GetBufferSize(blob), NULL, &cs);
	ID3D10Blob_Release(blob);
	if (err) ID3D10Blob_Release(err);
	return cs;
}

/* a structured buffer of n 4-byte elements: with data an SRV, without a UAV */
static ID3D11Buffer* gp_buffer(gp_ctx* g, void const* data, unsigned int n, ID3D11ShaderResourceView** srv, ID3D11UnorderedAccessView** uav)
{
	D3D11_BUFFER_DESC bd;
	D3D11_SUBRESOURCE_DATA sd;
	ID3D11Buffer* b = NULL;
	memset(&bd, 0, sizeof(bd));
	bd.ByteWidth = 4 * (n ? n : 1);
	bd.Usage = D3D11_USAGE_DEFAULT;
	bd.BindFlags = srv ? D3D11_BIND_SHADER_RESOURCE : D3D11_BIND_UNORDERED_ACCESS;
	bd.MiscFlags = D3D11_RESOURCE_MISC_BUFFER_STRUCTURED;
	bd.StructureByteStride = 4;
	sd.pSysMem = data; sd.SysMemPitch = 0; sd.SysMemSlicePitch = 0;
	ID3D11Device_CreateBuffer(g->dev, &bd, data ? &sd : NULL, &b);
	if (srv) ID3D11Device_CreateShaderResourceView(g->dev, (ID3D11Resource*)b, NULL, srv);
	if (uav) ID3D11Device_CreateUnorderedAccessView(g->dev, (ID3D11Resource*)b, NULL, uav);
	return b;
}

static void gp_read(gp_ctx* g, ID3D11Buffer* b, unsigned int n, float* out)
{
	D3D11_BUFFER_DESC bd;
	D3D11_MAPPED_SUBRESOURCE m;
	ID3D11Buffer* st = NULL;
	ID3D11Buffer_GetDesc(b, &bd);
	bd.Usage = D3D11_USAGE_STAGING; bd.BindFlags = 0; bd.CPUAccessFlags = D3D11_CPU_ACCESS_READ; bd.MiscFlags = 0;
	ID3D11Device_CreateBuffer(g->dev, &bd, NULL, &st);
	ID3D11DeviceContext_CopyResource(g->ctx, (ID3D11Resource*)st, (ID3D11Resource*)b);
	ID3D11DeviceContext_Map(g->ctx, (ID3D11Resource*)st, 0, D3D11_MAP_READ, 0, &m);
	memcpy(out, m.pData, 4 * n);
	ID3D11DeviceContext_Unmap(g->ctx, (ID3D11Resource*)st, 0);
	ID3D11Buffer_Release(st);
}

/* the kernels on x (dim 3 per vertex): energy, gradient, y = A u; the time of one run (after a warm-up) */
static double gp_run(gp_ctx* g, double const* x, unsigned int nv, unsigned int const* idx, unsigned int ne, unsigned int k, double alpha, double beta,
	double const* u, double* E, double* grad, double* y)
{
	float *xf = (float*)malloc(4 * 3 * nv), *uf = (float*)malloc(4 * nv), *sh = (float*)malloc(4 * ne), *gf = (float*)malloc(4 * 3 * nv), *yf = (float*)malloc(4 * nv);
	unsigned int *off = (unsigned int*)calloc(nv + 1, 4), *ent = (unsigned int*)malloc(4 * k * ne), *fill = (unsigned int*)calloc(nv + 1, 4), i, a, rep;
	ID3D11ShaderResourceView* srv[5];
	ID3D11UnorderedAccessView* uav[6];
	ID3D11Buffer* buf[11];
	ID3D11Buffer* cb = NULL;
	D3D11_BUFFER_DESC bd;
	D3D11_SUBRESOURCE_DATA sd;
	gp_params prm;
	double t = 0;
	for (i = 0; i < 3 * nv; i++) xf[i] = (float)x[i];
	for (i = 0; i < nv; i++) uf[i] = (float)u[i];
	for (i = 0; i < k * ne; i++) off[idx[i] + 1]++;
	for (i = 0; i < nv; i++) off[i + 1] += off[i];
	for (i = 0; i < ne; i++)
		for (a = 0; a < k; a++) ent[off[idx[k * i + a]] + fill[idx[k * i + a]]++] = k * i + a;
	buf[0] = gp_buffer(g, xf, 3 * nv, &srv[0], NULL);
	buf[1] = gp_buffer(g, idx, k * ne, &srv[1], NULL);
	buf[2] = gp_buffer(g, off, nv + 1, &srv[2], NULL);
	buf[3] = gp_buffer(g, ent, k * ne, &srv[3], NULL);
	buf[4] = gp_buffer(g, uf, nv, &srv[4], NULL);
	buf[5] = gp_buffer(g, NULL, ne, NULL, &uav[0]);
	buf[6] = gp_buffer(g, NULL, 9 * ne, NULL, &uav[1]);
	buf[7] = gp_buffer(g, NULL, 4 * ne, NULL, &uav[2]);
	buf[8] = gp_buffer(g, NULL, 9 * ne, NULL, &uav[3]);
	buf[9] = gp_buffer(g, NULL, 3 * nv, NULL, &uav[4]);
	buf[10] = gp_buffer(g, NULL, nv, NULL, &uav[5]);
	prm.element_count = ne; prm.vertex_count = nv; prm.k = k; prm.pad0 = 0;
	prm.alpha = (float)alpha; prm.beta = (float)beta;
	prm.ex = (float)(k == 2 ? 2 * ((beta - 1) / alpha - 1) + 1 : 2 * (1 - 2 / alpha) + 2); prm.pad1 = 0;
	memset(&bd, 0, sizeof(bd));
	bd.ByteWidth = sizeof(prm); bd.Usage = D3D11_USAGE_DEFAULT; bd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	sd.pSysMem = &prm; sd.SysMemPitch = 0; sd.SysMemSlicePitch = 0;
	ID3D11Device_CreateBuffer(g->dev, &bd, &sd, &cb);
	for (rep = 0; rep < 2; rep++)
	{
		double t0 = gp_now();
		ID3D11DeviceContext_CSSetShaderResources(g->ctx, 0, 5, srv);
		ID3D11DeviceContext_CSSetUnorderedAccessViews(g->ctx, 0, 6, uav, NULL);
		ID3D11DeviceContext_CSSetConstantBuffers(g->ctx, 0, 1, &cb);
		ID3D11DeviceContext_CSSetShader(g->ctx, g->element, NULL, 0);
		ID3D11DeviceContext_Dispatch(g->ctx, (ne + 63) / 64, 1, 1);
		ID3D11DeviceContext_CSSetShader(g->ctx, g->gather, NULL, 0);
		ID3D11DeviceContext_Dispatch(g->ctx, (nv + 63) / 64, 1, 1);
		gp_read(g, buf[5], ne, sh);
		gp_read(g, buf[9], 3 * nv, gf);
		gp_read(g, buf[10], nv, yf);
		t = gp_now() - t0;
	}
	*E = 0;
	for (i = 0; i < ne; i++) *E += (double)sh[i];
	for (i = 0; i < 3 * nv; i++) grad[i] = (double)gf[i];
	for (i = 0; y && i < nv; i++) y[i] = (double)yf[i];
	for (i = 0; i < 5; i++) ID3D11ShaderResourceView_Release(srv[i]);
	for (i = 0; i < 6; i++) ID3D11UnorderedAccessView_Release(uav[i]);
	for (i = 0; i < 11; i++) ID3D11Buffer_Release(buf[i]);
	ID3D11Buffer_Release(cb);
	free(xf); free(uf); free(sh); free(gf); free(yf); free(off); free(ent); free(fill);
	return t;
}

static double gp_rel(double const* a, double const* b, unsigned int n)
{
	double e = 0, m = 0;
	unsigned int i;
	for (i = 0; i < n; i++) { e = fmax(e, fabs(a[i] - b[i])); m = fmax(m, fabs(b[i])); }
	return e / m;
}

/* a lumpy trefoil of n edges */
static void gp_trefoil(unsigned int n, double* x, unsigned int* e)
{
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double t = 2 * 3.14159265358979323846 * i / n;
		x[3 * i] = sin(t) + 2 * sin(2 * t) + 0.13 * cos(5 * t); x[3 * i + 1] = cos(t) - 2 * cos(2 * t) + 0.07 * sin(4 * t); x[3 * i + 2] = -0.6 * sin(3 * t);
		e[2 * i] = i; e[2 * i + 1] = (i + 1) % n;
	}
}

/* an octahedron subdivided `level` times, lumpy; x holds 3 (2 + 4^(level+1)) */
static void gp_sphere(unsigned int level, double* x, unsigned int* tris, unsigned int* nv, unsigned int* nf)
{
	static double const ov[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	static unsigned int const of[8][3] = { { 0, 2, 4 }, { 2, 1, 4 }, { 1, 3, 4 }, { 3, 0, 4 }, { 2, 0, 5 }, { 1, 2, 5 }, { 3, 1, 5 }, { 0, 3, 5 } };
	unsigned int l, i, n = 6, f = 8, q, k;
	unsigned int* t2 = (unsigned int*)malloc(4 * 3 * 8 * ((size_t)1 << (2 * level)));
	memcpy(x, ov, sizeof(ov));
	memcpy(tris, of, sizeof(of));
	for (l = 0; l < level; l++)
	{
		unsigned int nf2 = 0, n0 = n;
		for (i = 0; i < f; i++)
		{
			unsigned int c[3], mid[3];
			for (k = 0; k < 3; k++) c[k] = tris[3 * i + k];
			for (k = 0; k < 3; k++)
			{
				unsigned int a = c[k], b = c[(k + 1) % 3], j;
				double m[3], len;
				for (q = 0; q < 3; q++) m[q] = 0.5 * (x[3 * a + q] + x[3 * b + q]);
				len = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
				for (q = 0; q < 3; q++) m[q] /= len;
				/* (only this level's midpoints can match) */
				for (j = n0; j < n && !(fabs(x[3 * j] - m[0]) < 1e-12 && fabs(x[3 * j + 1] - m[1]) < 1e-12 && fabs(x[3 * j + 2] - m[2]) < 1e-12); j++)
					;
				if (j == n) { for (q = 0; q < 3; q++) x[3 * n + q] = m[q]; n++; }
				mid[k] = j;
			}
			t2[3 * nf2] = c[0]; t2[3 * nf2 + 1] = mid[0]; t2[3 * nf2 + 2] = mid[2]; nf2++;
			t2[3 * nf2] = mid[0]; t2[3 * nf2 + 1] = c[1]; t2[3 * nf2 + 2] = mid[1]; nf2++;
			t2[3 * nf2] = mid[2]; t2[3 * nf2 + 1] = mid[1]; t2[3 * nf2 + 2] = c[2]; nf2++;
			t2[3 * nf2] = mid[0]; t2[3 * nf2 + 1] = mid[1]; t2[3 * nf2 + 2] = mid[2]; nf2++;
		}
		memcpy(tris, t2, 4 * 3 * (size_t)nf2);
		f = nf2;
	}
	for (i = 0; i < n; i++)
	{
		double r = 1 + 0.15 * sin(3 * x[3 * i]) * cos(2 * x[3 * i + 1]);
		for (q = 0; q < 3; q++) x[3 * i + q] *= r;
	}
	free(t2);
	*nv = n; *nf = f;
}

/* one case: parity (with the dense metric when small) and times */
static void gp_case(gp_ctx* g, char const* name, double const* x0, unsigned int nv, unsigned int const* idx, unsigned int ne, unsigned int k, int dense)
{
	double *u = (double*)malloc(8 * nv), *gg = (double*)malloc(8 * 3 * nv), *yg = (double*)malloc(8 * nv), *gc = (double*)calloc(3 * nv, 8), *yc = (double*)calloc(nv, 8);
	double Eg, Ec = 0, Eb = 0, tg, tc, tb, *gb = (double*)calloc(3 * nv, 8);
	double alpha = k == 2 ? 3 : 6, beta = k == 2 ? 6 : 12, t0;
	qaws_repulsive_curves* c = NULL;
	qaws_repulsive_mesh* m = NULL;
	double* x = (double*)malloc(8 * 3 * nv);
	double sens = 0;
	unsigned int i, j;
	for (i = 0; i < nv; i++) u[i] = sin(0.37 * i) + 0.5 * cos(1.3 * i);
	/* the GPU sees the positions rounded to float: the reference does too */
	for (i = 0; i < 3 * nv; i++) x[i] = (double)(float)x0[i];
	tg = gp_run(g, x, nv, idx, ne, k, alpha, beta, u, &Eg, gg, yg);
	if (k == 2)
	{
		qaws_repulsive_curves_desc d;
		memset(&d, 0, sizeof(d));
		d.dimension = 3; d.vertex_count = nv; d.edges = idx; d.edge_count = ne;
		qaws_repulsive_curves_create(&c, &d);
		t0 = gp_now(); qaws_internal_repulsive_curves_energy(&Ec, gc, c, x); tc = gp_now() - t0;
		t0 = gp_now(); qaws_internal_repulsive_curves_bh(&Eb, gb, c, x, 0.25); tb = gp_now() - t0;
		{
			/* how far that rounding alone moves the gradient (the problem's conditioning, in double) */
			double E0, *g0 = (double*)calloc(3 * nv, 8);
			qaws_internal_repulsive_curves_energy(&E0, g0, c, x0);
			sens = gp_rel(gc, g0, 3 * nv);
			free(g0);
		}
	}
	else
	{
		qaws_repulsive_mesh_desc d;
		memset(&d, 0, sizeof(d));
		d.vertex_count = nv; d.triangles = idx; d.triangle_count = ne;
		qaws_repulsive_mesh_create(&m, &d);
		t0 = gp_now(); qaws_internal_repulsive_mesh_energy(&Ec, gc, m, x); tc = gp_now() - t0;
		t0 = gp_now(); qaws_internal_repulsive_mesh_bh(&Eb, gb, m, x, 0.25); tb = gp_now() - t0;
		{
			double E0, *g0 = (double*)calloc(3 * nv, 8);
			qaws_internal_repulsive_mesh_energy(&E0, g0, m, x0);
			sens = gp_rel(gc, g0, 3 * nv);
			free(g0);
		}
	}
	printf("%s: GPU (float) against the CPU (double): energy %.2g, gradient %.2g", name, fabs(Eg - Ec) / Ec, gp_rel(gg, gc, 3 * nv));
	if (dense)
	{
		double* A = (double*)malloc(8 * (size_t)nv * nv);
		if (k == 2) qaws_internal_repulsive_curves_metric(A, c, x); else qaws_internal_repulsive_mesh_metric(A, m, x);
		for (i = 0; i < nv; i++)
			for (j = 0; j < nv; j++) yc[i] += A[(size_t)i * nv + j] * u[j];
		printf(", metric product %.2g", gp_rel(yg, yc, nv));
		free(A);
	}
	printf("\n    rounding the positions to float alone moves the double gradient by %.2g (the problem's conditioning: compare inputs alike)", sens);
	printf("\n    one energy and gradient: GPU (all pairs, with the metric product) %.4f s, CPU direct %.3f s, CPU Barnes-Hut (theta 0.25) %.3f s\n", tg, tc, tb);
	qaws_repulsive_curves_destroy(c); qaws_repulsive_mesh_destroy(m);
	free(u); free(gg); free(yg); free(gc); free(yc); free(gb); free(x);
}

int main(void)
{
	gp_ctx g;
	D3D_FEATURE_LEVEL fl;
	static unsigned int const curve_sizes[3] = { 300, 2000, 8000 };
	unsigned int s;
	IDXGIFactory1* fac = NULL;
	IDXGIAdapter1 *best = NULL, *a1 = NULL;
	DXGI_ADAPTER_DESC1 d1, bd1;
	SIZE_T mem = 0;
	UINT i;
	memset(&g, 0, sizeof(g));
	/* the adapter with the most dedicated memory (a discrete GPU over an integrated one) */
	if (SUCCEEDED(CreateDXGIFactory1(&IID_IDXGIFactory1, (void**)&fac)))
		for (i = 0; IDXGIFactory1_EnumAdapters1(fac, i, &a1) != DXGI_ERROR_NOT_FOUND; i++)
		{
			IDXGIAdapter1_GetDesc1(a1, &d1);
			if (!(d1.Flags & DXGI_ADAPTER_FLAG_SOFTWARE) && d1.DedicatedVideoMemory > mem)
			{
				if (best) IDXGIAdapter1_Release(best);
				best = a1; mem = d1.DedicatedVideoMemory; bd1 = d1;
			}
			else
				IDXGIAdapter1_Release(a1);
		}
	if (FAILED(D3D11CreateDevice((IDXGIAdapter*)best, best ? D3D_DRIVER_TYPE_UNKNOWN : D3D_DRIVER_TYPE_HARDWARE, NULL, 0, NULL, 0, D3D11_SDK_VERSION, &g.dev, &fl,
		&g.ctx)))
	{
		printf("no Direct3D 11 device\n");
		return 1;
	}
	if (best) printf("adapter: %ls\n", bd1.Description);
	g.element = gp_shader(&g, "CSElement");
	g.gather = gp_shader(&g, "CSGather");
	if (!g.element || !g.gather)
		return 1;
	for (s = 0; s < 3; s++)
	{
		unsigned int n = curve_sizes[s];
		double* x = (double*)malloc(8 * 3 * n);
		unsigned int* e = (unsigned int*)malloc(4 * 2 * n);
		char name[64];
		gp_trefoil(n, x, e);
		sprintf(name, "trefoil, %u edges", n);
		gp_case(&g, name, x, n, e, n, 2, n <= 2000);
		free(x); free(e);
	}
	for (s = 3; s <= 5; s++)
	{
		unsigned int nv, nf;
		double* x = (double*)malloc(8 * 3 * (2 + ((size_t)4 << (2 * s))));
		unsigned int* t = (unsigned int*)malloc(4 * 3 * 8 * ((size_t)1 << (2 * s)));
		char name[64];
		gp_sphere(s, x, t, &nv, &nf);
		sprintf(name, "lumpy sphere, %u faces", nf);
		gp_case(&g, name, x, nv, t, nf, 3, nf <= 2048);
		free(x); free(t);
	}
	return 0;
}
