/*
 * GPU parity of the evaluation cores (examples/eval_parity.hlsl): Direct3D
 * 11 compute, Windows only, outside the library (no runtime dependency).
 * Each family is evaluated on the GPU (float) and by the same cores compiled
 * as C (double) on the same float inputs; the largest difference of the
 * position and derivatives, relative to 1 + |value|, is printed. The C cores
 * equal the runtime objects (test 144), so this closes the chain GPU = C core
 * = runtime.
 *
 * Build (f64 library):
 *   cl /O2 /MD /DQAWS_SCALAR_IS_FLOAT=0 /Isrc/qaws examples/eval_gpu_parity.c d3d11.lib d3dcompiler.lib dxgi.lib
 * Run from the repository root (it compiles examples/eval_parity.hlsl). Exit code: the families over 1e-4.
 */

#define COBJMACROS
#include <windows.h>
#include <d3d11.h>
#include <d3dcompiler.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "qaws_platform.h"
#include "qaws_core_types.h"
#include "core/qaws_arc_core.h"
#include "core/qaws_bspline_basis_core.h"
#include "core/qaws_bspline_eval_core.h"
#include "core/qaws_cubic_poly_core.h"
#include "core/qaws_decasteljau_core.h"
#include "core/qaws_horner_core.h"
#include "core/qaws_nurbs_eval_core.h"
#include "core/qaws_rational_bezier_core.h"
#include "core/qaws_surface_core.h"

#define EP_SAMPLES 256
#define EP_DATA 2048

typedef struct ep_params
{
	unsigned int family, count, degree, n, v_degree, v_count, off_knots, off_extra;
} ep_params;

typedef struct ep_ctx
{
	ID3D11Device* dev;
	ID3D11DeviceContext* ctx;
	ID3D11ComputeShader* cs[3];             /* curves, Bezier surface, B-spline surface */
} ep_ctx;

static HRESULT __stdcall ep_open(ID3DInclude* self, D3D_INCLUDE_TYPE type, LPCSTR name, LPCVOID parent, LPCVOID* data, UINT* bytes)
{
	static char const* const dirs[4] = { "examples/", "src/qaws/", "src/qaws/core/", "" };
	unsigned int d;
	(void)self; (void)type; (void)parent;
	for (d = 0; d < 4; d++)
	{
		char path[512];
		FILE* fp;
		long n;
		char* buf;
		char const* base = strncmp(name, "../src/qaws/", 12) == 0 ? name + 12 : strncmp(name, "../", 3) == 0 ? name + 3 : name;
		snprintf(path, sizeof(path), "%s%s", dirs[d], base);
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

static HRESULT __stdcall ep_close(ID3DInclude* self, LPCVOID data)
{
	(void)self;
	free((void*)data);
	return S_OK;
}

static ID3DIncludeVtbl ep_include_vtbl = { ep_open, ep_close };
static ID3DInclude ep_include = { &ep_include_vtbl };

static ID3D11Buffer* ep_buffer(ep_ctx* g, void const* data, unsigned int n, ID3D11ShaderResourceView** srv, ID3D11UnorderedAccessView** uav)
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

/* one dispatch: data and params in, count * stride floats out */
static void ep_run(ep_ctx* g, ep_params const* p, float const* data, float const* params, unsigned int nparams, unsigned int stride, float* out)
{
	ID3D11ShaderResourceView *sd = NULL, *sp = NULL;
	ID3D11UnorderedAccessView* uo = NULL;
	ID3D11Buffer *bd = ep_buffer(g, data, EP_DATA, &sd, NULL), *bp = ep_buffer(g, params, nparams, &sp, NULL);
	ID3D11Buffer* bo = ep_buffer(g, NULL, p->count * stride, NULL, &uo);
	ID3D11Buffer *cb = NULL, *st = NULL;
	ID3D11ShaderResourceView* srvs[2];
	D3D11_BUFFER_DESC cbd;
	D3D11_SUBRESOURCE_DATA csd;
	D3D11_MAPPED_SUBRESOURCE m;
	memset(&cbd, 0, sizeof(cbd));
	cbd.ByteWidth = sizeof(ep_params);
	cbd.Usage = D3D11_USAGE_DEFAULT;
	cbd.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
	csd.pSysMem = p; csd.SysMemPitch = 0; csd.SysMemSlicePitch = 0;
	ID3D11Device_CreateBuffer(g->dev, &cbd, &csd, &cb);
	srvs[0] = sd; srvs[1] = sp;
	ID3D11DeviceContext_CSSetShader(g->ctx, g->cs[p->family == 7 ? 1 : p->family == 8 ? 2 : 0], NULL, 0);
	ID3D11DeviceContext_CSSetShaderResources(g->ctx, 0, 2, srvs);
	ID3D11DeviceContext_CSSetUnorderedAccessViews(g->ctx, 0, 1, &uo, NULL);
	ID3D11DeviceContext_CSSetConstantBuffers(g->ctx, 0, 1, &cb);
	ID3D11DeviceContext_Dispatch(g->ctx, (p->count + 63) / 64, 1, 1);
	ID3D11Buffer_GetDesc(bo, &cbd);
	cbd.Usage = D3D11_USAGE_STAGING; cbd.BindFlags = 0; cbd.CPUAccessFlags = D3D11_CPU_ACCESS_READ; cbd.MiscFlags = 0;
	ID3D11Device_CreateBuffer(g->dev, &cbd, NULL, &st);
	ID3D11DeviceContext_CopyResource(g->ctx, (ID3D11Resource*)st, (ID3D11Resource*)bo);
	ID3D11DeviceContext_Map(g->ctx, (ID3D11Resource*)st, 0, D3D11_MAP_READ, 0, &m);
	memcpy(out, m.pData, 4 * p->count * stride);
	ID3D11DeviceContext_Unmap(g->ctx, (ID3D11Resource*)st, 0);
	ID3D11Buffer_Release(st); ID3D11Buffer_Release(cb); ID3D11Buffer_Release(bo); ID3D11Buffer_Release(bp); ID3D11Buffer_Release(bd);
	ID3D11UnorderedAccessView_Release(uo); ID3D11ShaderResourceView_Release(sp); ID3D11ShaderResourceView_Release(sd);
}

static double ep_err(float const* gpu, double const* cpu, unsigned int n)
{
	double worst = 0;
	unsigned int i;
	for (i = 0; i < n; i++)
	{
		double e = fabs((double)gpu[i] - cpu[i]) / (1 + fabs(cpu[i]));
		worst = e > worst ? e : worst;
	}
	return worst;
}

static void ep_put3(double* o, qaws_vec3 v)
{
	o[0] = v.x; o[1] = v.y; o[2] = v.z;
}

/* the C cores on the same inputs (data as float, computed in double) */
static void ep_cpu_curve(ep_params const* p, float const* data, double t, double* o)
{
	qaws_scalar cp[QAWS_CORE_MAX_POINTS * 4], cp3[QAWS_CORE_MAX_POINTS * 3], lk[QAWS_CORE_MAX_POINTS * 2], all[QAWS_CORE_MAX_POINTS * 2];
	qaws_scalar lw[QAWS_CORE_MAX_POINTS], wcp[QAWS_CORE_MAX_POINTS * 3];
	qaws_eval_3d e;
	int j, d = (int)p->degree;
	for (j = 0; j < QAWS_CORE_MAX_POINTS * 4; j++) cp[j] = j < (int)(p->n * 4) ? data[j] : 0;
	for (j = 0; j < QAWS_CORE_MAX_POINTS * 3; j++) cp3[j] = cp[j];
	if (p->family == 0)
		e = qaws_decasteljau_eval_3d(cp3, d, t, 7);
	else if (p->family == 1)
		e = qaws_rational_bezier_eval_3d(cp, d, t, 7);
	else if (p->family == 2 || p->family == 3)
	{
		int span;
		for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++) all[j] = j < (int)(p->n + p->degree + 1) ? data[p->off_knots + j] : 0;
		span = qaws_find_span(all, d, (int)p->n, t);
		for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++) lk[j] = j < 2 * (d + 1) ? all[span - d + j] : 0;
		for (j = 0; j < QAWS_CORE_MAX_POINTS * 3; j++) wcp[j] = j < 3 * (d + 1) ? data[(span - d) * 3 + j] : 0;
		for (j = 0; j < QAWS_CORE_MAX_POINTS; j++) lw[j] = p->family == 3 && j <= d ? data[p->off_extra + span - d + j] : 1;
		e = p->family == 2 ? qaws_bspline_eval_3d(wcp, lk, d, d, t, 7) : qaws_nurbs_eval_3d(wcp, lk, lw, d, d, t, 7);
	}
	else if (p->family == 4)
	{
		qaws_vec3 a, b, c, dd;
		a.x = cp[0]; a.y = cp[1]; a.z = cp[2]; b.x = cp[3]; b.y = cp[4]; b.z = cp[5];
		c.x = cp[6]; c.y = cp[7]; c.z = cp[8]; dd.x = cp[9]; dd.y = cp[10]; dd.z = cp[11];
		e = qaws_cubic_eval_3d(a, b, c, dd, t, 7);
	}
	else if (p->family == 5)
		e = qaws_horner_eval_3d(cp3, d, t, 7);
	else
	{
		qaws_vec3 center, au, av;
		center.x = cp[0]; center.y = cp[1]; center.z = cp[2]; au.x = cp[3]; au.y = cp[4]; au.z = cp[5]; av.x = cp[6]; av.y = cp[7]; av.z = cp[8];
		e = qaws_arc_eval_full_3d(center, au, av, cp[9], cp[10], t, 7);
	}
	ep_put3(o, e.position); ep_put3(o + 3, e.d1); ep_put3(o + 6, e.d2);
}

static void ep_cpu_surface(ep_params const* p, float const* data, double u, double v, double* o)
{
	static qaws_scalar grid[QAWS_CORE_MAX_POINTS * QAWS_CORE_MAX_POINTS * 3];
	qaws_scalar uk[QAWS_CORE_MAX_POINTS * 2], vk[QAWS_CORE_MAX_POINTS * 2];
	qaws_surface_eval_3d s;
	int j;
	for (j = 0; j < (int)(p->n * p->v_count * 3); j++) grid[j] = data[j];
	if (p->family == 7)
		s = qaws_surface_bezier_eval_3d(grid, (int)p->degree, (int)p->v_degree, (int)p->v_count, u, v, 127);
	else
	{
		for (j = 0; j < QAWS_CORE_MAX_POINTS * 2; j++)
		{
			uk[j] = j < (int)(p->n + p->degree + 1) ? data[p->off_knots + j] : 0;
			vk[j] = j < (int)(p->v_count + p->v_degree + 1) ? data[p->off_extra + j] : 0;
		}
		s = qaws_surface_bspline_eval_3d(grid, uk, vk, (int)p->degree, (int)p->v_degree, (int)p->n, (int)p->v_count,
			qaws_find_span(uk, (int)p->degree, (int)p->n, u), qaws_find_span(vk, (int)p->v_degree, (int)p->v_count, v), u, v, 127);
	}
	ep_put3(o, s.position); ep_put3(o + 3, s.du); ep_put3(o + 6, s.dv); ep_put3(o + 9, s.normal);
	ep_put3(o + 12, s.duu); ep_put3(o + 15, s.duv); ep_put3(o + 18, s.dvv);
}

int main(void)
{
	static char const* const names[9] = { "Bezier (degree 5)", "rational Bezier (degree 4)", "B-spline (cubic, 9 points, a double knot)",
		"NURBS (cubic, 9 points)", "cubic", "polynomial (degree 4)", "arc", "Bezier surface (3 x 3)", "B-spline surface (3 x 2, 6 x 5 points)" };
	static float data[EP_DATA], params[2 * EP_SAMPLES], gpu[21 * EP_SAMPLES];
	static double cpu[21 * EP_SAMPLES];
	ep_ctx g;
	D3D_FEATURE_LEVEL fl;
	static char const* const entries[3] = { "CSCurve", "CSBezierSurface", "CSBsplineSurface" };
	unsigned int f, i, bad = 0;
	memset(&g, 0, sizeof(g));
	if (FAILED(D3D11CreateDevice(NULL, D3D_DRIVER_TYPE_HARDWARE, NULL, 0, NULL, 0, D3D11_SDK_VERSION, &g.dev, &fl, &g.ctx)))
	{
		printf("no Direct3D 11 device\n");
		return 1;
	}
	for (i = 0; i < 3; i++)
	{
		ID3DBlob *blob = NULL, *err = NULL;
		if (FAILED(D3DCompileFromFile(L"examples/eval_parity.hlsl", NULL, &ep_include, entries[i], "cs_5_0", D3DCOMPILE_OPTIMIZATION_LEVEL3, 0, &blob, &err)))
		{
			printf("compile %s failed: %s\n", entries[i], err ? (char const*)ID3D10Blob_GetBufferPointer(err) : "?");
			return 1;
		}
		ID3D11Device_CreateComputeShader(g.dev, ID3D10Blob_GetBufferPointer(blob), ID3D10Blob_GetBufferSize(blob), NULL, &g.cs[i]);
		ID3D10Blob_Release(blob);
	}
	for (f = 0; f < 9; f++)
	{
		ep_params p;
		unsigned int stride = f >= 7 ? 21 : 9;
		double worst;
		memset(&p, 0, sizeof(p));
		memset(data, 0, sizeof(data));
		p.family = f;
		p.count = EP_SAMPLES;
		for (i = 0; i < 300; i++)
			data[i] = (float)(0.9 * sin(1.37 * i + 0.4 * f) + 0.1 * i / 30.0);
		switch (f)
		{
		case 0: p.degree = 5; p.n = 6; break;
		case 1: p.degree = 4; p.n = 5; for (i = 0; i < 5; i++) data[4 * i + 3] = (float)(0.6 + 0.2 * i); break;
		case 2: case 3:
		{
			static float const k[13] = { 0, 0, 0, 0, 0.5f, 1.25f, 1.25f, 2, 3.5f, 4, 4, 4, 4 };
			p.degree = 3; p.n = 9; p.off_knots = 512; p.off_extra = 600;
			memcpy(data + 512, k, sizeof(k));
			for (i = 0; i < 9; i++) data[600 + i] = (float)(0.5 + 0.25 * (i % 4));
			break;
		}
		case 5: p.degree = 4; p.n = 5; break;
		case 4: p.n = 3; break;   /* (12 values: a, b, c, d) */
		case 6: p.n = 3; data[9] = 0.3f; data[10] = 4.0f; break;
		case 7: p.degree = 3; p.n = 4; p.v_degree = 3; p.v_count = 4; break;
		case 8:
		{
			static float const uk[10] = { 0, 0, 0, 0, 0.4f, 0.4f, 1, 1, 1, 1 }, vk[8] = { 0, 0, 0, 0.3f, 0.75f, 1, 1, 1 };
			p.degree = 3; p.n = 6; p.v_degree = 2; p.v_count = 5; p.off_knots = 512; p.off_extra = 600;
			memcpy(data + 512, uk, sizeof(uk));
			memcpy(data + 600, vk, sizeof(vk));
			break;
		}
		default: break;
		}
		for (i = 0; i < EP_SAMPLES; i++)
		{
			double t = (f == 2 || f == 3 ? 4.0 : 1.0) * i / (EP_SAMPLES - 1);
			if (f >= 7)
			{
				params[2 * i] = (float)((i % 16) / 15.0);
				params[2 * i + 1] = (float)((i / 16) / 15.0);
				ep_cpu_surface(&p, data, params[2 * i], params[2 * i + 1], cpu + stride * i);
			}
			else
			{
				params[i] = (float)t;
				ep_cpu_curve(&p, data, params[i], cpu + stride * i);
			}
		}
		ep_run(&g, &p, data, params, f >= 7 ? 2 * EP_SAMPLES : EP_SAMPLES, stride, gpu);
		worst = ep_err(gpu, cpu, stride * EP_SAMPLES);
		{
			/* a family evaluated on nothing compares zeros: refused */
			double size = 0;
			for (i = 0; i < stride * EP_SAMPLES; i++) size = fabs(cpu[i]) > size ? fabs(cpu[i]) : size;
			if (!(size > 0)) worst = HUGE_VAL;
		}
		bad += !(worst <= 1e-4);
		printf("%-45s GPU against the C core: %.2g%s\n", names[f], worst, worst <= 1e-4 ? "" : "  (over 1e-4, or no data)");
	}
	return (int)bad;
}
