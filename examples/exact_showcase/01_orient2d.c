/* ================================================================== */
/*  1. Orientation near a line: naive f64 against certified           */
/* ================================================================== */

#define OR_N 256

/* One map, run-length encoded per row: value -> color. */
static void or_map(svg* s, double x0, double y0, double cell, int N, signed char const* v, char const* const* colors, char const* title)
{
	int i, j;
	fprintf(s->f, "<text x=\"%.1f\" y=\"%.1f\" font-size=\"14\" font-weight=\"600\" fill=\"#24292f\">%s</text>\n", x0, y0 - 10, title);
	for (j = 0; j < N; j++)
	{
		int start = 0;
		for (i = 1; i <= N; i++)
			if (i == N || v[j * N + i] != v[j * N + start])
			{
				/* row j is drawn top-down from the largest y */
				fprintf(s->f, "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\" shape-rendering=\"crispEdges\"/>\n",
					x0 + start * cell, y0 + (N - 1 - j) * cell, (i - start) * cell + 0.5, cell + 0.5, colors[v[j * N + start] + 1]);
				start = i;
			}
	}
	fprintf(s->f, "<rect x=\"%.1f\" y=\"%.1f\" width=\"%.1f\" height=\"%.1f\" fill=\"none\" stroke=\"#8c959f\"/>\n", x0, y0, N * cell,
		N * cell);
}

static void demo_orient2d(void)
{
	static signed char naive[OR_N * OR_N], exact[OR_N * OR_N], wrong[OR_N * OR_N], path[OR_N * OR_N];
	static char const* const sign_colors[3] = { "#f2a900", "#24292f", "#0969da" };   /* -1, 0, +1 */
	static char const* const wrong_colors[3] = { "#f6f8fa", "#f6f8fa", "#cf222e" };
	static char const* const path_colors[3] = { "#f6f8fa", "#d0d7de", "#8250df" };
	double b[2] = { 12, 12 }, c[2] = { 24, 24 }, u = ldexp(1.0, -53);
	int i, j, n_wrong = 0, n_exact = 0, n_zero = 0;
	char buf[300];
	svg s;
	for (j = 0; j < OR_N; j++)
		for (i = 0; i < OR_N; i++)
		{
			double a[2], det;
			qaws_exact_sign sg;
			qaws_exact_path p;
			int k = j * OR_N + i;
			a[0] = 0.5 + i * u;
			a[1] = 0.5 + j * u;
			det = (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0]);
			naive[k] = (signed char)(det > 0 ? 1 : (det < 0 ? -1 : 0));
			qaws_exact_orient2d(a, b, c, &sg, &p);
			exact[k] = (signed char)sg;
			wrong[k] = (signed char)(naive[k] != exact[k]);
			path[k] = (signed char)(p == QAWS_EXACT_PATH_EXACT);
			n_wrong += wrong[k];
			n_exact += path[k];
			n_zero += exact[k] == 0;
		}
	svg_open(&s, "showcase/exact1_orient2d.svg", 1270, 470, "Certified orientation: points one ulp apart near a line",
		"a = (0.5 + i u, 0.5 + j u), u = 2^-53, 0 &lt;= i, j &lt; 256; sign of orient2d(a, (12, 12), (24, 24)). Blue: left, orange: right, black: on the line.");
	or_map(&s, 20, 100, 1.15, OR_N, naive, sign_colors, "naive f64 determinant");
	or_map(&s, 330, 100, 1.15, OR_N, exact, sign_colors, "qaws_exact_orient2d (certified)");
	or_map(&s, 640, 100, 1.15, OR_N, wrong, wrong_colors, "cells where f64 is wrong");
	or_map(&s, 950, 100, 1.15, OR_N, path, path_colors, "decided by the exact fallback");
	sprintf(buf, "%d of %d cells: f64 sign wrong on %d (%.1f %%); %d points exactly on the line; the f64 filter proved %d signs, the exact path decided %d",
		OR_N * OR_N, OR_N * OR_N, n_wrong, 100.0 * n_wrong / (OR_N * OR_N), n_zero, OR_N * OR_N - n_exact, n_exact);
	svg_text(&s, 20, 440, 13, "#57606a", "start", buf);
	svg_close(&s);
	printf("exact1_orient2d: f64 wrong on %d of %d cells, %d on the line, exact fallback %d\n", n_wrong, OR_N * OR_N, n_zero, n_exact);
}
