/*
 * qaws_api_tool - Qaws's API checks (docs/parameter_convention_plan.md)
 *
 *   qaws_api_tool survey [-v] [-order] [-allow file] [-types header... --] header...
 *       checks every function declaration against the parameter
 *       convention shared with Alwan, Suwar and Abaad:
 *         R1 the context last;
 *         R2 a stride right after its buffer;
 *         R3 outputs before inputs (plain outputs, in-out buffers, then
 *            the object handle, input buffers, sizes, options);
 *         R4 sizes after the buffers.
 *       Exit 1 on a violation not in the allow-list, or an allow-list
 *       entry that no longer applies (the list only shrinks).
 *
 *   qaws_api_tool names [-check file] -once header... -- header...
 *       writes qaws_names.h: every precision-bearing name (functions,
 *       types, struct tags declared by the headers after "--", less those
 *       of the once-only headers and the enums) mapped to QAWS_SFX(name).
 *
 * The headers are tokenized (comments, strings and preprocessor lines
 * dropped) and read statement by statement at file scope.
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  Tokens                                                            */
/* ------------------------------------------------------------------ */

typedef struct tok
{
	char s[96];
} tok;

typedef struct toks
{
	tok* t;
	size_t n, cap;
} toks;

static void push(toks* v, char const* s, size_t len)
{
	if (v->n == v->cap)
	{
		v->cap = v->cap ? 2 * v->cap : 4096;
		v->t = (tok*)realloc(v->t, sizeof(tok) * v->cap);
		if (!v->t) { fprintf(stderr, "out of memory\n"); exit(2); }
	}
	if (len >= sizeof(v->t[0].s)) len = sizeof(v->t[0].s) - 1;
	memcpy(v->t[v->n].s, s, len);
	v->t[v->n].s[len] = 0;
	v->n++;
}

static int is_id(char c) { return isalnum((unsigned char)c) || c == '_'; }

/* the file's tokens: identifiers, numbers, single punctuation */
static int tokenize(char const* path, toks* v)
{
	FILE* f = fopen(path, "rb");
	char* s;
	long n, i;
	int bol = 1;
	if (!f) { fprintf(stderr, "cannot read %s\n", path); return 0; }
	fseek(f, 0, SEEK_END); n = ftell(f); fseek(f, 0, SEEK_SET);
	s = (char*)malloc((size_t)n + 1);
	if (!s || fread(s, 1, (size_t)n, f) != (size_t)n) { fclose(f); free(s); return 0; }
	s[n] = 0;
	fclose(f);
	for (i = 0; i < n;)
	{
		char c = s[i];
		if (c == '\n') { bol = 1; i++; continue; }
		if (isspace((unsigned char)c)) { i++; continue; }
		if (c == '/' && s[i + 1] == '*') { char* e = strstr(s + i + 2, "*/"); i = e ? (long)(e - s) + 2 : n; continue; }
		if (c == '/' && s[i + 1] == '/') { while (i < n && s[i] != '\n') i++; continue; }
		if (bol && c == '#')
		{
			/* a preprocessor line, with its continuations */
			while (i < n && s[i] != '\n') { if (s[i] == '\\' && s[i + 1] == '\n') i++; else if (s[i] == '\\' && s[i + 1] == '\r') i += 2; i++; }
			continue;
		}
		bol = 0;
		if (c == '"' || c == '\'')
		{
			long j = i + 1;
			while (j < n && s[j] != c) { if (s[j] == '\\') j++; j++; }
			push(v, "\"\"", 2);
			i = j + 1;
			continue;
		}
		if (is_id(c))
		{
			long j = i;
			while (j < n && is_id(s[j])) j++;
			push(v, s + i, (size_t)(j - i));
			i = j;
			continue;
		}
		push(v, s + i, 1);
		i++;
	}
	free(s);
	return 1;
}

/* ------------------------------------------------------------------ */
/*  Statements at file scope                                          */
/* ------------------------------------------------------------------ */

/* calls fn(user, v, first, end, has_body) per statement: a declaration up
   to ';', or a definition up to its body's '}' (the body excluded) */
static void statements(toks const* v, void (*fn)(void*, toks const*, size_t, size_t, int), void* user)
{
	size_t i = 0, start = 0;
	int paren = 0, brace = 0;
	while (i < v->n)
	{
		char const* s = v->t[i].s;
		if (!strcmp(s, "(")) paren++;
		else if (!strcmp(s, ")")) paren--;
		else if (!strcmp(s, ";") && paren == 0 && brace == 0)
		{
			fn(user, v, start, i, 0);
			start = i + 1;
		}
		else if (!strcmp(s, "{") && paren == 0)
		{
			/* a body: a function's ends the statement (skipped); a struct's
			   or an enum's is part of a declaration going on to ';' */
			if (brace == 0 && i > 0 && !strcmp(v->t[i - 1].s, ")"))
			{
				size_t j = i, depth = 0;
				for (; j < v->n; j++)
				{
					if (!strcmp(v->t[j].s, "{")) depth++;
					else if (!strcmp(v->t[j].s, "}") && --depth == 0) break;
				}
				fn(user, v, start, i, 1);
				start = j + 1;
				i = j + 1;
				continue;
			}
			brace++;
		}
		else if (!strcmp(s, "}") && paren == 0)
		{
			if (brace > 0)
				brace--;
			else
				start = i + 1;   /* (the closing brace of extern "C") */
		}
		else if (!strcmp(s, "extern") && i + 1 < v->n && !strcmp(v->t[i + 1].s, "\"\"") && brace == 0)
		{
			/* extern "C" { : its brace is not a body */
			i += 2;
			if (i < v->n && !strcmp(v->t[i].s, "{")) i++;
			start = i;
			continue;
		}
		i++;
	}
}

/* ------------------------------------------------------------------ */
/*  The survey                                                        */
/* ------------------------------------------------------------------ */

enum kind { K_OUT, K_HANDLE, K_IN, K_KNOB, K_SIZE, K_VALUE, K_STRIDE, K_CTX };

typedef struct survey
{
	int verbose;
	int order;                       /* print each violator's spec (the order) instead */
	char const* file;
	char** allow;                    /* "name rule" lines */
	unsigned char* allow_used;
	size_t nallow;
	unsigned int fns, ok, rule[5], unexpected;
} survey;

static int scalarish(char const* s)
{
	static char const* const list[] = { "qaws_scalar", "qaws_f32", "qaws_f64", "double", "float", "int", "unsigned", "size_t", "char",
		"void", "long", "short", "uint8_t", "uint16_t", "uint32_t", "uint64_t", "int8_t", "int16_t", "int32_t", "int64_t",
		"qaws_vec2", "qaws_vec3", "qaws_vec4", "qaws_i64", "qaws_u64", NULL };
	int k;
	for (k = 0; list[k]; k++) if (!strcmp(s, list[k])) return 1;
	return 0;
}

static int contains(char const* s, char const* part) { return strstr(s, part) != NULL; }

/* the structs a public header gives a body: data, passed by pointer as
   buffers; every other struct is an opaque handle (collected first) */
static char** g_data;
static size_t g_ndata;

static int is_data(char const* s)
{
	size_t k;
	for (k = 0; k < g_ndata; k++) if (!strcmp(g_data[k], s)) return 1;
	return 0;
}

static void add_data(char const* s)
{
	size_t l = strlen(s);
	if (is_data(s)) return;
	g_data = (char**)realloc(g_data, sizeof(char*) * (g_ndata + 1));
	g_data[g_ndata] = (char*)malloc(l + 1);
	memcpy(g_data[g_ndata], s, l + 1);
	g_ndata++;
}

static void data_stmt(void* user, toks const* v, size_t a, size_t b, int body)
{
	size_t i;
	int has_body = 0, is_struct = 0;
	(void)user; (void)body;
	for (i = a; i < b; i++)
	{
		if (!strcmp(v->t[i].s, "{")) has_body = 1;
		if (!strcmp(v->t[i].s, "struct") || !strcmp(v->t[i].s, "union")) is_struct = 1;
	}
	if (!has_body || !is_struct) return;
	for (i = a; i + 1 < b; i++)
		if ((!strcmp(v->t[i].s, "struct") || !strcmp(v->t[i].s, "union")) && is_id(v->t[i + 1].s[0])) { add_data(v->t[i + 1].s); break; }
	if (is_id(v->t[b - 1].s[0])) add_data(v->t[b - 1].s);   /* (typedef struct { ... } name;) */
}

/* a parameter's facts: pointer depth, whether the pointee is const (a
   const before the first '*'), its type and name */
typedef struct param
{
	int ptr, cnst, scal, data, out, inout;
	char const* name;
	char const* type;
} param;

static param analyze(toks const* v, size_t a, size_t b)
{
	param p;
	size_t i;
	memset(&p, 0, sizeof(p));
	p.name = ""; p.type = "";
	int outer_const = 0;   /* (a const right before the last '*': the pointed-to pointers are const) */
	for (i = a; i < b; i++)
	{
		char const* s = v->t[i].s;
		if (!strcmp(s, "*") || !strcmp(s, "[")) { p.ptr++; outer_const = i > a && !strcmp(v->t[i - 1].s, "const") && p.ptr >= 2; }
		else if (!strcmp(s, "const")) { if (p.ptr == 0) p.cnst = 1; }
		else if (!strcmp(s, "QAWS_OUT")) p.out = 1;
		else if (!strcmp(s, "QAWS_INOUT")) p.inout = 1;
		else if (is_id(s[0]) && !isdigit((unsigned char)s[0]))
		{
			if (!strcmp(p.type, "") || !strcmp(p.type, "struct") || !strcmp(p.type, "unsigned") || !strcmp(p.type, "signed")) p.type = s;
			else p.name = s;
			if (scalarish(s)) p.scal = 1;
		}
	}
	if (!strcmp(p.name, "")) p.name = p.type;
	p.data = p.scal || is_data(p.type);
	/* T const** (a pointer handed back) writes: an output, as T* const*;
	   only T const* const* reads at every level */
	if (p.ptr >= 2 && p.cnst && !outer_const) p.cnst = 0;
	return p;
}

static int is_ctx(param const* p)
{
	return contains(p->type, "context") || contains(p->type, "_ctx") || contains(p->type, "executor") || contains(p->type, "allocator") ||
		!strcmp(p->name, "ctx");
}

static int is_size(param const* p)
{
	return !p->ptr && (contains(p->name, "count") || !strcmp(p->name, "n") || contains(p->name, "width") || contains(p->name, "height") ||
		!strcmp(p->name, "dim") || !strcmp(p->name, "dimension") || !strcmp(p->name, "size") || contains(p->name, "capacity"));
}

static int is_user_data(param const* p)
{
	return !strcmp(p->type, "void") && (!strcmp(p->name, "user") || !strcmp(p->name, "user_data") || contains(p->name, "_user") ||
		!strcmp(p->name, "impl") || !strcmp(p->name, "ptr"));
}

static int is_desc(param const* p)
{
	size_t l = strlen(p->name);
	if (l > 1 && p->name[l - 1] == 's' && p->ptr)
		return 0;   /* (descs, faces: an array of them, a buffer its count follows) */
	return contains(p->type, "desc") || contains(p->type, "options") || contains(p->type, "params") || contains(p->type, "settings") ||
		contains(p->type, "config") || contains(p->type, "targets") || contains(p->type, "spec");
}

/* a parameter tokens [a, b): its kind */
static enum kind classify(toks const* v, size_t a, size_t b)
{
	param p = analyze(v, a, b);
	if (is_ctx(&p)) return K_CTX;
	if (contains(p.name, "stride")) return K_STRIDE;
	if (p.out || p.inout) return K_OUT;
	if (p.ptr >= 2 && !p.cnst) return K_OUT;
	if (p.ptr && !p.cnst)
	{
		if (is_user_data(&p)) return K_VALUE;
		if (!strncmp(p.name, "out", 3) || p.data) return K_OUT;
		return K_HANDLE;
	}
	if (p.ptr) return p.data && !is_desc(&p) ? K_IN : K_KNOB;
	if (is_size(&p)) return K_SIZE;
	return K_VALUE;
}

/* a parameter's place in the convention's order: 0 plain outputs, 1
   reports, 2 in-out buffers, 3 object handles, 4 input buffers, 5 sizes,
   6 descriptors and options, 7 scalars, enums, callbacks and their data,
   8 the context (a stride takes its buffer's place, right after it) */
static int rank_of(toks const* v, size_t a, size_t b)
{
	param p = analyze(v, a, b);
	if (is_ctx(&p)) return 8;
	if (p.out) return 0;
	if (p.inout) return 2;
	if (p.ptr >= 2 && !p.cnst) return 0;
	if (p.ptr && !p.cnst)
	{
		if (is_user_data(&p)) return 7;
		if (!strncmp(p.name, "out", 3)) return 0;
		if (contains(p.type, "report") || contains(p.type, "result") || contains(p.type, "stats") || contains(p.type, "summary")) return 1;
		return p.data ? 2 : 3;
	}
	if (p.ptr)
	{
		if (is_desc(&p)) return 6;
		return p.data ? 4 : 3;
	}
	if (is_size(&p)) return 5;
	return 7;
}

static int allowed(survey* S, char const* fn, char const* rule)
{
	size_t k;
	char key[256];
	sprintf(key, "%s %s", fn, rule);
	for (k = 0; k < S->nallow; k++)
		if (!strcmp(S->allow[k], key)) { S->allow_used[k] = 1; return 1; }
	return 0;
}

static void survey_stmt(void* user, toks const* v, size_t a, size_t b, int body)
{
	survey* S = (survey*)user;
	size_t i, open = 0, close, p, q;
	enum kind kinds[64];
	int ranks[64];
	int nk = 0, bad[5] = { 0, 0, 0, 0, 0 };
	char const* fn = NULL;
	(void)body;
	if (b <= a) return;
	for (i = a; i < b; i++)
		if (!strcmp(v->t[i].s, "typedef") || !strcmp(v->t[i].s, "static")) return;   /* (inline helpers and types: not API) */
	/* name ( params ) at the top level, the name a qaws_ identifier */
	for (i = a + 1; i < b; i++)
		if (!strcmp(v->t[i].s, "(") && !strncmp(v->t[i - 1].s, "qaws_", 5)) { open = i; fn = v->t[i - 1].s; break; }
	if (!fn) return;
	{
		int depth = 0;
		for (close = open; close < b; close++)
		{
			if (!strcmp(v->t[close].s, "(")) depth++;
			else if (!strcmp(v->t[close].s, ")") && --depth == 0) break;
		}
	}
	if (close >= b) return;
	if (close == open + 1 || (close == open + 2 && !strcmp(v->t[open + 1].s, "void"))) { S->fns++; S->ok++; return; }
	for (p = open + 1; p < close && nk < 64;)
	{
		int depth = 0;
		for (q = p; q < close; q++)
		{
			if (!strcmp(v->t[q].s, "(")) depth++;
			else if (!strcmp(v->t[q].s, ")")) depth--;
			else if (!strcmp(v->t[q].s, ",") && depth == 0) break;
		}
		ranks[nk] = rank_of(v, p, q);
		kinds[nk++] = classify(v, p, q);
		p = q + 1;
	}
	for (i = 1; i < (size_t)nk; i++)
		if (kinds[i] == K_STRIDE) ranks[i] = ranks[i - 1];   /* (a stride stays with its buffer) */
	S->fns++;
	/* R1 */
	for (i = 0; i < (size_t)nk; i++) if (kinds[i] == K_CTX && i + 1 != (size_t)nk) bad[0] = 1;
	/* R2 */
	for (i = 0; i < (size_t)nk; i++) if (kinds[i] == K_STRIDE && (i == 0 || (kinds[i - 1] != K_IN && kinds[i - 1] != K_OUT))) bad[1] = 1;
	/* R3: no output after a non-output */
	{
		int seen = 0;
		for (i = 0; i < (size_t)nk; i++)
		{
			if (kinds[i] == K_OUT && seen) bad[2] = 1;
			if (kinds[i] != K_OUT && kinds[i] != K_CTX && kinds[i] != K_STRIDE) seen = 1;
		}
	}
	/* R4: no size before the last buffer */
	{
		int last = -1;
		for (i = 0; i < (size_t)nk; i++) if (kinds[i] == K_IN || kinds[i] == K_OUT) last = (int)i;
		for (i = 0; (int)i < last; i++) if (kinds[i] == K_SIZE) bad[3] = 1;
	}
	{
		static char const* const names[4] = { "R1", "R2", "R3", "R4" };
		int any = 0;
		for (i = 0; i < 4; i++)
			if (bad[i])
			{
				any = 1;
				S->rule[i]++;
				if (!allowed(S, fn, names[i]))
				{
					S->unexpected++;
					printf("%s: %s: %s\n", S->file, fn, names[i]);
				}
				else if (S->verbose)
					printf("%s: %s: %s (allowed)\n", S->file, fn, names[i]);
			}
		if (!any) S->ok++;
		if (any && S->order)
		{
			/* the convention's order: a stable sort by rank, as old indices */
			int idx[64], k, m;
			for (k = 0; k < nk; k++) idx[k] = k;
			for (k = 1; k < nk; k++)
				for (m = k; m > 0 && ranks[idx[m - 1]] > ranks[idx[m]]; m--) { int t = idx[m]; idx[m] = idx[m - 1]; idx[m - 1] = t; }
			printf("%s", fn);
			for (k = 0; k < nk; k++) printf(" %d", idx[k]);
			printf("\n");
		}
	}
}

static int load_allow(char const* path, survey* S)
{
	FILE* f = fopen(path, "r");
	char line[512];
	if (!f) { fprintf(stderr, "cannot read %s\n", path); return 0; }
	while (fgets(line, sizeof(line), f))
	{
		size_t l = strlen(line);
		while (l && (line[l - 1] == '\n' || line[l - 1] == '\r' || line[l - 1] == ' ')) line[--l] = 0;
		if (!l || line[0] == '#') continue;
		S->allow = (char**)realloc(S->allow, sizeof(char*) * (S->nallow + 1));
		S->allow[S->nallow] = (char*)malloc(l + 1);
		memcpy(S->allow[S->nallow], line, l + 1);
		S->nallow++;
	}
	fclose(f);
	S->allow_used = (unsigned char*)calloc(S->nallow + 1, 1);
	return 1;
}

static int run_survey(int argc, char** argv)
{
	survey S;
	int i, stale = 0;
	size_t k;
	memset(&S, 0, sizeof(S));
	/* first the data structs: the bodies the public headers give (the
	   surveyed ones and those after -types, read for their types only) */
	for (i = 0; i < argc; i++)
	{
		toks v = { NULL, 0, 0 };
		if (argv[i][0] == '-') { if (!strcmp(argv[i], "-allow")) i++; continue; }
		if (contains(argv[i], "internal")) continue;
		if (!tokenize(argv[i], &v)) return 2;
		statements(&v, data_stmt, NULL);
		free(v.t);
	}
	for (i = 0; i < argc; i++)
	{
		if (!strcmp(argv[i], "-types")) { while (i + 1 < argc && strcmp(argv[i + 1], "--")) i++; i++; continue; }
		if (!strcmp(argv[i], "-v")) { S.verbose = 1; continue; }
		if (!strcmp(argv[i], "-order")) { S.order = 1; continue; }
		if (!strcmp(argv[i], "-allow") && i + 1 < argc) { if (!load_allow(argv[++i], &S)) return 2; continue; }
		{
			toks v = { NULL, 0, 0 };
			char const* base = strrchr(argv[i], '/');
			char const* base2 = strrchr(argv[i], '\\');
			if (base2 && (!base || base2 > base)) base = base2;
			S.file = base ? base + 1 : argv[i];
			if (!tokenize(argv[i], &v)) return 2;
			statements(&v, survey_stmt, &S);
			free(v.t);
		}
	}
	for (k = 0; k < S.nallow; k++)
		if (!S.allow_used[k]) { printf("stale allow-list entry: %s\n", S.allow[k]); stale++; }
	printf("%u functions, %u compliant; R1 %u, R2 %u, R3 %u, R4 %u; %u not allowed, %d stale allowances\n", S.fns, S.ok, S.rule[0], S.rule[1],
		S.rule[2], S.rule[3], S.unexpected, stale);
	return S.unexpected || stale ? 1 : 0;
}

/* ------------------------------------------------------------------ */
/*  The names                                                         */
/* ------------------------------------------------------------------ */

typedef struct names
{
	char** s;
	size_t n;
	int collect_once;     /* collecting the once-only headers' names (excluded) */
	char** once;
	size_t nonce;
} names;

static int has(char** list, size_t n, char const* s)
{
	size_t k;
	for (k = 0; k < n; k++) if (!strcmp(list[k], s)) return 1;
	return 0;
}

static void add(char*** list, size_t* n, char const* s)
{
	size_t l = strlen(s);
	if (strncmp(s, "qaws_", 5) || has(*list, *n, s)) return;
	*list = (char**)realloc(*list, sizeof(char*) * (*n + 1));
	(*list)[*n] = (char*)malloc(l + 1);
	memcpy((*list)[*n], s, l + 1);
	(*n)++;
}

static void names_stmt(void* user, toks const* v, size_t a, size_t b, int body)
{
	names* N = (names*)user;
	char*** list = N->collect_once ? &N->once : &N->s;
	size_t* n = N->collect_once ? &N->nonce : &N->n;
	size_t i;
	int is_enum = 0, is_typedef = 0;
	(void)body;
	for (i = a; i < b; i++)
	{
		if (!strcmp(v->t[i].s, "enum")) is_enum = 1;
		if (!strcmp(v->t[i].s, "typedef")) is_typedef = 1;
	}
	if (is_enum)
	{
		/* enums are precision-free: their names go to the once list */
		for (i = a; i < b; i++)
			if (!strncmp(v->t[i].s, "qaws_", 5) && (i == a || strcmp(v->t[i - 1].s, "=")))
				add(&N->once, &N->nonce, v->t[i].s);
		return;
	}
	/* struct tags (struct X), the typedef'd name, function names (X before
	   a top-level '('), a function pointer's name ((*X)) */
	for (i = a; i < b; i++)
	{
		char const* s = v->t[i].s;
		if (strncmp(s, "qaws_", 5)) continue;
		if (i > a && !strcmp(v->t[i - 1].s, "struct")) { add(list, n, s); continue; }
		if (i + 1 < b && !strcmp(v->t[i + 1].s, "(") && (i == a || strcmp(v->t[i - 1].s, "(")))
		{
			/* a function: only at depth 0 */
			size_t j;
			int depth = 0;
			for (j = a; j < i; j++)
			{
				if (!strcmp(v->t[j].s, "(") || !strcmp(v->t[j].s, "{")) depth++;
				else if (!strcmp(v->t[j].s, ")") || !strcmp(v->t[j].s, "}")) depth--;
			}
			if (depth == 0) add(list, n, s);
			continue;
		}
		if (i > a + 1 && !strcmp(v->t[i - 1].s, "*") && !strcmp(v->t[i - 2].s, "(") && is_typedef) { add(list, n, s); continue; }
		if (is_typedef && i + 1 == b) add(list, n, s);
		if (is_typedef && i + 1 < b && (!strcmp(v->t[i + 1].s, "[")) && !strcmp(v->t[b - 1].s, "]")) add(list, n, s);
	}
}

static int cmp(void const* p, void const* q) { return strcmp(*(char* const*)p, *(char* const*)q); }

/* the map as text (malloc'd) */
static char* names_text(names const* N)
{
	size_t k, cap = 4096, len = 0;
	char* out = (char*)malloc(cap);
	char line[400];
	static char const head[] = "/* generated by tools/qaws_api_tool names: every precision-bearing name of\n"
		"   Qaws mapped to its form for the current precision (QAWS_SFX, set by\n"
		"   qaws_precision.h to _f32 or _f64). Do not edit: regenerate. */\n\n"
		"#ifndef QAWS_NAMES_H\n#define QAWS_NAMES_H\n\n";
	static char const tail[] = "\n#endif /* QAWS_NAMES_H */\n";
	for (k = 0; k <= N->n + 1 && out; k++)
	{
		size_t l;
		if (k == 0) { strcpy(line, head); }
		else if (k == N->n + 1) { strcpy(line, tail); }
		else
		{
			char const* s = N->s[k - 1];
			if (has(N->once, N->nonce, s) || !strcmp(s, "qaws_scalar") || !strcmp(s, "qaws_f32") || !strcmp(s, "qaws_f64"))
				continue;
			sprintf(line, "#define %s QAWS_SFX(%s)\n", s, s);
		}
		l = strlen(line);
		if (len + l + 1 > cap)
		{
			cap = 2 * (len + l + 1);
			out = (char*)realloc(out, cap);
			if (!out) break;
		}
		memcpy(out + len, line, l + 1);
		len += l;
	}
	return out;
}

/* the file's contents with CR dropped (malloc'd; NULL if unreadable) */
static char* read_text(char const* path)
{
	FILE* f = fopen(path, "rb");
	char* s;
	long n, i, j;
	if (!f) return NULL;
	fseek(f, 0, SEEK_END); n = ftell(f); fseek(f, 0, SEEK_SET);
	s = (char*)malloc((size_t)n + 1);
	if (!s || fread(s, 1, (size_t)n, f) != (size_t)n) { fclose(f); free(s); return NULL; }
	fclose(f);
	for (i = j = 0; i < n; i++) if (s[i] != '\r') s[j++] = s[i];
	s[j] = 0;
	return s;
}

/* names [-check file] -once header... -- header...: the map on stdout, or
   exit 1 when the file differs from it (CR ignored) */
static int run_names(int argc, char** argv)
{
	names N;
	int i, status = 0;
	char const* check = NULL;
	char* text;
	memset(&N, 0, sizeof(N));
	for (i = 0; i < argc; i++)
	{
		toks v = { NULL, 0, 0 };
		if (!strcmp(argv[i], "-check") && i + 1 < argc) { check = argv[++i]; continue; }
		if (!strcmp(argv[i], "-once")) { N.collect_once = 1; continue; }
		if (!strcmp(argv[i], "--")) { N.collect_once = 0; continue; }
		if (!tokenize(argv[i], &v)) return 2;
		statements(&v, names_stmt, &N);
		free(v.t);
	}
	qsort(N.s, N.n, sizeof(char*), cmp);
	text = names_text(&N);
	if (!text) return 2;
	if (check)
	{
		char* old = read_text(check);
		status = !old || strcmp(old, text) != 0;
		if (status) printf("%s is not the map of the headers: regenerate it (qaws_api_tool names)\n", check);
		free(old);
	}
	else
		fputs(text, stdout);
	free(text);
	return status;
}

int main(int argc, char** argv)
{
	if (argc >= 2 && !strcmp(argv[1], "survey"))
		return run_survey(argc - 2, argv + 2);
	if (argc >= 2 && !strcmp(argv[1], "names"))
		return run_names(argc - 2, argv + 2);
	fprintf(stderr, "usage: qaws_api_tool survey [-v] [-order] [-allow file] [-types header... --] header...\n"
		"       qaws_api_tool names [-check file] -once header... -- header...\n");
	return 2;
}
