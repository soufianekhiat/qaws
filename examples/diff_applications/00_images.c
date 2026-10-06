/* ================================================================== */
/*  Images                                                            */
/* ================================================================== */

typedef struct image
{
	int w, h;
	float* rgb;   /* w * h * 3, 0..1 */
	float* lum;   /* w * h */
} image;

static int read_token(FILE* f, char* buf, int cap)
{
	int c, n = 0;
	do
	{
		c = fgetc(f);
		if (c == '#')
			while (c != '\n' && c != EOF)
				c = fgetc(f);
	} while (c == ' ' || c == '\n' || c == '\r' || c == '\t');
	while (c != EOF && c != ' ' && c != '\n' && c != '\r' && c != '\t' && n < cap - 1)
	{
		buf[n++] = (char)c;
		c = fgetc(f);
	}
	buf[n] = 0;
	return n;
}

static int image_load_ppm(char const* path, image* img)
{
	FILE* f = fopen(path, "rb");
	char tok[32];
	int maxv, i, n;
	unsigned char* raw;
	memset(img, 0, sizeof(*img));
	if (!f)
		return 0;
	read_token(f, tok, 32);
	if (strcmp(tok, "P6") != 0)
	{
		fclose(f);
		return 0;
	}
	read_token(f, tok, 32); img->w = atoi(tok);
	read_token(f, tok, 32); img->h = atoi(tok);
	read_token(f, tok, 32); maxv = atoi(tok);
	n = img->w * img->h;
	raw = (unsigned char*)malloc((size_t)n * 3);
	img->rgb = (float*)malloc(sizeof(float) * (size_t)n * 3);
	img->lum = (float*)malloc(sizeof(float) * (size_t)n);
	if (fread(raw, 1, (size_t)n * 3, f) != (size_t)n * 3)
	{
		fclose(f);
		free(raw);
		return 0;
	}
	fclose(f);
	for (i = 0; i < n * 3; i++)
		img->rgb[i] = raw[i] / (float)maxv;
	for (i = 0; i < n; i++)
		img->lum[i] = 0.299f * img->rgb[3 * i] + 0.587f * img->rgb[3 * i + 1] + 0.114f * img->rgb[3 * i + 2];
	free(raw);
	return 1;
}

static void image_free(image* img)
{
	free(img->rgb);
	free(img->lum);
}

/* Separable Gaussian blur of a w x h field. */
static void blur(float* f, int w, int h, double sigma)
{
	int r = (int)ceil(3 * sigma), x, y, k;
	float* tmp = (float*)malloc(sizeof(float) * (size_t)w * h);
	float* ker = (float*)malloc(sizeof(float) * (size_t)(2 * r + 1));
	double sum = 0;
	for (k = -r; k <= r; k++)
		sum += ker[k + r] = (float)exp(-k * k / (2 * sigma * sigma));
	for (k = 0; k <= 2 * r; k++)
		ker[k] = (float)(ker[k] / sum);
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		{
			double a = 0;
			for (k = -r; k <= r; k++)
			{
				int xx = x + k < 0 ? 0 : (x + k >= w ? w - 1 : x + k);
				a += ker[k + r] * f[y * w + xx];
			}
			tmp[y * w + x] = (float)a;
		}
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		{
			double a = 0;
			for (k = -r; k <= r; k++)
			{
				int yy = y + k < 0 ? 0 : (y + k >= h ? h - 1 : y + k);
				a += ker[k + r] * tmp[yy * w + x];
			}
			f[y * w + x] = (float)a;
		}
	free(tmp);
	free(ker);
}

/* Bilinear sample with its spatial gradient (pixel units, clamped). */
static double sample(float const* f, int w, int h, double x, double y, double* dx, double* dy)
{
	int x0, y0;
	double fx, fy, a, b, c, d;
	if (x < 0) x = 0;
	if (y < 0) y = 0;
	if (x > w - 1.001) x = w - 1.001;
	if (y > h - 1.001) y = h - 1.001;
	x0 = (int)x;
	y0 = (int)y;
	fx = x - x0;
	fy = y - y0;
	a = f[y0 * w + x0];
	b = f[y0 * w + x0 + 1];
	c = f[(y0 + 1) * w + x0];
	d = f[(y0 + 1) * w + x0 + 1];
	if (dx) *dx = (1 - fy) * (b - a) + fy * (d - c);
	if (dy) *dy = (1 - fx) * (c - a) + fx * (d - b);
	return (1 - fy) * ((1 - fx) * a + fx * b) + fy * ((1 - fx) * c + fx * d);
}

