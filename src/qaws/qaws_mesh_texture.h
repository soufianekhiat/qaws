#include "qaws_platform.h"   /* (QAWS_SCALAR_IS_FLOAT: the precision this pass declares) */
#if QAWS_SCALAR_IS_FLOAT
#ifndef QAWS_MESH_TEXTURE_H_F32
#define QAWS_MESH_TEXTURE_H_F32
#define QAWS_MESH_TEXTURE_H_BODY
#endif
#else
#ifndef QAWS_MESH_TEXTURE_H_F64
#define QAWS_MESH_TEXTURE_H_F64
#define QAWS_MESH_TEXTURE_H_BODY
#endif
#endif
#ifdef QAWS_MESH_TEXTURE_H_BODY
#undef QAWS_MESH_TEXTURE_H_BODY
#define QAWS_MESH_TEXTURE_H

#include "qaws_mesh.h"

/*
 * Texture attributes: an attribute held on a grid of texels over a mesh's
 * (u, v), at a resolution of its own rather than one value per vertex,
 * corner or face. The texels are the caller's array (no file is read or
 * written); a texture keeps the attribute's rule, so it is baked and read
 * by it (colours in their space, normals unit, quaternions on one side,
 * angles unwrapped, NEAREST labels never mixed).
 *
 * The (u, v) of the mesh: a corner (or vertex) attribute of 2 components,
 * "uv" unless named (qaws_mesh_3d_parameterize makes one). It is carried
 * through Booleans, remeshing and decimation as any corner attribute, so a
 * texture stays valid on the meshes made from it. With a face attribute of
 * texture ids (texture_id, one component), only the faces of one id are
 * baked or read: a Boolean of two textured meshes keeps each face's id.
 *
 * Texel (i, j) is at (j width + i) components in the array, its centre at
 * ((i + 0.5) / width, (j + 0.5) / height). Reading is bilinear between the
 * four nearest centres (NEAREST: the nearest texel), clamped to the edge
 * texels or repeated (repeat_u, repeat_v). Baking writes each texel whose
 * centre lies in a triangle the value there (the triangle's corner values
 * blended by the rule at the centre's barycentric weights; a triangle
 * across the edge of a repeated axis unwrapped to its first corner's
 * copy), then grows the charts by gutter rings of texels, each the mean of
 * its written neighbours (so bilinear reading and mipmaps do not pull in
 * texels of no chart); the others are left as they were.
 */

typedef struct qaws_texture
{
	qaws_scalar* texels;                /* width x height x components, row by row from v = 0 */
	unsigned int width;
	unsigned int height;
	unsigned int components;            /* 1 to 16 */
	qaws_attribute_rule rule;           /* how texels blend (baking sets it to the attribute's) */
	int repeat_u;                       /* nonzero: repeated along u (else clamped) */
	int repeat_v;
} qaws_texture;

typedef struct qaws_texture_desc
{
	char const* uv;                     /* the 2-component corner or vertex attribute of (u, v); NULL: "uv" */
	unsigned int gutter;                /* baking: rings of texels grown past the charts; 0: none */
	char const* texture_id;             /* NULL, or a 1-component face attribute: only the faces whose value is id */
	unsigned int id;
} qaws_texture_desc;

/* The vertex, corner or face attribute name baked into texture (its size and components given, its texels written, its
   rule set). out_coverage (may be NULL): width x height bytes, 1 a texel in a triangle, 2 in the gutter, 0 untouched.
   INVALID_ARGUMENT for no such attribute or uv, an edge attribute, components that differ, no texels or size; desc NULL:
   "uv", no gutter, every face. */
qaws_status qaws_mesh_3d_bake_texture(unsigned char* out_coverage, qaws_texture* texture, qaws_mesh_3d const* mesh, char const* name,
	qaws_texture_desc const* desc);

/* count points' values (components each) read from texture at uv (2 each) */
qaws_status qaws_texture_sample(qaws_scalar* out_values, qaws_texture const* texture, qaws_scalar const* uv, unsigned int count);

/* The corner attribute name of mesh read from texture at each corner's (u, v) (desc as for baking: with texture_id only
   those faces' corners, the others kept). A new attribute (the texture's components and rule, zeros elsewhere) unless
   one of that name is a corner attribute of as many components already (INVALID_ARGUMENT for another of that name). */
qaws_status qaws_mesh_3d_sample_texture(qaws_mesh_3d* mesh, qaws_texture const* texture, char const* name, qaws_texture_desc const* desc);

#endif /* QAWS_MESH_TEXTURE_H */
