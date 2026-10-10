#ifndef QAWS_INTERNAL_CLIP_H
#define QAWS_INTERNAL_CLIP_H

#include "../qaws_clip.h"
#include "../qaws_clip64.h"

/* Moves everything of src (paths, open paths, curves, vertices) to the end
   of dst and destroys src; parents keep pointing at their own paths. */
qaws_status qaws_internal_clip_result_append(qaws_clip_result* dst, qaws_clip_result* src);
qaws_status qaws_internal_clip64_result_append(qaws_clip64_result* dst, qaws_clip64_result* src);

/* An empty result. */
qaws_status qaws_internal_clip_result_empty(qaws_clip_result** out);
qaws_status qaws_internal_clip64_result_empty(qaws_clip64_result** out);

#endif /* QAWS_INTERNAL_CLIP_H */
