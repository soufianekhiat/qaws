#include "qaws_internal_span.h"

unsigned int qaws_internal_find_span(
	qaws_scalar* out_local_t,
	qaws_curve const* curve,
	qaws_scalar parameter)
{
	unsigned int span_count = curve->span_count;
	qaws_scalar const* boundaries = curve->span_boundaries;
	unsigned int span;

	if (parameter <= boundaries[0]) {
		span = 0;
	} else if (parameter >= boundaries[span_count]) {
		span = span_count - 1;
	} else {
		unsigned int lo = 0;
		unsigned int hi = span_count - 1;

		while (lo < hi) {
			unsigned int mid = lo + (hi - lo) / 2;
			if (parameter < boundaries[mid + 1]) {
				hi = mid;
			} else {
				lo = mid + 1;
			}
		}
		span = lo;
	}

	if (out_local_t) {
		qaws_scalar span_start = boundaries[span];
		qaws_scalar span_end = boundaries[span + 1];
		qaws_scalar span_len = span_end - span_start;
		qaws_scalar local_t;

		if (span_len < (qaws_scalar)1e-10) {
			local_t = (qaws_scalar)0;
		} else {
			local_t = (parameter - span_start) / span_len;
		}

		if (local_t < (qaws_scalar)0) local_t = (qaws_scalar)0;
		else if (local_t > (qaws_scalar)1) local_t = (qaws_scalar)1;

		*out_local_t = local_t;
	}

	return span;
}
