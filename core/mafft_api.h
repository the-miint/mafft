#ifndef MAFFT_API_H
#define MAFFT_API_H

#include <stddef.h>
#include <stdint.h>

/* ---- Error codes (negative = error, 0 = success) ---- */
#define MAFFT_OK                    0
#define MAFFT_ERR_NOMEM            -1
#define MAFFT_ERR_INVALID_INPUT    -2
#define MAFFT_ERR_INTERNAL         -3
#define MAFFT_ERR_CANCELLED        -4
#define MAFFT_ERR_OUTPUT_TOO_LARGE -5

/* ---- Strategy constants ---- */
#define MAFFT_STRATEGY_AUTO      0
#define MAFFT_STRATEGY_FFTNS2    1
#define MAFFT_STRATEGY_FFTNSI    2
#define MAFFT_STRATEGY_LINSI     3
#define MAFFT_STRATEGY_GINSI     4
#define MAFFT_STRATEGY_EINSI     5
#define MAFFT_STRATEGY_PARTTREE  6

/* ---- Sequence type ---- */
#define MAFFT_SEQ_AUTO     0
#define MAFFT_SEQ_DNA      1
#define MAFFT_SEQ_RNA      2
#define MAFFT_SEQ_PROTEIN  3

/* ---- Config ---- */
typedef struct {
	size_t  struct_size;        /* set by mafft_config_init(). DO NOT set manually. */

	int     strategy;           /* MAFFT_STRATEGY_* constant */
	int     seqtype;            /* MAFFT_SEQ_* constant. AUTO scans input chars. */

	int     max_iterate;        /* reserved: iterative strategies (Phase 6) */
	int     retree;             /* reserved: FFT-NS tree rebuilding (Phase 5) */

	int     partsize;           /* PartTree bucket size (default 50) */
	int     groupsize;          /* PartTree group size. -1 = njob+1 (single group). */

	double  gap_open;           /* gap opening penalty. 0.0 = strategy default. */
	double  gap_extend;         /* reserved: progressive strategies (Phase 5) */
	double  offset;             /* offset value. 0.0 = strategy default. */

	int     n_threads;          /* 0 = single-threaded */
	int64_t seed;               /* reserved: deterministic iteration (Phase 6) */

	/* Custom allocator for output only.  Context and internal working
	 * memory always use system malloc/free.  Must set both or neither. */
	void *(*alloc_fn)(size_t size, void *alloc_ud);
	void  (*free_fn)(void *ptr, void *alloc_ud);
	void   *alloc_ud;

	int   (*progress_cb)(const char *stage, double frac, void *ud);
	void   *progress_ud;
	void  (*log_cb)(const char *msg, void *ud);
	void   *log_ud;
} mafft_config_t;

/* ---- Opaque context ---- */
typedef struct mafft_ctx mafft_ctx_t;

/* ---- Output (SOA, single backing allocation) ---- */
typedef struct {
	int          n_seqs;
	int          aligned_len;
	const char **names;
	const char **seqs;
	void        *_base;
} mafft_output_t;

/* ---- Stats ---- */
typedef struct {
	int     n_iterations;   /* reserved: always 0 until dvtditr exposes count */
	int     strategy_used;  /* actual strategy (for MAFFT_STRATEGY_AUTO) */
	double  elapsed_secs;
} mafft_stats_t;

/* ---- Lifecycle ---- */
void           mafft_config_init(mafft_config_t *cfg);
mafft_ctx_t   *mafft_create(const mafft_config_t *cfg);
void           mafft_destroy(mafft_ctx_t *ctx);

/* ---- Core computation ---- */
int  mafft_align(mafft_ctx_t *ctx,
                 const char **names, const char **seqs, int n_seqs,
                 mafft_output_t **out, mafft_stats_t *stats);

/* ---- Output cleanup ---- */
void mafft_output_free(mafft_output_t *out);

/* ---- Log access ---- */
/* Returns captured log text from the most recent mafft_align() call on ctx.
 * The pointer is valid until the next mafft_align() on the same context or
 * until mafft_destroy().  Returns "" if no log was captured.  Do not free. */
const char *mafft_ctx_log(const mafft_ctx_t *ctx);

/* ---- Error reporting ---- */
const char *mafft_strerror(int code);
const char *mafft_last_error(const mafft_ctx_t *ctx);

#endif /* MAFFT_API_H */
