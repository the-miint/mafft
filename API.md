# MAFFT C Library API

This document covers the C library interface for embedding MAFFT alignment
into other programs.  The library is built as `libmafft_parttree.a` and its
public header is `core/mafft_api.h`.

## Building

```bash
cd core
make libmafft_parttree.a
```

Link against it with `-lm -lpthread`:

```bash
gcc -o myprogram myprogram.c -Lcore -lmafft_parttree -lm -lpthread
```

Include the header:

```c
#include "mafft_api.h"
```

## Quick start

```c
#include <stdio.h>
#include "mafft_api.h"

int main(void)
{
    mafft_config_t cfg;
    mafft_config_init(&cfg);
    cfg.strategy = MAFFT_STRATEGY_PARTTREE;
    cfg.seqtype  = MAFFT_SEQ_DNA;

    mafft_ctx_t *ctx = mafft_create(&cfg);
    if (!ctx) { fprintf(stderr, "OOM\n"); return 1; }

    const char *names[] = { "s1", "s2", "s3" };
    const char *seqs[]  = {
        "ACGTACGTACGT",
        "ACGAACGTACGT",
        "ACGTACGAACGT"
    };

    mafft_output_t *out = NULL;
    mafft_stats_t stats;
    int rc = mafft_align(ctx, names, seqs, 3, &out, &stats);

    if (rc == MAFFT_OK)
    {
        for (int i = 0; i < out->n_seqs; i++)
            printf(">%s\n%s\n", out->names[i], out->seqs[i]);
        mafft_output_free(out);
    }
    else
    {
        fprintf(stderr, "%s: %s\n",
                mafft_strerror(rc), mafft_last_error(ctx));
    }

    mafft_destroy(ctx);
    return rc;
}
```

## Config

Initialize with `mafft_config_init()` which sets `struct_size` and all
defaults.  Modify fields as needed, then pass to `mafft_create()`.

```c
mafft_config_t cfg;
mafft_config_init(&cfg);
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `struct_size` | `size_t` | `sizeof(mafft_config_t)` | ABI version check. Do not set manually. |
| `strategy` | `int` | `MAFFT_STRATEGY_AUTO` | Alignment strategy (see below). |
| `seqtype` | `int` | `MAFFT_SEQ_AUTO` | `MAFFT_SEQ_DNA`, `_RNA`, `_PROTEIN`, or `_AUTO`. |
| `max_iterate` | `int` | 0 | Max refinement iterations. 0 = strategy default. |
| `retree` | `int` | 2 | Guide-tree rebuilding cycles. |
| `partsize` | `int` | 50 | PartTree partition bucket size. |
| `groupsize` | `int` | -1 | PartTree group size. -1 = njob+1 (single group). |
| `gap_open` | `double` | 0.0 | Gap opening penalty. 0.0 = strategy default. |
| `gap_extend` | `double` | 0.0 | Gap extension penalty. 0.0 = strategy default. |
| `offset` | `double` | 0.0 | Offset value. 0.0 = strategy default. |
| `n_threads` | `int` | 0 | Thread count. 0 = single-threaded. |
| `seed` | `int64_t` | 0 | Random seed for determinism. |
| `alloc_fn` | function pointer | `NULL` | Custom allocator for output. NULL = malloc. Must return 16-byte aligned. |
| `free_fn` | function pointer | `NULL` | Custom deallocator for output. NULL = free. |
| `alloc_ud` | `void *` | `NULL` | User data passed to alloc_fn/free_fn. |
| `progress_cb` | function pointer | `NULL` | Progress callback. Return non-zero to cancel. |
| `progress_ud` | `void *` | `NULL` | User data for progress callback. |
| `log_cb` | function pointer | `NULL` | Log message callback. |
| `log_ud` | `void *` | `NULL` | User data for log callback. |

## Strategies

| Constant | CLI equivalent | Description |
|----------|----------------|-------------|
| `MAFFT_STRATEGY_AUTO` | `mafft --auto` | Auto-select based on input size. |
| `MAFFT_STRATEGY_FFTNS2` | `mafft --retree 2` | Fast progressive (FFT-NS-2). |
| `MAFFT_STRATEGY_FFTNSI` | `mafft --retree 2 --maxiterate 2` | Progressive + refinement. |
| `MAFFT_STRATEGY_LINSI` | `mafft --localpair --maxiterate 1000` | Most accurate for small datasets. |
| `MAFFT_STRATEGY_GINSI` | `mafft --globalpair --maxiterate 1000` | Global pairwise + refinement. Requires LAST. |
| `MAFFT_STRATEGY_EINSI` | `mafft --genafpair --maxiterate 1000` | For sequences with large gaps. Requires LAST. |
| `MAFFT_STRATEGY_PARTTREE` | `mafft --parttree` | Fast for large datasets (10k+ seqs). |

**Implemented:** PARTTREE, FFTNS2, FFTNSI, AUTO, LINSI, GINSI, EINSI.

AUTO selects FFTNSI for small inputs (n<500, len<10k), FFTNS2 for medium
(n<200k), and PARTTREE for large datasets.

### External tool requirements

LINSI, GINSI, and EINSI require the [LAST aligner](https://gitlab.com/mcfrith/last)
(`lastdb` and `lastal` binaries) to be installed and available in `PATH`.
If these tools are not found, `mafft_align()` returns `MAFFT_ERR_INVALID_INPUT`
with an error message naming the missing tool and suggesting alternatives.

PARTTREE, FFTNS2, FFTNSI, and AUTO have no external dependencies.

## Context lifecycle

```c
mafft_ctx_t *ctx = mafft_create(&cfg);  /* NULL on OOM or bad config */
/* ... use ctx ... */
mafft_destroy(ctx);                      /* safe to pass NULL */
```

Contexts are lightweight.  Multiple contexts can coexist, but `mafft_align()`
is serialized by an internal mutex -- only one alignment runs at a time.
Callers can safely call `mafft_align()` from different threads without
external synchronization.

## Alignment

```c
int mafft_align(mafft_ctx_t *ctx,
                const char **names, const char **seqs, int n_seqs,
                mafft_output_t **out, mafft_stats_t *stats);
```

| Parameter | Description |
|-----------|-------------|
| `ctx` | Context from `mafft_create()`. |
| `names` | Array of `n_seqs` sequence name strings. |
| `seqs` | Array of `n_seqs` unaligned sequence strings. |
| `n_seqs` | Number of sequences. Must be >= 2. |
| `out` | On success, receives a pointer to the output struct. |
| `stats` | Optional (may be NULL). Receives timing and strategy info. |

Input strings are not modified.  The function copies them internally.

### Output struct

```c
typedef struct {
    int          n_seqs;       /* number of aligned sequences */
    int          aligned_len;  /* length of each aligned sequence (with gaps) */
    const char **names;        /* array of n_seqs name strings */
    const char **seqs;         /* array of n_seqs aligned sequence strings */
    void        *_base;        /* internal -- do not touch */
} mafft_output_t;
```

All data is packed into a single allocation.  Free with `mafft_output_free()`.

### Stats struct

```c
typedef struct {
    int     n_iterations;   /* refinement iterations performed */
    int     strategy_used;  /* actual strategy (for MAFFT_STRATEGY_AUTO) */
    double  elapsed_secs;   /* wall-clock time */
} mafft_stats_t;
```

## Error codes

| Code | Constant | Description |
|------|----------|-------------|
| 0 | `MAFFT_OK` | Success. |
| -1 | `MAFFT_ERR_NOMEM` | Out of memory. |
| -2 | `MAFFT_ERR_INVALID_INPUT` | Bad input (too few seqs, illegal chars, etc.). |
| -3 | `MAFFT_ERR_INTERNAL` | Internal error. |
| -4 | `MAFFT_ERR_CANCELLED` | Cancelled via progress callback. |
| -5 | `MAFFT_ERR_OUTPUT_TOO_LARGE` | Aligned length exceeded internal buffer. |

Use `mafft_strerror(code)` for a static category string, and
`mafft_last_error(ctx)` for a detailed message from the last failed call.

## Custom memory allocation

Supply `alloc_fn` and `free_fn` in the config to control output memory
allocation.  Internal working memory still uses `malloc`/`free`.

```c
cfg.alloc_fn = my_alloc;   /* must return 16-byte aligned */
cfg.free_fn  = my_free;
cfg.alloc_ud = my_context;
```

`mafft_output_free()` uses the allocator that was active when the output was
created.  No context reference is needed at free time.

## Thread safety

`mafft_align()` acquires an internal mutex.  Multiple threads can call it
concurrently -- calls are serialized automatically.  `mafft_create()`,
`mafft_destroy()`, `mafft_output_free()`, `mafft_strerror()`, and
`mafft_last_error()` are safe to call from any thread without locking.

## Legacy API

The original `splittbfast_library()` and `disttbfast()` functions in
`core/mafft.h` remain functional but are deprecated.  New code should use
`mafft_api.h`.
