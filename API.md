# MAFFT C Library API

Embed MAFFT multiple sequence alignment into C/C++ programs.  The library
is built as `libmafft_parttree.a` with public header `core/mafft_api.h`.

## Building

```bash
cd core
make libmafft_parttree.a
```

Link with `-lm -lpthread`:

```bash
gcc -o myprogram myprogram.c -Lcore -lmafft_parttree -lm -lpthread
```

## Quick start

```c
#include <stdio.h>
#include "mafft_api.h"

int main(void)
{
    /* Initialize config with defaults */
    mafft_config_t cfg;
    mafft_config_init(&cfg);

    /* Create context */
    mafft_ctx_t *ctx = mafft_create(&cfg);

    /* Input sequences */
    const char *names[] = { "Human", "Mouse", "Rat" };
    const char *seqs[]  = {
        "MKFLILLFNILCLFPVLAADNHGVS",
        "MKFLVLLFNILCLFPVLAADNHGVS",
        "MKFLILLFNILCLFPVLAADNHGVQ"
    };

    /* Align */
    mafft_output_t *out = NULL;
    mafft_stats_t stats;
    int rc = mafft_align(ctx, names, seqs, 3, &out, &stats);

    if (rc == MAFFT_OK)
    {
        printf("Strategy used: %d\n", stats.strategy_used);
        printf("Time: %.3f sec\n", stats.elapsed_secs);
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

## Lifecycle

```
mafft_config_init(&cfg)    Set defaults (must call before modifying cfg)
mafft_create(&cfg)         Create context (NULL on error)
mafft_align(ctx, ...)      Run alignment (thread-safe, serialized by mutex)
mafft_output_free(out)     Free output (safe to pass NULL)
mafft_ctx_log(ctx)         Get captured log text
mafft_destroy(ctx)         Destroy context (safe to pass NULL)
```

Contexts are lightweight.  Multiple contexts can coexist.  `mafft_align()`
acquires an internal mutex -- only one alignment runs at a time.  Callers
can safely call `mafft_align()` from different threads without external
synchronization; calls block until the mutex is available.

## Config

Initialize with `mafft_config_init()` which sets `struct_size` and all
defaults.  Modify fields as needed, then pass to `mafft_create()`.

```c
mafft_config_t cfg;
mafft_config_init(&cfg);
cfg.strategy = MAFFT_STRATEGY_FFTNS2;
cfg.seqtype  = MAFFT_SEQ_PROTEIN;
```

### Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `struct_size` | `size_t` | `sizeof(mafft_config_t)` | ABI version check. Do not set manually. |
| `strategy` | `int` | `MAFFT_STRATEGY_AUTO` | Alignment strategy (see [Strategies](#strategies)). |
| `seqtype` | `int` | `MAFFT_SEQ_AUTO` | Sequence type (see [Sequence type](#sequence-type)). |
| `max_iterate` | `int` | 0 | Max refinement iterations. 0 = use strategy default (2 for FFTNSI, 1000 for LINSI/GINSI/EINSI). |
| `retree` | `int` | 2 | Guide-tree rebuilding cycles (FFTNS2/FFTNSI). |
| `partsize` | `int` | 50 | PartTree partition bucket size. |
| `groupsize` | `int` | -1 | PartTree group size. -1 = njob+1 (single group). |
| `gap_open` | `double` | 0.0 | Gap opening penalty. 0.0 = strategy default (-1.53). |
| `gap_extend` | `double` | 0.0 | Gap extension penalty. 0.0 = strategy default. |
| `offset` | `double` | 0.0 | Offset value. 0.0 = strategy default. |
| `n_threads` | `int` | 0 | Thread count. 0 = single-threaded. |
| `seed` | `int64_t` | 0 | Random seed (reserved for future deterministic mode). |
| `alloc_fn` | `void *(*)(size_t, void *)` | `NULL` | Custom allocator for output. Must return 16-byte aligned. |
| `free_fn` | `void (*)(void *, void *)` | `NULL` | Custom deallocator for output. |
| `alloc_ud` | `void *` | `NULL` | User data passed to alloc_fn/free_fn. |
| `progress_cb` | `int (*)(const char *, double, void *)` | `NULL` | Progress callback. Return non-zero to cancel. |
| `progress_ud` | `void *` | `NULL` | User data for progress callback. |
| `log_cb` | `void (*)(const char *, void *)` | `NULL` | Log message callback (see [Log capture](#log-capture)). |
| `log_ud` | `void *` | `NULL` | User data for log callback. |

`alloc_fn` and `free_fn` must be both set or both NULL.  They control only
the output allocation (`mafft_output_t`).  Internal working memory and the
context itself always use system `malloc`/`free`.

## Strategies

| Constant | CLI equivalent | Speed | Accuracy | External deps |
|----------|----------------|-------|----------|---------------|
| `MAFFT_STRATEGY_AUTO` | `mafft --auto` | varies | varies | None |
| `MAFFT_STRATEGY_PARTTREE` | `mafft --parttree` | Very fast | Low | None |
| `MAFFT_STRATEGY_FFTNS2` | `mafft --retree 2` | Fast | Medium | None |
| `MAFFT_STRATEGY_FFTNSI` | `mafft --retree 2 --maxiterate 2` | Medium | Medium-High | None |
| `MAFFT_STRATEGY_LINSI` | `mafft --localpair --maxiterate 1000` | Slow | Highest | LAST |
| `MAFFT_STRATEGY_GINSI` | `mafft --globalpair --maxiterate 1000` | Slow | Very High | LAST |
| `MAFFT_STRATEGY_EINSI` | `mafft --genafpair --maxiterate 1000` | Slow | High (domains) | LAST |

### Choosing a strategy

- **< 500 sequences, < 10k residues each**: FFTNSI (good balance)
- **500 -- 200k sequences**: FFTNS2 (fast progressive)
- **> 200k sequences**: PARTTREE (scalable)
- **< 200 sequences, highest accuracy needed**: LINSI (requires LAST)
- **Sequences with large unalignable regions**: EINSI (requires LAST)
- **Don't know**: `MAFFT_STRATEGY_AUTO` (selects FFTNSI, FFTNS2, or PARTTREE)

### Auto strategy selection

AUTO selects based on input size:

| Condition | Strategy selected |
|-----------|-------------------|
| n_seqs < 500 AND max_len < 10,000 | FFTNSI |
| n_seqs < 200,000 | FFTNS2 |
| n_seqs >= 200,000 | PARTTREE |

AUTO never selects LINSI/GINSI/EINSI (external tool dependency).
`stats->strategy_used` reports the actual strategy chosen.

### External tool requirements

LINSI, GINSI, and EINSI require the
[LAST aligner](https://gitlab.com/mcfrith/last) (`lastdb` and `lastal`)
in `PATH`.  If not found, `mafft_align()` returns `MAFFT_ERR_INVALID_INPUT`
with a message naming the missing tool and suggesting FFTNSI/FFTNS2.

PARTTREE, FFTNS2, FFTNSI, and AUTO have no external dependencies.

## Sequence type

| Constant | Description |
|----------|-------------|
| `MAFFT_SEQ_AUTO` | Auto-detect by scanning input characters. DNA if >= 80% are A/C/G/T/U/N. |
| `MAFFT_SEQ_DNA` | DNA sequences. |
| `MAFFT_SEQ_RNA` | RNA sequences (treated as DNA internally). |
| `MAFFT_SEQ_PROTEIN` | Protein sequences. |

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
The context is reusable after both success and error.

### Output

```c
typedef struct {
    int          n_seqs;       /* number of aligned sequences */
    int          aligned_len;  /* length of each aligned sequence (with gaps) */
    const char **names;        /* array of n_seqs name strings */
    const char **seqs;         /* array of n_seqs aligned sequence strings */
    void        *_base;        /* internal -- do not touch */
} mafft_output_t;
```

All data is packed into a single allocation (SOA layout).  Free with
`mafft_output_free(out)`.  If a custom allocator was active when the
output was created, `mafft_output_free()` uses it automatically -- no
context reference needed at free time.

### Stats

```c
typedef struct {
    int     n_iterations;   /* reserved: always 0 (dvtditr does not expose count) */
    int     strategy_used;  /* actual strategy (useful when AUTO was requested) */
    double  elapsed_secs;   /* wall-clock computation time (excludes mutex wait) */
} mafft_stats_t;
```

## Error handling

### Error codes

| Code | Constant | Description |
|------|----------|-------------|
| 0 | `MAFFT_OK` | Success. |
| -1 | `MAFFT_ERR_NOMEM` | Out of memory. |
| -2 | `MAFFT_ERR_INVALID_INPUT` | Bad input, unsupported config, or missing external tool. |
| -3 | `MAFFT_ERR_INTERNAL` | Internal error (temp directory, file I/O). |
| -4 | `MAFFT_ERR_CANCELLED` | Cancelled via progress callback. |
| -5 | `MAFFT_ERR_OUTPUT_TOO_LARGE` | Aligned length exceeded internal buffer. |

### Error reporting

```c
const char *mafft_strerror(int code);          /* static category string */
const char *mafft_last_error(const mafft_ctx_t *ctx);  /* detailed message */
```

`mafft_last_error()` returns a detailed message from the most recent failed
call on `ctx`.  Returns `""` (never NULL) if no error occurred.  Both
functions are safe to call from any thread.

## Log capture

All diagnostic output (progress messages, version info, parameter summaries)
is captured internally.  Nothing is written to the caller's stderr or stdout.

### Retrieving logs

```c
const char *mafft_ctx_log(const mafft_ctx_t *ctx);
```

Returns the full captured log from the most recent `mafft_align()` call.
The pointer is valid until the next `mafft_align()` on the same context or
until `mafft_destroy()`.  Returns `""` if no log was captured.  Do not free.

### Log callback

Set `cfg.log_cb` to receive log messages line-by-line after alignment
completes:

```c
void my_log(const char *msg, void *ud)
{
    fprintf(stderr, "[mafft] %s\n", msg);
}

cfg.log_cb = my_log;
cfg.log_ud = NULL;  /* or your context pointer */
```

The callback receives each non-empty line of the captured log.  Messages
are temporary pointers valid only for the duration of the callback -- copy
if you need to retain them.  The callback is invoked under the internal
mutex, so it must not call `mafft_align()`.

When `log_cb` is NULL (default), logs are still captured and available via
`mafft_ctx_log()`.

## Thread safety

| Function | Thread-safe | Notes |
|----------|-------------|-------|
| `mafft_config_init()` | Yes | Pure function, no shared state. |
| `mafft_create()` | Yes | Allocates independent context. |
| `mafft_destroy()` | Yes | Only touches its own context. |
| `mafft_align()` | Yes | Serialized by internal mutex. Blocks if another call is running. |
| `mafft_output_free()` | Yes | Only touches its own allocation. |
| `mafft_ctx_log()` | Caller must ensure no concurrent `mafft_align()` on same ctx. |
| `mafft_strerror()` | Yes | Returns static string. |
| `mafft_last_error()` | Caller must ensure no concurrent `mafft_align()` on same ctx. |

MAFFT internally uses global mutable state.  True concurrent alignment
(multiple `mafft_align()` calls running simultaneously) is not supported.
The internal mutex serializes all calls.  The opaque context API is designed
to allow future migration to per-context state without breaking callers.

## Custom memory allocation

```c
cfg.alloc_fn = my_alloc;   /* must return 16-byte aligned pointer */
cfg.free_fn  = my_free;
cfg.alloc_ud = my_arena;
```

The custom allocator controls only the output allocation (`mafft_output_t`
and its packed data).  Internal working memory and the context itself use
system `malloc`/`free`.

`mafft_output_free()` uses the allocator that was active when the output
was created.  The allocator info is stored in a hidden header inside the
allocation, so no context reference is needed at free time.

Both `alloc_fn` and `free_fn` must be set, or both must be NULL.
`mafft_create()` returns NULL if only one is set.

## Iterative strategies and temp files

FFTNSI, LINSI, GINSI, and EINSI use a two-stage pipeline:

1. **Progressive alignment** (disttbfast or tbfast) -- produces initial
   alignment and writes intermediate files (`hat2`, `hat3`)
2. **Iterative refinement** (dvtditr) -- reads intermediate files and
   refines the alignment

The library manages this automatically using a private temp directory
created via `mkdtemp()` in `$TMPDIR` (falls back to `/tmp`).  All
intermediate files are cleaned up after the call returns, including on
error paths.

Note: `mafft_align()` temporarily changes the process working directory
(`chdir`) while holding the mutex.  Other threads performing relative-path
file I/O may be affected during this window.

## Legacy API

The original functions in `core/mafft.h` remain functional but are
deprecated:

```c
/* Deprecated -- use mafft_api.h instead */
int splittbfast_library(int ngui, int lgui, char **namegui, char **seqgui,
                        int argc, char **argv, int (*callback)(int, int, char*));
int disttbfast(...);
int tbfast_library(...);
int dvtditr_library(...);
const char *mafft_get_log(void);
void mafft_clear_log(void);
```

These provide direct access to individual alignment engines but lack
structured config, thread safety, log capture, and error reporting.
