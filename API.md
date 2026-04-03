# MAFFT C Library API

This document covers the C library interface for embedding MAFFT's PartTree
alignment algorithm into other programs.  The library is built as
`libmafft_parttree.a` and its public header is `core/mafft.h`.

## Building

```bash
cd core
make libmafft_parttree.a
```

This produces a static library.  Link against it with `-lm -lpthread`:

```bash
gcc -o myprogram myprogram.c -Lcore -lmafft_parttree -lm -lpthread
```

Include the header:

```c
#include "mafft.h"
```

## Functions

### splittbfast_library

```c
int splittbfast_library(
    int    ngui,
    int    lgui,
    char **namegui,
    char **seqgui,
    int    argc,
    char **argv,
    int  (*callback)(int, int, char *)
);
```

Run a PartTree multiple sequence alignment.  Sequences are aligned in-place:
on success each `seqgui[i]` is replaced with its aligned (gap-inserted)
version.

**Parameters**

| Parameter  | Description |
|------------|-------------|
| `ngui`     | Number of input sequences.  Must be >= 2. |
| `lgui`     | Maximum buffer length of each `seqgui[i]` entry.  The aligned sequences may be longer than the originals due to gap insertion; if the aligned length exceeds `lgui` the function returns `GUI_LENGTHOVER` instead of overflowing the buffer. |
| `namegui`  | Array of `ngui` sequence name strings. |
| `seqgui`   | Array of `ngui` sequence strings (unaligned on input, aligned on output).  Each buffer must be at least `lgui + 1` bytes. |
| `argc`     | Length of `argv`. |
| `argv`     | Algorithm arguments passed as a string array, the same flags that the internal `splittbfast` binary accepts (see [Arguments](#arguments) below).  `argv[0]` should be `"splittbfast"`. |
| `callback` | Optional progress callback, or `NULL`.  See [Callback](#callback). |

**Return value**

| Value             | Meaning |
|-------------------|---------|
| `0`               | Success.  `seqgui` contains the aligned sequences. |
| `GUI_ERROR` (1)   | Fatal error (bad input, sequence too short, illegal characters, etc.). |
| `GUI_LENGTHOVER` (2) | Aligned length exceeds `lgui`.  Retry with a larger buffer. |
| `GUI_CANCEL` (3)  | Cancelled via callback. |

**Log capture**

When `ngui > 0` (library mode), all diagnostic output that MAFFT would
normally print to stderr and stdout is captured internally.  Nothing is
written to the caller's stderr or stdout.  Use `mafft_get_log()` to
retrieve the captured messages after the call returns.

---

### disttbfast

```c
int disttbfast(
    int    ngui,
    int    lgui,
    char **namegui,
    char **seqgui,
    int    argc,
    char **argv,
    int  (*callback)(int, int, char *)
);
```

Run a progressive (FFT-NS / distance-based) multiple sequence alignment.
Same calling convention as `splittbfast_library`.

> **Note:** `disttbfast` does not yet capture log output.  Diagnostic
> messages may still appear on stderr.

---

### mafft_get_log

```c
const char *mafft_get_log(void);
```

Return a pointer to the captured log output from the most recent
`splittbfast_library()` call.  The string is null-terminated.  Returns
`""` (empty string, never `NULL`) if no log has been captured.

The returned pointer is valid until the next call to `splittbfast_library()`
or `mafft_clear_log()`.  Do not `free()` it.

---

### mafft_clear_log

```c
void mafft_clear_log(void);
```

Free the internal log buffer.  After this call, `mafft_get_log()` returns
`""`.  This is optional -- `splittbfast_library()` clears the previous log
automatically before each run.

---

## Arguments

The `argv` array mirrors the flags accepted by the internal `splittbfast`
binary.  Common combinations:

| argv entry   | Meaning |
|-------------|---------|
| `"-D"`      | DNA mode. |
| `"-f"`, `"-1.53"` | Gap open penalty. |
| `"-Q"`, `"100"`   | SP score factor. |
| `"-h"`, `"0"`     | Offset (gap extension modifier). |
| `"-p"`, `"50"`    | Partition size (PartTree bucket size). |
| `"-s"`, `"-1"`    | Group size.  `-1` means `njob + 1` (single group, full alignment). |

A minimal DNA PartTree invocation:

```c
char *argv[] = {
    "splittbfast",
    "-D",           /* DNA */
    "-f", "-1.53",  /* gap open */
    "-Q", "100",    /* spfactor */
    "-h", "0",      /* offset */
    "-p", "50",     /* partsize */
    "-s", "-1"      /* groupsize = njob+1 */
};
int argc = 12;
```

---

## Callback

The optional callback allows the caller to monitor progress and cancel a
running alignment.

```c
int my_callback(int reserved, int percent, char *stage);
```

| Parameter   | Description |
|-------------|-------------|
| `reserved`  | Currently always 0. |
| `percent`   | Progress estimate, 0 -- 100. |
| `stage`     | Human-readable phase name (e.g. `"Distance matrix"`, `"Guide tree"`, `"Progressive alignment"`). |

**Return value:** return non-zero to cancel the alignment (the library
function will return `GUI_CANCEL`).  Return 0 to continue.

The callback is currently invoked by `disttbfast` at major phase
transitions.  `splittbfast_library` accepts the parameter for API
compatibility but does not call it.

---

## Error codes

Defined in `mafft.h`:

```c
#define GUI_ERROR      1   /* fatal error         */
#define GUI_LENGTHOVER 2   /* output buffer too small */
#define GUI_CANCEL     3   /* cancelled via callback  */
```

---

## Thread safety

The library uses global mutable state internally.  It is **not**
thread-safe.  Do not call any library function from multiple threads
concurrently.  Sequential calls from the same thread are safe -- all
internal state is reset between invocations.

---

## Complete example

```c
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mafft.h"

int main(void)
{
    int i;
    int n = 3;
    int l = 10000;  /* max aligned length */

    /* Allocate name and sequence arrays */
    char **seq  = (char **)calloc(n, sizeof(char *));
    char **name = (char **)calloc(n, sizeof(char *));
    for (i = 0; i < n; i++) seq[i]  = calloc(l + 1, sizeof(char));
    for (i = 0; i < n; i++) name[i] = calloc(100, sizeof(char));

    /* Fill in sequences */
    strcpy(name[0], "s1");  strcpy(seq[0], "ACGTACGTACGT");
    strcpy(name[1], "s2");  strcpy(seq[1], "ACGAACGTACGT");
    strcpy(name[2], "s3");  strcpy(seq[2], "ACGTACGAACGT");

    /* Build argv for DNA PartTree alignment */
    int argc = 12;
    char *argv[] = {
        "splittbfast",
        "-D", "-f", "-1.53", "-Q", "100",
        "-h", "0", "-p", "50", "-s", "-1"
    };

    /* Run alignment -- nothing is printed to stderr/stdout */
    int res = splittbfast_library(n, l, name, seq, argc, argv, NULL);

    if (res == 0)
    {
        /* seq[0..n-1] now contain aligned sequences */
        for (i = 0; i < n; i++)
            printf(">%s\n%s\n", name[i], seq[i]);
    }
    else
    {
        fprintf(stderr, "Alignment failed (code %d)\n", res);
        /* Check log for diagnostic detail */
        fprintf(stderr, "%s\n", mafft_get_log());
    }

    /* Optional: inspect diagnostic log */
    const char *log = mafft_get_log();
    if (log[0] != '\0')
    {
        /* log contains progress messages, version info, etc. */
    }

    mafft_clear_log();

    /* Cleanup */
    for (i = 0; i < n; i++) free(seq[i]);
    free(seq);
    for (i = 0; i < n; i++) free(name[i]);
    free(name);

    return res;
}
```

Compile and run:

```bash
gcc -o example example.c -Lcore -lmafft_parttree -lm -lpthread
./example
```
