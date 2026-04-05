#define _POSIX_C_SOURCE 200809L
#include "mltaln.h"
#include "mafft_api.h"
#include <string.h>
#include <time.h>
#include <pthread.h>

/* ------------------------------------------------------------------ */
/* Internal helpers                                                    */
/* ------------------------------------------------------------------ */

static pthread_mutex_t mafft_global_lock = PTHREAD_MUTEX_INITIALIZER;

struct mafft_ctx {
	mafft_config_t config;
	char           last_error[1024];
	/* Log buffer is allocated with malloc (not custom allocator).
	 * The custom allocator (alloc_fn/free_fn) is only used for
	 * mafft_output_t allocations returned to the caller. */
	char          *log_buf;
	size_t         log_len;
};

static void *mafft_alloc(mafft_ctx_t *ctx, size_t size)
{
	if( ctx->config.alloc_fn )
		return ctx->config.alloc_fn( size, ctx->config.alloc_ud );
	return malloc( size );
}

static void set_error(mafft_ctx_t *ctx, const char *fmt, ...)
{
	va_list args;
	va_start( args, fmt );
	vsnprintf( ctx->last_error, sizeof(ctx->last_error), fmt, args );
	va_end( args );
}

/* ---- fd-level log capture helpers ---- */

typedef struct {
	int saved_stderr;
	int saved_stdout;
	FILE *tmpfp;
} fd_capture_t;

static void fd_capture_start(fd_capture_t *cap)
{
	cap->saved_stderr = -1;
	cap->saved_stdout = -1;
	cap->tmpfp = tmpfile();
	if( cap->tmpfp )
	{
		int tmpfd = fileno( cap->tmpfp );
		int se = dup( STDERR_FILENO );
		int so = dup( STDOUT_FILENO );
		if( se < 0 || so < 0 )
		{
			/* fd table exhausted -- abort capture entirely */
			if( se >= 0 ) close( se );
			if( so >= 0 ) close( so );
			fclose( cap->tmpfp );
			cap->tmpfp = NULL;
			return;
		}
		cap->saved_stderr = se;
		cap->saved_stdout = so;
		dup2( tmpfd, STDERR_FILENO );
		dup2( tmpfd, STDOUT_FILENO );
	}
}

static char *fd_capture_end(fd_capture_t *cap)
{
	char *buf = NULL;

	if( cap->saved_stderr >= 0 )
	{
		long log_size;
		size_t nread;
		fflush( stderr );
		fflush( stdout );

		log_size = ftell( cap->tmpfp );
		if( log_size > 0 )
		{
			buf = (char *)malloc( log_size + 1 );
			if( buf )
			{
				rewind( cap->tmpfp );
				nread = fread( buf, 1, log_size, cap->tmpfp );
				buf[nread] = '\0';
			}
		}

		dup2( cap->saved_stderr, STDERR_FILENO );
		dup2( cap->saved_stdout, STDOUT_FILENO );
		close( cap->saved_stderr );
		close( cap->saved_stdout );
		fclose( cap->tmpfp );
	}
	return buf;
}

/* Store captured log in context and deliver to callback.
 * The callback receives temporary pointers into ctx->log_buf that are
 * valid only for the duration of the callback (in-place NUL splitting).
 * Must be called under mafft_global_lock. */
static void deliver_log(mafft_ctx_t *ctx, char *captured)
{
	free( ctx->log_buf );
	ctx->log_buf = captured;
	ctx->log_len = captured ? strlen( captured ) : 0;

	if( ctx->config.log_cb && ctx->log_buf )
	{
		char *p = ctx->log_buf;
		char *end;
		while( *p )
		{
			end = strchr( p, '\n' );
			if( end )
			{
				*end = '\0';
				if( *p )
					ctx->config.log_cb( p, ctx->config.log_ud );
				*end = '\n';
				p = end + 1;
			}
			else
			{
				if( *p )
					ctx->config.log_cb( p, ctx->config.log_ud );
				break;
			}
		}
	}
}

/* ---- Hidden header for self-describing output allocations ---- */

typedef struct {
	void (*free_fn)(void *ptr, void *ud);
	void  *alloc_ud;
} mafft_output_hidden_t;

/* Round up to guarantee mafft_output_t alignment after the hidden header. */
#define ALIGN_UP(x, a) ( ((x) + (a) - 1) & ~((a) - 1) )
#define OUTPUT_HIDDEN_SIZE \
	ALIGN_UP(sizeof(mafft_output_hidden_t), _Alignof(mafft_output_t))

/* Pack aligned sequences into a single-allocation SOA output struct.
 *
 * Layout inside _base:
 *   mafft_output_hidden_t  (free_fn + alloc_ud, padded for alignment)
 *   mafft_output_t         (returned pointer)
 *   const char *names[n]   (pointer array)
 *   const char *seqs[n]    (pointer array)
 *   name strings (null-terminated, packed)
 *   seq  strings (null-terminated, packed)
 */
static mafft_output_t *pack_output(mafft_ctx_t *ctx,
                                   char **names, char **seqs, int n)
{
	int i;
	size_t total;
	size_t name_bytes = 0;
	size_t seq_bytes = 0;
	int aligned_len;
	char *base;
	char *cursor;
	mafft_output_hidden_t *hidden;
	mafft_output_t *out;
	const char **name_ptrs;
	const char **seq_ptrs;

	if( n <= 0 ) return NULL;

	aligned_len = (int)strlen( seqs[0] );
	for( i = 0; i < n; i++ )
	{
		name_bytes += strlen( names[i] ) + 1;
		seq_bytes  += strlen( seqs[i] ) + 1;
	}

	total = OUTPUT_HIDDEN_SIZE
	      + sizeof(mafft_output_t)
	      + sizeof(const char *) * n
	      + sizeof(const char *) * n
	      + name_bytes
	      + seq_bytes;

	base = (char *)mafft_alloc( ctx, total );
	if( !base ) return NULL;

	hidden = (mafft_output_hidden_t *)base;
	hidden->free_fn  = ctx->config.free_fn;
	hidden->alloc_ud = ctx->config.alloc_ud;

	out = (mafft_output_t *)(base + OUTPUT_HIDDEN_SIZE);
	out->_base = base;
	out->n_seqs = n;
	out->aligned_len = aligned_len;

	name_ptrs = (const char **)((char *)out + sizeof(mafft_output_t));
	seq_ptrs  = name_ptrs + n;
	out->names = name_ptrs;
	out->seqs  = seq_ptrs;

	cursor = (char *)(seq_ptrs + n);

	for( i = 0; i < n; i++ )
	{
		size_t len = strlen( names[i] ) + 1;
		memcpy( cursor, names[i], len );
		name_ptrs[i] = cursor;
		cursor += len;
	}
	for( i = 0; i < n; i++ )
	{
		size_t len = strlen( seqs[i] ) + 1;
		memcpy( cursor, seqs[i], len );
		seq_ptrs[i] = cursor;
		cursor += len;
	}

	return out;
}

/* ---- Sequence type auto-detection ---- */

static int detect_seqtype(const char **seqs, int n)
{
	int i;
	int dna_chars = 0, total_chars = 0;

	for( i = 0; i < n; i++ )
	{
		const char *s = seqs[i];
		for( ; *s; s++ )
		{
			char c = *s;
			if( c >= 'a' && c <= 'z' ) c -= 32; /* uppercase */
			if( c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
			    c == 'U' || c == 'N' || c == '-' )
				dna_chars++;
			total_chars++;
		}
	}
	if( total_chars == 0 ) return MAFFT_SEQ_DNA;
	return ( dna_chars * 100 / total_chars >= 80 )
	       ? MAFFT_SEQ_DNA : MAFFT_SEQ_PROTEIN;
}

/* ---- Argv builder for splittbfast ---- */

static int build_parttree_argv(const mafft_config_t *cfg, int resolved_seqtype,
                               char ***argv_out, int *argc_out)
{
	char **av = NULL;
	int ac = 0;
	int cap = 20;
	int rc = MAFFT_OK;
	char tmp[64];

	av = (char **)calloc( cap, sizeof(char *) );
	if( !av ) return MAFFT_ERR_NOMEM;

#define PUSH_ARG(s) do { \
	if( ac >= cap ) { \
		char **newav; \
		cap *= 2; \
		newav = (char **)realloc( av, cap * sizeof(char *) ); \
		if( !newav ) { rc = MAFFT_ERR_NOMEM; goto argv_fail; } \
		av = newav; \
	} \
	av[ac] = strdup(s); \
	if( !av[ac] ) { rc = MAFFT_ERR_NOMEM; goto argv_fail; } \
	ac++; \
} while(0)

#define PUSH_ARG_FMT(fmt, val) do { \
	snprintf(tmp, sizeof(tmp), fmt, val); \
	PUSH_ARG(tmp); \
} while(0)

	PUSH_ARG( "splittbfast" );

	/* Sequence type -- always explicit, never NOTSPECIFIED */
	if( resolved_seqtype == MAFFT_SEQ_DNA || resolved_seqtype == MAFFT_SEQ_RNA )
		PUSH_ARG( "-D" );
	else
		PUSH_ARG( "-P" );

	/* Gap open penalty */
	PUSH_ARG( "-f" );
	if( cfg->gap_open != 0.0 )
		PUSH_ARG_FMT( "%.2f", cfg->gap_open );
	else
		PUSH_ARG( "-1.53" );

	/* SP score factor */
	PUSH_ARG( "-Q" );
	PUSH_ARG( "100" );

	/* Offset */
	PUSH_ARG( "-h" );
	if( cfg->offset != 0.0 )
		PUSH_ARG_FMT( "%.2f", cfg->offset );
	else
		PUSH_ARG( "0" );

	/* Partition size */
	PUSH_ARG( "-p" );
	if( cfg->partsize > 0 )
		PUSH_ARG_FMT( "%d", cfg->partsize );
	else
		PUSH_ARG( "50" );

	/* Group size: config default is -1 (njob+1).
	 * 0 is a valid user value (classsize=0), passed through. */
	PUSH_ARG( "-s" );
	PUSH_ARG_FMT( "%d", cfg->groupsize );

	/* Reorder output by similarity */
	PUSH_ARG( "-x" );

	/* Thread count */
	if( cfg->n_threads > 1 )
	{
		PUSH_ARG( "-C" );
		PUSH_ARG_FMT( "%d", cfg->n_threads );
	}

#undef PUSH_ARG
#undef PUSH_ARG_FMT

	*argv_out = av;
	*argc_out = ac;
	return MAFFT_OK;

argv_fail:
	{
		int j;
		for( j = 0; j < ac; j++ )
			free( av[j] );
		free( av );
	}
	*argv_out = NULL;
	*argc_out = 0;
	return rc;
}

/* Build internal argv for disttbfast (FFT-NS-2). */
static int build_fftns2_argv(const mafft_config_t *cfg, int resolved_seqtype,
                             char ***argv_out, int *argc_out)
{
	char **av = NULL;
	int ac = 0;
	int cap = 30;
	int rc = MAFFT_OK;
	char tmp[64];

	av = (char **)calloc( cap, sizeof(char *) );
	if( !av ) return MAFFT_ERR_NOMEM;

#define PUSH_ARG(s) do { \
	if( ac >= cap ) { \
		char **newav; \
		cap *= 2; \
		newav = (char **)realloc( av, cap * sizeof(char *) ); \
		if( !newav ) { rc = MAFFT_ERR_NOMEM; goto argv_fail; } \
		av = newav; \
	} \
	av[ac] = strdup(s); \
	if( !av[ac] ) { rc = MAFFT_ERR_NOMEM; goto argv_fail; } \
	ac++; \
} while(0)

#define PUSH_ARG_FMT(fmt, val) do { \
	snprintf(tmp, sizeof(tmp), fmt, val); \
	PUSH_ARG(tmp); \
} while(0)

	PUSH_ARG( "disttbfast" );

	/* Sequence type */
	if( resolved_seqtype == MAFFT_SEQ_DNA || resolved_seqtype == MAFFT_SEQ_RNA )
		PUSH_ARG( "-D" );
	else
		PUSH_ARG( "-P" );

	/* Retree cycles */
	PUSH_ARG( "-E" );
	PUSH_ARG_FMT( "%d", cfg->retree > 0 ? cfg->retree : 2 );

	/* Gap open penalty */
	PUSH_ARG( "-f" );
	if( cfg->gap_open != 0.0 )
		PUSH_ARG_FMT( "%.2f", cfg->gap_open );
	else
		PUSH_ARG( "-1.53" );

	/* Offset */
	PUSH_ARG( "-h" );
	if( cfg->offset != 0.0 )
		PUSH_ARG_FMT( "%.2f", cfg->offset );
	else
		PUSH_ARG( "0" );

	/* K-tuple size (6 for both DNA and protein in MAFFT's default) */
	PUSH_ARG( "-W" );
	PUSH_ARG( "6" );

	/* Thread count */
	if( cfg->n_threads > 1 )
	{
		PUSH_ARG( "-C" );
		PUSH_ARG_FMT( "%d", cfg->n_threads );
	}

#undef PUSH_ARG
#undef PUSH_ARG_FMT

	*argv_out = av;
	*argc_out = ac;
	return MAFFT_OK;

argv_fail:
	{
		int j;
		for( j = 0; j < ac; j++ )
			free( av[j] );
		free( av );
	}
	*argv_out = NULL;
	*argc_out = 0;
	return rc;
}

static void free_argv(char **av, int ac)
{
	int i;
	if( !av ) return;
	for( i = 0; i < ac; i++ )
		free( av[i] );
	free( av );
}

/* ------------------------------------------------------------------ */
/* Public API                                                          */
/* ------------------------------------------------------------------ */

void mafft_config_init(mafft_config_t *cfg)
{
	memset( cfg, 0, sizeof(*cfg) );
	cfg->struct_size  = sizeof(mafft_config_t);
	cfg->strategy     = MAFFT_STRATEGY_AUTO;
	cfg->seqtype      = MAFFT_SEQ_AUTO;
	cfg->max_iterate  = 0;   /* Phase 6: used by iterative strategies */
	cfg->retree       = 2;   /* Phase 5: used by FFT-NS strategies */
	cfg->partsize     = 50;
	cfg->groupsize    = -1;
	cfg->gap_open     = 0.0; /* 0.0 = use strategy-specific default */
	cfg->gap_extend   = 0.0; /* Phase 5: used by progressive strategies */
	cfg->offset       = 0.0; /* 0.0 = use strategy-specific default */
	cfg->n_threads    = 0;
	cfg->seed         = 0;   /* Phase 6: used for deterministic iteration */
	cfg->alloc_fn     = NULL;
	cfg->free_fn      = NULL;
	cfg->alloc_ud     = NULL;
	cfg->progress_cb  = NULL;
	cfg->progress_ud  = NULL;
	cfg->log_cb       = NULL;
	cfg->log_ud       = NULL;
}

mafft_ctx_t *mafft_create(const mafft_config_t *cfg)
{
	mafft_ctx_t *ctx;

	if( !cfg ) return NULL;
	if( cfg->struct_size != sizeof(mafft_config_t) ) return NULL;

	/* alloc_fn and free_fn must be both set or both NULL */
	if( (cfg->alloc_fn == NULL) != (cfg->free_fn == NULL) ) return NULL;

	ctx = (mafft_ctx_t *)calloc( 1, sizeof(mafft_ctx_t) );
	if( !ctx ) return NULL;

	ctx->config = *cfg;
	ctx->last_error[0] = '\0';
	ctx->log_buf = NULL;
	ctx->log_len = 0;
	mafft_library_mode = 1;

	return ctx;
}

void mafft_destroy(mafft_ctx_t *ctx)
{
	if( !ctx ) return;
	/* ctx and log_buf are always allocated with system malloc/calloc,
	 * not the custom allocator.  The custom allocator is only used
	 * for mafft_output_t allocations. */
	free( ctx->log_buf );
	free( ctx );
	/* mafft_library_mode is process-global.  Clear it so CLI code
	 * paths in the same process resume normal exit() behavior. */
	mafft_library_mode = 0;
}

int mafft_align(mafft_ctx_t *ctx,
                const char **names, const char **seqs, int n_seqs,
                mafft_output_t **out, mafft_stats_t *stats)
{
	int i, rc;
	int internal_argc = 0;
	char **internal_argv = NULL;
	char **work_names = NULL;
	char **work_seqs = NULL;
	size_t lgui;
	size_t max_input_len = 0;
	mafft_output_t *result = NULL;
	struct timespec t_start, t_end;
	int strategy;
	int resolved_seqtype;

	if( !ctx || !names || !seqs || !out || n_seqs < 2 )
	{
		if( ctx ) set_error( ctx, "Invalid arguments: need ctx, names, seqs, out, and n_seqs >= 2" );
		return MAFFT_ERR_INVALID_INPUT;
	}

	*out = NULL;
	if( stats ) memset( stats, 0, sizeof(*stats) );

	/* Resolve strategy */
	strategy = ctx->config.strategy;
	if( strategy == MAFFT_STRATEGY_AUTO )
		strategy = MAFFT_STRATEGY_PARTTREE; /* Phase 7 will add full auto logic */

	if( strategy != MAFFT_STRATEGY_PARTTREE && strategy != MAFFT_STRATEGY_FFTNS2 )
	{
		set_error( ctx, "Strategy %d not yet implemented", strategy );
		return MAFFT_ERR_INVALID_INPUT;
	}

	/* Resolve sequence type -- must be explicit before calling engine */
	resolved_seqtype = ctx->config.seqtype;
	if( resolved_seqtype == MAFFT_SEQ_AUTO )
		resolved_seqtype = detect_seqtype( seqs, n_seqs );

	pthread_mutex_lock( &mafft_global_lock );

	/* Timer starts after lock acquired -- measures computation, not wait */
	clock_gettime( CLOCK_MONOTONIC, &t_start );

	/* Compute max input length for buffer sizing (overflow-safe) */
	for( i = 0; i < n_seqs; i++ )
	{
		size_t len = strlen( seqs[i] );
		if( len > max_input_len ) max_input_len = len;
	}
	lgui = max_input_len * 3;
	if( lgui < 10000 ) lgui = 10000;
	if( lgui > (size_t)INT_MAX )
	{
		set_error( ctx, "Sequences too long for internal buffer (max_len=%zu)", max_input_len );
		rc = MAFFT_ERR_INVALID_INPUT;
		goto cleanup;
	}

	/* Build working copies (splittbfast modifies in-place) */
	work_names = (char **)calloc( n_seqs, sizeof(char *) );
	work_seqs  = (char **)calloc( n_seqs, sizeof(char *) );
	if( !work_names || !work_seqs )
	{
		set_error( ctx, "Out of memory allocating work arrays" );
		rc = MAFFT_ERR_NOMEM;
		goto cleanup;
	}
	for( i = 0; i < n_seqs; i++ )
	{
		work_names[i] = (char *)calloc( strlen(names[i]) + 1, sizeof(char) );
		work_seqs[i]  = (char *)calloc( lgui + 1, sizeof(char) );
		if( !work_names[i] || !work_seqs[i] )
		{
			set_error( ctx, "Out of memory copying sequence %d", i );
			rc = MAFFT_ERR_NOMEM;
			goto cleanup;
		}
		strcpy( work_names[i], names[i] );
		strcpy( work_seqs[i], seqs[i] );
	}

	/* Build argv and call engine */
	if( strategy == MAFFT_STRATEGY_PARTTREE )
	{
		rc = build_parttree_argv( &ctx->config, resolved_seqtype,
		                          &internal_argv, &internal_argc );
		if( rc != MAFFT_OK )
		{
			set_error( ctx, "Failed to build internal argv" );
			goto cleanup;
		}

		/* splittbfast_library() handles its own fd-level capture */
		rc = splittbfast_library( n_seqs, (int)lgui, work_names, work_seqs,
		                          internal_argc, internal_argv, NULL );

		/* Harvest log from splittbfast's internal buffer */
		{
			const char *log = mafft_get_log();
			char *captured = NULL;
			if( log && log[0] != '\0' )
			{
				captured = (char *)malloc( strlen(log) + 1 );
				if( captured ) strcpy( captured, log );
			}
			mafft_clear_log();
			deliver_log( ctx, captured );
		}
	}
	else if( strategy == MAFFT_STRATEGY_FFTNS2 )
	{
		fd_capture_t cap;
		char *captured;

		rc = build_fftns2_argv( &ctx->config, resolved_seqtype,
		                        &internal_argv, &internal_argc );
		if( rc != MAFFT_OK )
		{
			set_error( ctx, "Failed to build internal argv" );
			goto cleanup;
		}

		/* disttbfast has no built-in capture, so we wrap it */
		fd_capture_start( &cap );

		rc = disttbfast( n_seqs, (int)lgui, work_names, work_seqs,
		                 internal_argc, internal_argv, NULL );

		captured = fd_capture_end( &cap );
		deliver_log( ctx, captured );
	}

	/* Translate internal return code */
	if( rc == GUI_LENGTHOVER )
	{
		set_error( ctx, "Aligned length exceeds buffer (lgui=%d)", (int)lgui );
		rc = MAFFT_ERR_OUTPUT_TOO_LARGE;
		goto cleanup;
	}
	else if( rc == GUI_CANCEL )
	{
		set_error( ctx, "Alignment cancelled via callback" );
		rc = MAFFT_ERR_CANCELLED;
		goto cleanup;
	}
	else if( rc != 0 )
	{
		set_error( ctx, "Internal alignment error (code %d)", rc );
		rc = MAFFT_ERR_INTERNAL;
		goto cleanup;
	}

	/* Pack output */
	result = pack_output( ctx, work_names, work_seqs, n_seqs );
	if( !result )
	{
		set_error( ctx, "Out of memory packing output" );
		rc = MAFFT_ERR_NOMEM;
		goto cleanup;
	}

	*out = result;
	rc = MAFFT_OK;
	ctx->last_error[0] = '\0';

cleanup:
	free_argv( internal_argv, internal_argc );
	if( work_names )
	{
		for( i = 0; i < n_seqs; i++ )
			free( work_names[i] );
		free( work_names );
	}
	if( work_seqs )
	{
		for( i = 0; i < n_seqs; i++ )
			free( work_seqs[i] );
		free( work_seqs );
	}

	clock_gettime( CLOCK_MONOTONIC, &t_end );

	pthread_mutex_unlock( &mafft_global_lock );

	if( stats )
	{
		stats->strategy_used = strategy;
		stats->n_iterations  = 0; /* PartTree does not iterate */
		stats->elapsed_secs  = (t_end.tv_sec - t_start.tv_sec)
		                     + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;
	}

	return rc;
}

void mafft_output_free(mafft_output_t *out)
{
	mafft_output_hidden_t *hidden;
	if( !out ) return;
	hidden = (mafft_output_hidden_t *)out->_base;
	if( hidden->free_fn )
		hidden->free_fn( out->_base, hidden->alloc_ud );
	else
		free( out->_base );
}

const char *mafft_ctx_log(const mafft_ctx_t *ctx)
{
	if( !ctx || !ctx->log_buf ) return "";
	return ctx->log_buf;
}

const char *mafft_strerror(int code)
{
	switch( code )
	{
		case MAFFT_OK:                   return "Success";
		case MAFFT_ERR_NOMEM:            return "Out of memory";
		case MAFFT_ERR_INVALID_INPUT:    return "Invalid input";
		case MAFFT_ERR_INTERNAL:         return "Internal error";
		case MAFFT_ERR_CANCELLED:        return "Cancelled";
		case MAFFT_ERR_OUTPUT_TOO_LARGE: return "Output too large";
		default:                          return "Unknown error";
	}
}

const char *mafft_last_error(const mafft_ctx_t *ctx)
{
	if( !ctx ) return "";
	return ctx->last_error[0] ? ctx->last_error : "";
}
