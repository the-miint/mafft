#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include "mafft_api.h"

static int n_pass = 0;
static int n_fail = 0;

#define CHECK(cond, msg) do { \
	if (!(cond)) { \
		fprintf(stderr, "FAIL: %s (line %d)\n", msg, __LINE__); \
		n_fail++; \
	} else { \
		n_pass++; \
	} \
} while(0)

/* ---- Test data ---- */

static const char *dna_names[] = { "s1", "s2", "s3" };
static const char *dna_seqs[]  = {
	"ACGTACGTACGT",
	"ACGAACGTACGT",
	"ACGTACGAACGT"
};

static const char *protein_names[] = { "p1", "p2", "p3" };
static const char *protein_seqs[]  = {
	"MKFLILLFNILCLFPVLAADNHGVS",
	"MKFLVLLFNILCLFPVLAADNHGVS",
	"MKFLILLFNILCLFPVLAADNHGVQ"
};

/* ---- Native comparison helper ---- */

static int run_native_parttree(const char **names, const char **seqs, int n,
                               const char *seqtype_flag,
                               char ***out_seqs, int *out_n)
{
	char tmpfile_in[256], tmpfile_out[256];
	FILE *fp;
	int i;
	char line[65536];
	int count = 0;
	int cap = 0;
	char **result = NULL;
	char *cur_seq = NULL;

	snprintf(tmpfile_in, sizeof(tmpfile_in), "/tmp/mafft_test_in_%d.fa", (int)getpid());
	snprintf(tmpfile_out, sizeof(tmpfile_out), "/tmp/mafft_test_out_%d.fa", (int)getpid());

	fp = fopen(tmpfile_in, "w");
	if (!fp) return -1;
	for (i = 0; i < n; i++)
		fprintf(fp, ">%s\n%s\n", names[i], seqs[i]);
	fclose(fp);

	{
		char cmd[1024];
		snprintf(cmd, sizeof(cmd),
			"./splittbfast %s -f -1.53 -Q 100 -h 0 -p 50 -s -1 -x"
			" < %s > %s 2>/dev/null",
			seqtype_flag, tmpfile_in, tmpfile_out);
		if (system(cmd) != 0)
		{
			remove(tmpfile_in);
			remove(tmpfile_out);
			return -1;
		}
	}

	fp = fopen(tmpfile_out, "r");
	if (!fp) { remove(tmpfile_in); return -1; }

	while (fgets(line, sizeof(line), fp))
	{
		int len = strlen(line);
		while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
			line[--len] = '\0';

		if (line[0] == '>')
		{
			if (cur_seq)
			{
				if (count >= cap)
				{
					char **tmp;
					cap = cap ? cap * 2 : 16;
					tmp = (char **)realloc(result, cap * sizeof(char *));
					if (!tmp) { free(cur_seq); break; }
					result = tmp;
				}
				result[count++] = cur_seq;
				cur_seq = NULL;
			}
			cur_seq = calloc(1, 1);
		}
		else if (cur_seq && len > 0)
		{
			size_t old_len = strlen(cur_seq);
			char *tmp = (char *)realloc(cur_seq, old_len + len + 1);
			if (!tmp) break;
			cur_seq = tmp;
			memcpy(cur_seq + old_len, line, len + 1);
		}
	}
	if (cur_seq)
	{
		if (count >= cap)
		{
			char **tmp;
			cap = cap ? cap * 2 : 16;
			tmp = (char **)realloc(result, cap * sizeof(char *));
			if (tmp) result = tmp;
		}
		if (result) result[count++] = cur_seq;
		else free(cur_seq);
	}

	fclose(fp);
	remove(tmpfile_in);
	remove(tmpfile_out);

	*out_seqs = result;
	*out_n = count;
	return 0;
}

/* Compare library output vs native as sets (order may differ) */
static int compare_as_sets(const mafft_output_t *lib, char **native, int native_n)
{
	int i, j;
	if (lib->n_seqs != native_n) return 0;
	for (i = 0; i < native_n; i++)
	{
		int found = 0;
		for (j = 0; j < lib->n_seqs; j++)
		{
			if (strcmp(native[i], lib->seqs[j]) == 0) { found = 1; break; }
		}
		if (!found) return 0;
	}
	return 1;
}

/* ---- Tests ---- */

static void test_abi_size(void)
{
	mafft_config_t cfg;
	mafft_config_init(&cfg);
	CHECK(cfg.struct_size == sizeof(mafft_config_t),
	      "config.struct_size matches sizeof(mafft_config_t)");
}

static void test_config_defaults(void)
{
	mafft_config_t cfg;
	mafft_config_init(&cfg);

	CHECK(cfg.strategy == MAFFT_STRATEGY_AUTO, "default strategy is AUTO");
	CHECK(cfg.seqtype == MAFFT_SEQ_AUTO, "default seqtype is AUTO");
	CHECK(cfg.retree == 2, "default retree is 2");
	CHECK(cfg.partsize == 50, "default partsize is 50");
	CHECK(cfg.groupsize == -1, "default groupsize is -1");
	CHECK(cfg.gap_open == 0.0, "default gap_open is 0.0");
	CHECK(cfg.n_threads == 0, "default n_threads is 0");
	CHECK(cfg.seed == 0, "default seed is 0");
	CHECK(cfg.alloc_fn == NULL, "default alloc_fn is NULL");
	CHECK(cfg.free_fn == NULL, "default free_fn is NULL");
	CHECK(cfg.log_cb == NULL, "default log_cb is NULL");
	CHECK(cfg.progress_cb == NULL, "default progress_cb is NULL");
}

static void test_create_destroy(void)
{
	mafft_config_t cfg;
	mafft_config_init(&cfg);

	mafft_ctx_t *ctx = mafft_create(&cfg);
	CHECK(ctx != NULL, "mafft_create returns non-NULL");
	mafft_destroy(ctx);

	ctx = mafft_create(NULL);
	CHECK(ctx == NULL, "mafft_create(NULL) returns NULL");

	cfg.struct_size = 1;
	ctx = mafft_create(&cfg);
	CHECK(ctx == NULL, "mafft_create rejects bad struct_size");
}

static void test_half_paired_allocator(void)
{
	mafft_config_t cfg;
	mafft_config_init(&cfg);

	/* alloc_fn set, free_fn NULL */
	cfg.alloc_fn = (void *(*)(size_t, void *))malloc;
	cfg.free_fn = NULL;
	mafft_ctx_t *ctx = mafft_create(&cfg);
	CHECK(ctx == NULL, "mafft_create rejects alloc_fn without free_fn");

	/* free_fn set, alloc_fn NULL */
	mafft_config_init(&cfg);
	cfg.alloc_fn = NULL;
	cfg.free_fn = (void (*)(void *, void *))free;
	ctx = mafft_create(&cfg);
	CHECK(ctx == NULL, "mafft_create rejects free_fn without alloc_fn");
}

static void test_strerror(void)
{
	CHECK(strcmp(mafft_strerror(MAFFT_OK), "Success") == 0,
	      "strerror(OK)");
	CHECK(strcmp(mafft_strerror(MAFFT_ERR_NOMEM), "Out of memory") == 0,
	      "strerror(NOMEM)");
	CHECK(strcmp(mafft_strerror(MAFFT_ERR_INVALID_INPUT), "Invalid input") == 0,
	      "strerror(INVALID_INPUT)");
	CHECK(strcmp(mafft_strerror(MAFFT_ERR_INTERNAL), "Internal error") == 0,
	      "strerror(INTERNAL)");
	CHECK(strcmp(mafft_strerror(MAFFT_ERR_CANCELLED), "Cancelled") == 0,
	      "strerror(CANCELLED)");
	CHECK(strcmp(mafft_strerror(MAFFT_ERR_OUTPUT_TOO_LARGE), "Output too large") == 0,
	      "strerror(OUTPUT_TOO_LARGE)");
	CHECK(strcmp(mafft_strerror(999), "Unknown error") == 0,
	      "strerror(unknown)");
}

static void test_last_error_null(void)
{
	const char *err = mafft_last_error(NULL);
	CHECK(err != NULL && err[0] == '\0', "mafft_last_error(NULL) returns empty string");
}

static void test_bad_input(void)
{
	mafft_config_t cfg;
	mafft_config_init(&cfg);
	cfg.strategy = MAFFT_STRATEGY_PARTTREE;
	cfg.seqtype = MAFFT_SEQ_DNA;
	mafft_ctx_t *ctx = mafft_create(&cfg);
	mafft_output_t *out = NULL;
	int rc;

	rc = mafft_align(ctx, dna_names, dna_seqs, 1, &out, NULL);
	CHECK(rc == MAFFT_ERR_INVALID_INPUT, "n_seqs=1 returns INVALID_INPUT");
	CHECK(out == NULL, "output is NULL on error");
	CHECK(strlen(mafft_last_error(ctx)) > 0, "last_error is set on error");

	rc = mafft_align(ctx, dna_names, NULL, 3, &out, NULL);
	CHECK(rc == MAFFT_ERR_INVALID_INPUT, "NULL seqs returns INVALID_INPUT");

	mafft_destroy(ctx);
}

static void test_parttree_align(void)
{
	mafft_config_t cfg;
	mafft_config_init(&cfg);
	cfg.strategy = MAFFT_STRATEGY_PARTTREE;
	cfg.seqtype = MAFFT_SEQ_DNA;
	mafft_ctx_t *ctx = mafft_create(&cfg);
	mafft_output_t *out = NULL;
	mafft_stats_t stats;
	int rc, i;

	rc = mafft_align(ctx, dna_names, dna_seqs, 3, &out, &stats);
	CHECK(rc == MAFFT_OK, "parttree align returns OK");
	CHECK(out != NULL, "output is non-NULL");

	if (out)
	{
		CHECK(out->n_seqs == 3, "output has 3 sequences");
		CHECK(out->aligned_len > 0, "aligned_len > 0");
		CHECK(out->names != NULL, "names array is non-NULL");
		CHECK(out->seqs != NULL, "seqs array is non-NULL");
		CHECK(out->_base != NULL, "_base is non-NULL");

		for (i = 0; i < out->n_seqs; i++)
			CHECK((int)strlen(out->seqs[i]) == out->aligned_len,
			      "seq length matches aligned_len");

		CHECK((char *)out->names >= (char *)out->_base
		   && (char *)out->names <  (char *)out->_base + 1000000,
		      "names array inside _base allocation");

		mafft_output_free(out);
	}

	CHECK(stats.strategy_used == MAFFT_STRATEGY_PARTTREE,
	      "stats.strategy_used is PARTTREE");
	CHECK(stats.elapsed_secs >= 0.0, "elapsed_secs >= 0");
	CHECK(stats.n_iterations == 0, "no iterations for PARTTREE");

	mafft_destroy(ctx);
}

static void test_seq_auto_detect(void)
{
	mafft_config_t cfg;
	mafft_config_init(&cfg);
	cfg.strategy = MAFFT_STRATEGY_PARTTREE;
	cfg.seqtype = MAFFT_SEQ_AUTO;

	mafft_ctx_t *ctx = mafft_create(&cfg);
	mafft_output_t *out = NULL;
	int rc;

	/* DNA input with AUTO seqtype */
	rc = mafft_align(ctx, dna_names, dna_seqs, 3, &out, NULL);
	CHECK(rc == MAFFT_OK, "SEQ_AUTO with DNA input returns OK");
	if (out) mafft_output_free(out);
	out = NULL;

	/* Protein input with AUTO seqtype */
	rc = mafft_align(ctx, protein_names, protein_seqs, 3, &out, NULL);
	CHECK(rc == MAFFT_OK, "SEQ_AUTO with protein input returns OK");
	if (out)
	{
		CHECK(out->n_seqs == 3, "SEQ_AUTO protein: 3 output sequences");
		mafft_output_free(out);
	}

	mafft_destroy(ctx);
}

static void test_parttree_vs_native(void)
{
	mafft_config_t cfg;
	mafft_config_init(&cfg);
	cfg.strategy = MAFFT_STRATEGY_PARTTREE;
	cfg.seqtype = MAFFT_SEQ_DNA;
	mafft_ctx_t *ctx = mafft_create(&cfg);
	mafft_output_t *out = NULL;
	char **native_seqs = NULL;
	int native_n = 0;
	int rc, i;

	rc = mafft_align(ctx, dna_names, dna_seqs, 3, &out, NULL);
	if (rc != MAFFT_OK || !out)
	{
		CHECK(0, "library align failed, cannot compare to native");
		mafft_destroy(ctx);
		return;
	}

	rc = run_native_parttree(dna_names, dna_seqs, 3, "-D",
	                         &native_seqs, &native_n);
	if (rc != 0 || native_n == 0)
	{
		CHECK(0, "native splittbfast failed, cannot compare");
		mafft_output_free(out);
		mafft_destroy(ctx);
		return;
	}

	CHECK(native_n == out->n_seqs, "native and library produce same seq count");
	CHECK(compare_as_sets(out, native_seqs, native_n),
	      "DNA: library output matches native splittbfast");

	for (i = 0; i < native_n; i++) free(native_seqs[i]);
	free(native_seqs);
	mafft_output_free(out);
	mafft_destroy(ctx);
}

static void test_protein_vs_native(void)
{
	mafft_config_t cfg;
	mafft_config_init(&cfg);
	cfg.strategy = MAFFT_STRATEGY_PARTTREE;
	cfg.seqtype = MAFFT_SEQ_PROTEIN;
	mafft_ctx_t *ctx = mafft_create(&cfg);
	mafft_output_t *out = NULL;
	char **native_seqs = NULL;
	int native_n = 0;
	int rc, i;

	rc = mafft_align(ctx, protein_names, protein_seqs, 3, &out, NULL);
	if (rc != MAFFT_OK || !out)
	{
		CHECK(0, "protein library align failed");
		mafft_destroy(ctx);
		return;
	}

	rc = run_native_parttree(protein_names, protein_seqs, 3, "-P",
	                         &native_seqs, &native_n);
	if (rc != 0 || native_n == 0)
	{
		CHECK(0, "native protein splittbfast failed");
		mafft_output_free(out);
		mafft_destroy(ctx);
		return;
	}

	CHECK(native_n == out->n_seqs, "protein: native and library same seq count");
	CHECK(compare_as_sets(out, native_seqs, native_n),
	      "protein: library output matches native splittbfast");

	for (i = 0; i < native_n; i++) free(native_seqs[i]);
	free(native_seqs);
	mafft_output_free(out);
	mafft_destroy(ctx);
}

/* ---- Concurrency test ---- */

typedef struct {
	int passed;
	int n_seqs_out;
	int aligned_len_out;
} thread_arg_t;

static void *align_thread(void *arg)
{
	thread_arg_t *ta = (thread_arg_t *)arg;
	mafft_config_t cfg;
	mafft_config_init(&cfg);
	cfg.strategy = MAFFT_STRATEGY_PARTTREE;
	cfg.seqtype = MAFFT_SEQ_DNA;

	mafft_ctx_t *ctx = mafft_create(&cfg);
	if (!ctx) { ta->passed = 0; return NULL; }

	mafft_output_t *out = NULL;
	int rc = mafft_align(ctx, dna_names, dna_seqs, 3, &out, NULL);

	if (rc == MAFFT_OK && out && out->n_seqs == 3)
	{
		int i;
		ta->passed = 1;
		ta->n_seqs_out = out->n_seqs;
		ta->aligned_len_out = out->aligned_len;
		/* Verify all sequences have matching lengths */
		for (i = 0; i < out->n_seqs; i++)
		{
			if ((int)strlen(out->seqs[i]) != out->aligned_len)
				ta->passed = 0;
		}
	}
	else
		ta->passed = 0;

	if (out) mafft_output_free(out);
	mafft_destroy(ctx);
	return NULL;
}

static void test_concurrency(void)
{
	pthread_t t1, t2;
	thread_arg_t arg1 = {0, 0, 0};
	thread_arg_t arg2 = {0, 0, 0};

	pthread_create(&t1, NULL, align_thread, &arg1);
	pthread_create(&t2, NULL, align_thread, &arg2);
	pthread_join(t1, NULL);
	pthread_join(t2, NULL);

	CHECK(arg1.passed, "thread 1 completed with correct output");
	CHECK(arg2.passed, "thread 2 completed with correct output");
	CHECK(arg1.n_seqs_out == 3, "thread 1: 3 output sequences");
	CHECK(arg2.n_seqs_out == 3, "thread 2: 3 output sequences");
	CHECK(arg1.aligned_len_out > 0, "thread 1: positive aligned_len");
	CHECK(arg2.aligned_len_out > 0, "thread 2: positive aligned_len");
}

/* ---- Custom allocator test ---- */

static int alloc_count = 0;

static void *counting_alloc(size_t size, void *ud)
{
	(void)ud;
	alloc_count++;
	return malloc(size);
}

static void counting_free(void *ptr, void *ud)
{
	(void)ud;
	if (ptr) alloc_count--;
	free(ptr);
}

static void test_custom_allocator(void)
{
	mafft_config_t cfg;
	mafft_config_init(&cfg);
	cfg.strategy = MAFFT_STRATEGY_PARTTREE;
	cfg.seqtype = MAFFT_SEQ_DNA;
	cfg.alloc_fn = counting_alloc;
	cfg.free_fn = counting_free;

	alloc_count = 0;

	mafft_ctx_t *ctx = mafft_create(&cfg);
	mafft_output_t *out = NULL;
	int rc = mafft_align(ctx, dna_names, dna_seqs, 3, &out, NULL);

	CHECK(rc == MAFFT_OK, "custom allocator: align succeeds");
	CHECK(alloc_count == 1, "custom allocator: exactly 1 outstanding allocation");

	if (out) mafft_output_free(out);
	CHECK(alloc_count == 0, "custom allocator: 0 allocations after output_free");

	mafft_destroy(ctx);
}

/* ---- Main ---- */

int main(void)
{
	test_abi_size();
	test_config_defaults();
	test_create_destroy();
	test_half_paired_allocator();
	test_strerror();
	test_last_error_null();
	test_bad_input();
	test_parttree_align();
	test_seq_auto_detect();
	test_parttree_vs_native();
	test_protein_vs_native();
	test_concurrency();
	test_custom_allocator();

	fprintf(stdout, "\n%d passed, %d failed\n", n_pass, n_fail);
	if (n_fail > 0)
		fprintf(stdout, "FAIL\n");
	else
		fprintf(stdout, "PASS\n");

	return n_fail > 0 ? 1 : 0;
}
