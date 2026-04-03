#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include "mafft.h"

/*
 * Test that splittbfast_library() captures all log output internally
 * and does not leak anything to stderr or stdout.
 *
 * Strategy: redirect stderr to a pipe before calling the library.
 * After the call, drain the pipe — it must be empty.  Meanwhile,
 * mafft_get_log() must return a non-empty string containing the
 * expected diagnostic messages.
 */
int main(void)
{
	int i;
	int res;
	int n = 3;
	int l = 10000;
	int pass = 1;

	char **seq = (char **)calloc(n, sizeof(char *));
	char **name = (char **)calloc(n, sizeof(char *));
	for (i = 0; i < n; i++) seq[i] = calloc(l + 1, sizeof(char));
	for (i = 0; i < n; i++) name[i] = calloc(100, sizeof(char));

	strcpy(name[0], "s1");
	strcpy(name[1], "s2");
	strcpy(name[2], "s3");
	strcpy(seq[0], "ACGTACGTACGT");
	strcpy(seq[1], "ACGAACGTACGT");
	strcpy(seq[2], "ACGTACGAACGT");

	/* argv matching: mafft --quiet --parttree (DNA, default penalties) */
	int argc = 12;
	char **argv = (char **)calloc(argc, sizeof(char *));
	for (i = 0; i < argc; i++) argv[i] = calloc(100, sizeof(char));
	strcpy(argv[0], "splittbfast");
	strcpy(argv[1], "-D");          /* DNA mode */
	strcpy(argv[2], "-f");
	strcpy(argv[3], "-1.53");       /* gap open penalty */
	strcpy(argv[4], "-Q");
	strcpy(argv[5], "100");         /* spfactor */
	strcpy(argv[6], "-h");
	strcpy(argv[7], "0");           /* aof */
	strcpy(argv[8], "-p");
	strcpy(argv[9], "50");          /* partsize */
	strcpy(argv[10], "-s");
	strcpy(argv[11], "-1");         /* groupsize = njob+1 */

	/* --- set up a pipe to detect any stderr leakage --- */
	int pipefd[2];
	if (pipe(pipefd) != 0)
	{
		perror("pipe");
		return 1;
	}
	int orig_stderr = dup(STDERR_FILENO);
	dup2(pipefd[1], STDERR_FILENO);
	close(pipefd[1]);

	res = splittbfast_library(n, l, name, seq, argc, argv, NULL);

	/* restore stderr so our own fprintf works */
	fflush(stderr);
	dup2(orig_stderr, STDERR_FILENO);
	close(orig_stderr);

	/* drain the pipe — anything here is leakage */
	{
		char leak_buf[4096];
		/* set read end non-blocking so we don't hang */
		int flags = fcntl(pipefd[0], F_GETFL);
		fcntl(pipefd[0], F_SETFL, flags | O_NONBLOCK);
		ssize_t leaked = read(pipefd[0], leak_buf, sizeof(leak_buf) - 1);
		close(pipefd[0]);

		if (leaked > 0)
		{
			leak_buf[leaked] = '\0';
			fprintf(stderr, "FAIL: %zd bytes leaked to stderr:\n%s\n",
					leaked, leak_buf);
			pass = 0;
		}
	}

	/* --- check the library return code --- */
	if (res != 0)
	{
		fprintf(stderr, "FAIL: splittbfast_library returned %d\n", res);
		pass = 0;
	}

	/* --- verify alignment correctness --- */
	if (res == 0)
	{
		int len0 = strlen(seq[0]);
		for (i = 1; i < n; i++)
		{
			if ((int)strlen(seq[i]) != len0)
			{
				fprintf(stderr, "FAIL: seq[%d] length %d != seq[0] length %d\n",
						i, (int)strlen(seq[i]), len0);
				pass = 0;
			}
		}
	}

	/* --- verify log capture --- */
	{
		const char *log = mafft_get_log();
		if (log[0] == '\0')
		{
			fprintf(stderr, "FAIL: mafft_get_log() returned empty string\n");
			pass = 0;
		}
		else
		{
			/* spot-check for known messages */
			if (!strstr(log, "Gap Penalty"))
			{
				fprintf(stderr, "FAIL: log missing 'Gap Penalty' message\n");
				pass = 0;
			}
			if (!strstr(log, "Done."))
			{
				fprintf(stderr, "FAIL: log missing 'Done.' message\n");
				pass = 0;
			}
		}
		mafft_clear_log();

		/* after clear, should be empty */
		log = mafft_get_log();
		if (log[0] != '\0')
		{
			fprintf(stderr, "FAIL: mafft_clear_log() did not clear buffer\n");
			pass = 0;
		}
	}

	if (pass)
		fprintf(stdout, "PASS: alignment correct, log captured, no stderr leakage\n");

	for (i = 0; i < n; i++) free(seq[i]);
	free(seq);
	for (i = 0; i < n; i++) free(name[i]);
	free(name);
	for (i = 0; i < argc; i++) free(argv[i]);
	free(argv);

	return pass ? 0 : 1;
}
