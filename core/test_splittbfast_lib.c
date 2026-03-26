#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mafft.h"

int main(void)
{
	int i;
	int argc;
	char **argv;
	char **seq;
	char **name;
	int res;
	int n, l;

	n = 3;
	l = 10000;
	seq = (char **)calloc(n, sizeof(char *));
	name = (char **)calloc(n, sizeof(char *));
	for (i = 0; i < n; i++) seq[i] = calloc(l + 1, sizeof(char));
	for (i = 0; i < n; i++) name[i] = calloc(100, sizeof(char));

	strcpy(name[0], "s1");
	strcpy(name[1], "s2");
	strcpy(name[2], "s3");

	strcpy(seq[0], "ACGTACGTACGT");
	strcpy(seq[1], "ACGAACGTACGT");
	strcpy(seq[2], "ACGTACGAACGT");

	/* Construct argv matching: mafft --quiet --parttree defaults */
	argc = 13;
	argv = (char **)calloc(argc, sizeof(char *));
	for (i = 0; i < argc; i++) argv[i] = calloc(100, sizeof(char));
	strcpy(argv[0], "splittbfast");
	strcpy(argv[1], "-D");          /* DNA */
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
	strcpy(argv[12], "-X");
	/* -X takes the next arg as treealg param */
	/* actually, let me re-read the arguments parsing... */

	/* Simpler: just use minimal required args */
	/* Reset and use fewer args */
	for (i = 0; i < argc; i++) free(argv[i]);
	free(argv);

	argc = 11;
	argv = (char **)calloc(argc, sizeof(char *));
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
	/* -s takes next arg, but we have no more args... */
	/* Need to add the groupsize value */

	for (i = 0; i < argc; i++) free(argv[i]);
	free(argv);

	argc = 12;
	argv = (char **)calloc(argc, sizeof(char *));
	for (i = 0; i < argc; i++) argv[i] = calloc(100, sizeof(char));
	strcpy(argv[0], "splittbfast");
	strcpy(argv[1], "-D");
	strcpy(argv[2], "-f");
	strcpy(argv[3], "-1.53");
	strcpy(argv[4], "-Q");
	strcpy(argv[5], "100");
	strcpy(argv[6], "-h");
	strcpy(argv[7], "0");
	strcpy(argv[8], "-p");
	strcpy(argv[9], "50");
	strcpy(argv[10], "-s");
	strcpy(argv[11], "-1");

	fprintf(stderr, "Calling splittbfast_library...\n");
	res = splittbfast_library(n, l, name, seq, argc, argv, NULL);
	fprintf(stderr, "Result: %d\n", res);

	if (res == 0)
	{
		fprintf(stdout, "Aligned output:\n");
		for (i = 0; i < n; i++)
			fprintf(stdout, ">%s\n%s\n", name[i], seq[i]);

		/* Verify all sequences have the same length */
		int len0 = strlen(seq[0]);
		int all_same = 1;
		for (i = 1; i < n; i++)
		{
			if ((int)strlen(seq[i]) != len0)
			{
				fprintf(stderr, "ERROR: seq[%d] length %d != seq[0] length %d\n",
						i, (int)strlen(seq[i]), len0);
				all_same = 0;
			}
		}
		if (all_same)
			fprintf(stdout, "OK: All aligned sequences have length %d\n", len0);
	}
	else
	{
		fprintf(stderr, "FAILED: splittbfast_library returned %d\n", res);
	}

	for (i = 0; i < n; i++) free(seq[i]);
	free(seq);
	for (i = 0; i < n; i++) free(name[i]);
	free(name);
	for (i = 0; i < argc; i++) free(argv[i]);
	free(argv);

	return res;
}
