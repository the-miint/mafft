#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mafft.h"

int main(void)
{
	int i;
	int res;
	int n = 3;
	int l = 10000;

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

	fprintf(stderr, "Calling splittbfast_library...\n");
	res = splittbfast_library(n, l, name, seq, argc, argv, NULL);

	if (res == 0)
	{
		/* Verify all sequences have the same aligned length */
		int len0 = strlen(seq[0]);
		int ok = 1;
		for (i = 1; i < n; i++)
		{
			if ((int)strlen(seq[i]) != len0)
			{
				fprintf(stderr, "FAIL: seq[%d] length %d != seq[0] length %d\n",
						i, (int)strlen(seq[i]), len0);
				ok = 0;
			}
		}
		if (ok)
			fprintf(stdout, "PASS: All aligned sequences have length %d\n", len0);
	}
	else
	{
		fprintf(stderr, "FAIL: splittbfast_library returned %d\n", res);
	}

	for (i = 0; i < n; i++) free(seq[i]);
	free(seq);
	for (i = 0; i < n; i++) free(name[i]);
	free(name);
	for (i = 0; i < argc; i++) free(argv[i]);
	free(argv);

	return res;
}
