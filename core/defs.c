#include "mltaln.h"
#include "dp.h"

TLS int commonAlloc1 = 0;
TLS int commonAlloc2 = 0;
TLS int **commonIP = NULL;
TLS int **commonJP = NULL;
int nthread = 1;
int nthreadpair = 1;
int randomseed = 0;
int parallelizationstrategy = BAATARI1;


char modelname[500];
int njob, nlenmax;
int amino_n[0x100];
char amino_grp[0x100];
//int amino_dis[0x100][0x100];
int **amino_dis = NULL;
double **n_disLN = NULL;
//double amino_dis_consweight_multi[0x100][0x100];
double **amino_dis_consweight_multi = NULL;
int **n_dis = NULL;
int **n_disFFT = NULL;
double **n_dis_consweight_multi = NULL;
unsigned char amino[0x100];
double polarity[0x100];
double volume[0x100];
int ribosumdis[37][37];

int ppid;
double thrinter;
double fastathreshold;
int pslocal, ppslocal;
int constraint;
int divpairscore;
int fmodel; // 1-> fmodel 0->default -1->raw
int nblosum; // 45, 50, 62, 80
int kobetsubunkatsu;
int bunkatsu;
int dorp = NOTSPECIFIED; // arguments de shitei suruto, tbfast -> pairlocalalign no yobidashi de futsugou
int niter;
int contin;
int calledByXced;
int devide;
int scmtd;
int weight;
int utree;
int tbutree;
int refine;
int check;
double cut;
int cooling;
int trywarp = 0;
int penalty, ppenalty, penaltyLN;
int penalty_dist, ppenalty_dist;
int RNApenalty, RNAppenalty;
int RNApenalty_ex, RNAppenalty_ex;
int penalty_ex, ppenalty_ex, penalty_exLN;
int penalty_EX, ppenalty_EX;
int penalty_OP, ppenalty_OP;
int penalty_shift, ppenalty_shift;
double penalty_shift_factor = 1000.0;
int RNAthr, RNApthr;
int offset, poffset, offsetLN, offsetFFT;
int scoremtx;
int TMorJTT;
char use_fft;
char force_fft;
int nevermemsave;
int fftscore;
int fftWinSize;
int fftThreshold;
int fftRepeatStop;
int fftNoAnchStop;
int divWinSize;
int divThreshold;
int disp;
int outgap = 1;
char alg;
int cnst;
int mix;
int tbitr;
int tbweight;
int tbrweight;
int disopt;
int pamN;
int checkC;
double geta2;
int treemethod;
int kimuraR;
char *swopt;
int fftkeika;
int score_check;
int makedistmtx;
char *inputfile;
char *addfile;
int addprofile = 1;
int rnakozo;
char rnaprediction;
int scoreout = 0;
int spscoreout = 0;
int outnumber = 0;
int legacygapcost = 0;
double minimumweight = 0.0005;
int nwildcard = 0;

char *signalSM;
FILE *prep_g;
FILE *trap_g;
char **seq_g;
char **res_g;

double consweight_multi = 1.0;
double consweight_rna = 0.0;
char RNAscoremtx = 'n';

TLS char *newgapstr = "-";

int nalphabets = 26;
int nscoredalphabets = 20;

double specificityconsideration = 0.0;
int ndistclass = 10;
int maxdistclass = -1;

int gmsg = 0;
int mafft_library_mode = 0;

double sueff_global = SUEFF;

double lenfaca, lenfacb, lenfacc, lenfacd;
int maxl, tsize;

/* usetmpfile globals (declared here so initglobalvariables can reset them) */
int compacttree = 0;
int lhlimit = INT_MAX;
int specifictarget = 0;
int nadd = 0;
int usenaivescoreinsteadofalignmentscore = 0;
int nthreadreadlh = 1;


void initglobalvariables()
{
	/* TLS DP arrays */
	commonAlloc1 = 0;
	commonAlloc2 = 0;
	commonIP = NULL;
	commonJP = NULL;

	/* Threading */
	nthread = 1;
	nthreadpair = 1;
	randomseed = 0;
	parallelizationstrategy = BAATARI1;

	/* Sequence identity */
	modelname[0] = '\0';
	njob = 0;
	nlenmax = 0;

	/* Scoring matrices -- freed by freeconstants(), NULLed here.
	 * Static lookup arrays (amino_n, amino_grp, amino, polarity,
	 * volume, ribosumdis) are intentionally NOT reset: they are
	 * unconditionally overwritten by makematrix()/constants() on
	 * each alignment call before use. */
	amino_dis = NULL;
	n_disLN = NULL;
	amino_dis_consweight_multi = NULL;
	n_dis = NULL;
	n_disFFT = NULL;
	n_dis_consweight_multi = NULL;

	/* Miscellaneous parameters */
	ppid = 0;
	thrinter = 0.0;
	fastathreshold = 0.0;
	pslocal = 0;
	ppslocal = 0;
	constraint = 0;
	divpairscore = 0;
	fmodel = 0;
	nblosum = 62;
	kobetsubunkatsu = 0;
	bunkatsu = 0;
	dorp = NOTSPECIFIED;
	niter = 0;
	contin = 0;
	calledByXced = 0;
	devide = 0;
	scmtd = 5;
	weight = 3;
	utree = 1;
	tbutree = 1;
	refine = 0;
	check = 1;
	cut = 0.0;
	cooling = 0;
	trywarp = 0;

	/* Penalties */
	penalty = 0;
	ppenalty = -1530;
	penaltyLN = 0;
	penalty_dist = 0;
	ppenalty_dist = 0;
	RNApenalty = 0;
	RNAppenalty = 0;
	RNApenalty_ex = 0;
	RNAppenalty_ex = 0;
	penalty_ex = 0;
	ppenalty_ex = NOTSPECIFIED;
	penalty_exLN = 0;
	penalty_EX = 0;
	ppenalty_EX = 0;
	penalty_OP = 0;
	ppenalty_OP = 0;
	penalty_shift = 0;
	ppenalty_shift = 0;
	penalty_shift_factor = 1000.0;
	RNAthr = 0;
	RNApthr = 0;
	offset = 0;
	poffset = -123;
	offsetLN = 0;
	offsetFFT = 0;

	/* Scoring / algorithm selection */
	scoremtx = 1;
	TMorJTT = JTT;
	use_fft = 0;
	force_fft = 0;
	nevermemsave = 0;
	fftscore = 1;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
	divWinSize = 0;
	divThreshold = 0;
	disp = 0;
	outgap = 1;
	alg = 'A';
	cnst = 0;
	mix = 0;
	tbitr = 0;
	tbweight = 0;
	tbrweight = 3;
	disopt = 0;
	pamN = NOTSPECIFIED;
	checkC = 0;
	geta2 = GETA2;
	treemethod = 'X';
	kimuraR = NOTSPECIFIED;
	swopt = NULL;
	fftkeika = 0;
	score_check = 0;
	makedistmtx = 0;

	/* I/O */
	inputfile = NULL;
	addfile = NULL;
	addprofile = 1;
	rnakozo = 0;
	rnaprediction = 0;
	scoreout = 0;
	spscoreout = 0;
	outnumber = 0;
	legacygapcost = 0;
	minimumweight = 0.0005;
	nwildcard = 0;

	signalSM = NULL;
	prep_g = NULL;
	trap_g = NULL;
	seq_g = NULL;
	res_g = NULL;

	/* Weights */
	consweight_multi = 1.0;
	consweight_rna = 0.0;
	RNAscoremtx = 'n';

	newgapstr = "-";

	nalphabets = 26;
	nscoredalphabets = 20;

	specificityconsideration = 0.0;
	ndistclass = 10;
	maxdistclass = -1;

	gmsg = 0;

	sueff_global = SUEFF;

	lenfaca = 0.0;
	lenfacb = 0.0;
	lenfacc = 0.0;
	lenfacd = 0.0;
	maxl = 0;
	tsize = 0;

	/* usetmpfile globals */
	compacttree = 0;
	lhlimit = INT_MAX;
	specifictarget = 0;
	nadd = 0;
	usenaivescoreinsteadofalignmentscore = 0;
	nthreadreadlh = 1;
}
