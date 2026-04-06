/*
 * Deprecated: use mafft_api.h for new code.
 *
 * This header exposes the original low-level library entry points.
 * They remain functional but the mafft_api.h interface is preferred:
 * it provides structured config, opaque context, log capture, and
 * thread-safe access via an internal mutex.
 */
extern int disttbfast( int ngui, int lgui, char **namegui, char **seqgui, int argc, char **argv, int (*callback)(int, int, char*));
extern int tbfast_library( int ngui, int lgui, char **namegui, char **seqgui, int argc, char **argv, int (*callback)(int, int, char*));
extern int dvtditr_library( int ngui, int lgui, char **namegui, char **seqgui, int argc, char **argv, int (*callback)(int, int, char*));
extern int splittbfast_library( int ngui, int lgui, char **namegui, char **seqgui, int argc, char **argv, int (*callback)(int, int, char*));
extern const char *mafft_get_log(void);
extern void mafft_clear_log(void);
#define GUI_ERROR 1
#define GUI_LENGTHOVER 2
#define GUI_CANCEL 3
