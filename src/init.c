#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void isMon(void *, void *, void *, void *);
extern void permNext(void *, void *);
extern void weightedMean(void *, void *, void *, void *, void *, void *, void *, void *);
extern void wmonreg(void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"isMon",        (DL_FUNC) &isMon,        4},
  {"permNext",     (DL_FUNC) &permNext,     2},
  {"weightedMean", (DL_FUNC) &weightedMean, 8},
  {"wmonreg",      (DL_FUNC) &wmonreg,      3},
  {NULL, NULL, 0}
};

void R_init_smacof(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
