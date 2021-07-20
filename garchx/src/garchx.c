#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void GARCHXRECURSION(void *, void *, void *, void *, void *, void *);
extern void GARCHXRECURSIONSIM(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* describe routines */
static const R_CMethodDef CEntries[] = {
  {"GARCHXRECURSION", (DL_FUNC) &GARCHXRECURSION, 6},
  {"GARCHXRECURSIONSIM", (DL_FUNC) &GARCHXRECURSIONSIM, 12},
  {NULL, NULL, 0}
};

/* register routines */
void R_init_garchx(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
