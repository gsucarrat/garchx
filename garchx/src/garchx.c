#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void GARCHXRECURSION(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"GARCHXRECURSION", (DL_FUNC) &GARCHXRECURSION, 6},
  {NULL, NULL, 0}
};

void R_init_garchx(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
