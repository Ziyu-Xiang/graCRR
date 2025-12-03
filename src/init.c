/* src/init.c */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* --- extern declarations --- */

/* bootstrap_individual.c */
extern SEXP bootstrap_individual(SEXP B_,
                                 SEXP epsilon_hat_,
                                 SEXP X_,
                                 SEXP W_hat_,
                                 SEXP S_hat_sqrt_inverse_,
                                 SEXP Z_,
                                 SEXP Param);

/* fabs_gra_bic.c */
extern SEXP LFabs(SEXP Y,
                  SEXP X,
                  SEXP Weight,
                  SEXP Epsilon,
                  SEXP Lam_min,
                  SEXP Xi,
                  SEXP Stoping,
                  SEXP Iter,
                  SEXP Param,
                  SEXP Edge,
                  SEXP EdgeRow,
                  SEXP EdgeCol,
                  SEXP Degree,
                  SEXP Max_S,
                  SEXP Lambda2,
                  SEXP H);

/* --- registration table --- */

static const R_CallMethodDef CallEntries[] = {
    {"bootstrap_individual", (DL_FUNC) &bootstrap_individual, 7},
    {"LFabs",                (DL_FUNC) &LFabs,                16},
    {NULL, NULL, 0}
};

/* --- init --- */

void R_init_graCRR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
