#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _TCGP_fun_dev(SEXP, SEXP);
extern SEXP _TCGP_fun_mul(SEXP, SEXP);
extern SEXP _TCGP_G_thres_cpp(SEXP, SEXP);
extern SEXP _TCGP_loglikelihood_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TCGP_post_e_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TCGP_post_tau_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TCGP_rcpparma_bothproducts(SEXP);
extern SEXP _TCGP_rcpparma_hello_world();
extern SEXP _TCGP_rcpparma_innerproduct(SEXP);
extern SEXP _TCGP_rcpparma_outerproduct(SEXP);
extern SEXP _TCGP_sample_c_l_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TCGP_TCGP_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _TCGP_sample_thres_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_TCGP_fun_dev",               (DL_FUNC) &_TCGP_fun_dev,                2},
    {"_TCGP_fun_mul",               (DL_FUNC) &_TCGP_fun_mul,                2},
    {"_TCGP_G_thres_cpp",           (DL_FUNC) &_TCGP_G_thres_cpp,            2},
    {"_TCGP_loglikelihood_cpp",     (DL_FUNC) &_TCGP_loglikelihood_cpp,     10},
    {"_TCGP_post_e_cpp",            (DL_FUNC) &_TCGP_post_e_cpp,            13},
    {"_TCGP_post_tau_cpp",          (DL_FUNC) &_TCGP_post_tau_cpp,          12},
    {"_TCGP_rcpparma_bothproducts", (DL_FUNC) &_TCGP_rcpparma_bothproducts,  1},
    {"_TCGP_rcpparma_hello_world",  (DL_FUNC) &_TCGP_rcpparma_hello_world,   0},
    {"_TCGP_rcpparma_innerproduct", (DL_FUNC) &_TCGP_rcpparma_innerproduct,  1},
    {"_TCGP_rcpparma_outerproduct", (DL_FUNC) &_TCGP_rcpparma_outerproduct,  1},
    {"_TCGP_sample_c_l_cpp",        (DL_FUNC) &_TCGP_sample_c_l_cpp,        14},
    {"_TCGP_TCGP_fit",      (DL_FUNC) &_TCGP_TCGP_fit,      20},
    {"_TCGP_sample_thres_cpp",      (DL_FUNC) &_TCGP_sample_thres_cpp,      12},
    {NULL, NULL, 0}
};
void R_init_TCGP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
