#ifndef __GEIG_H
#define __GEIG_H

#include <Python.h>

#include <complex.h>

/*
 * LAPACK prototypes
 */

void dsygvd_(const int *itype,
             const char *jobz,
             const char *uplo,
             const int *N,
             double *A,
             const int *lda,
             double *B,
             const int *ldb,
             double *w,
             double *work,
             const int *lwork,
             int *iwork,
             const int *liwork,
             int *info);

void zhegvd_(const int *itype,
             const char *jobz,
             const char *uplo,
             const int *N,
             double complex *A,
             const int *lda,
             double complex *B,
             const int *ldb,
             double *w,
             double complex *work,
             const int *lwork,
             double *rwork,
             const int *lrwork,
             int *iwork,
             const int *liwork,
             int *info);


/*
 * Python interface
 */

PyObject *py_geig(PyObject *self, PyObject *args);

PyObject *py_free_geig_workspace(PyObject *self, PyObject *args);

#endif
