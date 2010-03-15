#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL HOTBIT_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include "geig.h"


static int liwork = 0;
static int lrwork = 0;
static int lcwork = 0;
static int *iwork = NULL;
static double *rwork = NULL;
static double complex *cwork = NULL;


/*
 * Call dsygvd, allocate the appropriate workspace,
 * and allocate the appropriate return buffer.
 */
double *
geigr(int N, double *A, double *B) {
    int itype = 1;
    int info;
    char jobz = 'V';
    char uplo = 'L';
    double *w;


    /*
     * Check if work buffers are large enough for this eigenvalue problem.
     * Otherwise, reallocate.
     */
  
    if (lrwork < 1+6*N+2*N*N) {
        if (rwork)
            free(rwork);

        lrwork = 1+6*N+2*N*N;
        rwork = (double*) malloc(lrwork * sizeof(double));

        if (!rwork) {
            PyErr_SetString(PyExc_RuntimeError, "Error in eigensolver while allocating \"rwork\".");
            return NULL;
        }
    }

    if (liwork < 3+5*N) {
        if (iwork)
            free(iwork);

        liwork = 3+5*N;
        iwork = (int*) malloc(liwork * sizeof(int));

        if (!iwork) {
            PyErr_SetString(PyExc_RuntimeError, "Error in eigensolver while allocating \"iwork\".");
            return NULL;
        }
    }


    /*
     * Allocate eigenvalue array
     */

    w = (double*) malloc(N * sizeof (double));
    if (!w) {
        PyErr_SetString(PyExc_RuntimeError, "Error in eigensolver while allocating \"w\".");
        return NULL;
    }

    dsygvd_(&itype, &jobz, &uplo, &N, A, &N, B, &N, 
            w, rwork, &lrwork, iwork, &liwork,
            &info);

    if (info) {
        free(w);
        return NULL;
    }

    return w;
}



/*
 * Call zhygvd, allocate the appropriate workspace,
 * and allocate the appropriate return buffer.
 */
double *
geigc(int N, double complex *A, double complex *B) {
    int itype = 1;
    int info;
    char jobz = 'V';
    char uplo = 'L';
    double *w;


    /*
     * Check if work buffers are large enough for this eigenvalue problem.
     * Otherwise, reallocate.
     */

    if (lcwork < 2*N+N*N) {
        if (cwork)
            free(cwork);

        lcwork = 2*N+N*N;
        cwork = (double complex*) malloc(lcwork * sizeof(double complex));

        if (!cwork) {
            PyErr_SetString(PyExc_RuntimeError, "Error in eigensolver while allocating \"cwork\".");
            return NULL;
        }
    }

    if (lrwork < 1+5*N+2*N*N) {
        if (rwork)
            free(rwork);

        lrwork = 1+5*N+2*N*N;
        rwork = (double*) malloc(lrwork * sizeof(double));

        if (!rwork) {
            PyErr_SetString(PyExc_RuntimeError, "Error in eigensolver while allocating \"rwork\".");
            return NULL;
        }
    }

    if (liwork < 3+5*N) {
        if (iwork)
            free(iwork);

        liwork = 3+5*N;
        iwork = (int*) malloc(liwork * sizeof(int));

        if (!iwork) {
            PyErr_SetString(PyExc_RuntimeError, "Error in eigensolver while allocating \"iwork\".");
            return NULL;
        }
    }


    /*
     * Allocate eigenvalue array
     */

    w = (double*) malloc(N * sizeof (double));
    if (!w) {
        PyErr_SetString(PyExc_RuntimeError, "Error in eigensolver while allocating \"w\".");
        return NULL;
    }

    zhegvd_(&itype, &jobz, &uplo, &N, A, &N, B, &N,
            w, cwork, &lcwork, rwork, &lrwork, iwork, &liwork, 
            &info);

    if (info) {
        free(w);
        return NULL;
    }

    return w;
}


PyObject *
py_geig(PyObject *self, PyObject *args) {
    int N;
    npy_intp *A_dims, *B_dims;
    npy_intp ev_dims[1];
    double *ev;
    PyObject *po_A, *po_B;
    PyObject *po_eV;

    int typenum;
    PyObject *po_for_A, *po_for_B;

    if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &po_A, &PyArray_Type, &po_B))
        return NULL;

    A_dims = PyArray_DIMS(po_A);
    B_dims = PyArray_DIMS(po_B);

    if (A_dims[0] != B_dims[0] || A_dims[1] != B_dims[1] || A_dims[0] != A_dims[1]) {
        PyErr_SetString(PyExc_TypeError,
                        "The sizes of matrix A and matrix B need to be identical, "
                        "and the matrices need to be square.");
        return NULL;
    }

    typenum = NPY_DOUBLE;
    if (PyArray_ISCOMPLEX(po_A) || PyArray_ISCOMPLEX(po_B)) {
        typenum = NPY_CDOUBLE;
    }

    N = A_dims[0];

    /*
     * Create a Fortran-ordered A and B matrix.
     * The A matrix will additionally contain the eigenvectors on return.
     */
    po_for_A = PyArray_FROMANY(po_A, typenum, 2, 2, NPY_INOUT_FARRAY);
    if (!po_for_A)
        return NULL;
    po_for_B = PyArray_FROMANY(po_B, typenum, 2, 2, NPY_IN_FARRAY);
    if (!po_for_B)
        return NULL;

    /*
     * Solve the eigenvalue problem
     */
    if (typenum == NPY_DOUBLE) {
        ev = geigr(N, PyArray_DATA(po_for_A),  PyArray_DATA(po_for_B));
    } else {
        ev = geigc(N, PyArray_DATA(po_for_A),  PyArray_DATA(po_for_B));
    }
    if (!ev)
        return NULL;

    /*
     * Create a Python array containing the eigenvectors.
     * Note that we must use PyArray_New, otherwise the data won't be freed
     * when the array object goes out of scope.
     */
    ev_dims[0] = N;
    po_eV = PyArray_New(&PyArray_Type, 1, ev_dims, NPY_DOUBLE, NULL, ev, NPY_OWNDATA, 0, NULL);

    /* Release the Fortran-ordered B matrix */
    Py_DECREF(po_for_B);

    return Py_BuildValue("OO", po_eV, po_for_A);
}


/*
 * Free workspace arrays
 */
PyObject *
py_free_geig_workspace(PyObject *self, PyObject *args)
{
    if (iwork)
        free(iwork);
    if (rwork)
        free(rwork);
    if (cwork)
        free(cwork);

    liwork = 0;
    lrwork = 0;
    lcwork = 0;

    Py_RETURN_NONE;
}

