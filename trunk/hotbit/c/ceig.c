#include <Python.h>
#include <arrayobject.h>
#include <stdlib.h>
#include <complex.h>

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


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

double *ceigr_(int N, double *A, double *B) {
    int *iwork;
    int itype = 1;
    int liwork = 3+5*N;
    int lwork = 1+6*N+2*N*N;
    int info, i, j, m;
    char jobz = 'V';
    char uplo = 'L';
    double *w, *work, k;

    w = (double*)malloc(N * sizeof (double));
    work = (double*)malloc(lwork * sizeof (double));
    iwork = (int*)malloc(liwork * sizeof (int));

    // A and B are real and Hermitian, so no need to transpose.
    dsygvd_(&itype, &jobz, &uplo, &N, A, &N, B, &N, w, work, &lwork, iwork, &liwork, &info);

    free(work);
    free(iwork);

    // The resulting matrix of eigenvectors must be transposed, because
    // we need the C-ordered matrix to return.
    for (i=0; i<N; i++) {
        for (j=i+1; j<N; j++) {
            k = A[i*N+j];
            A[i*N+j] = A[j*N+i];
            A[j*N+i] = k;
        }
    }

    if (info != 0) {
        m = min(N, 3);
        printf("***********************************************************\n");
        printf("Parameter INFO from dsygvd: %i\n",info);
        printf("matrix size N=%i\n",N);
        if (info < 0) {
            printf("The info-th argument for dsygvd had an illegal value.\n");
        }
        else if (info > 0 && info <= N ) {
            printf("The algorithm failed to compute an eigenvalue while working\n");
            printf("on the submatrix lying in rows and columns INFO/(N+1) through mod(INFO,N+1)\n");
        }
        else if (info > N) {
            printf("If INFO = N + i, for 1 <= i <= N, then the leading minor of order i of \n");
            printf("S is not positive definite. The factorization of S could not be completed and\n");
            printf("no eigenvalues or eigenvectors were computed. (S=overlap matrix)\n");
        }
        printf("***********************************************************\n");
        printf("H=\n");
        for (i=0; i<m; i++) {
            for (j=0; j<m; j++) {
                printf("%f ", A[i*N+j]);
            }
            printf("\n");
        }
        printf("... and so on.\n");
        printf("S=\n");
        for (i=0; i<m; i++) {
            for (j=0; j<m; j++) {
                printf("%f ", B[i*N+j]);
            }
            printf("\n");
        }
        printf("... and so on.\n");
        //stop 'Error in LAPACK diagonalization routine.'
    }
    return w;
}

static PyObject *ceigr(PyObject *self, PyObject *args) {
    int N;
    npy_intp *A_dim, *ev_dim;
    double *A, *B, *ev;
    PyObject *po_A, *po_B;
    PyArrayObject *pao_eV, *pao_A;

    if (!PyArg_ParseTuple(args, "OO", &po_A, &po_B))
        return NULL;

    A = (double *)(((PyArrayObject *)po_A)->data);
    B = (double *)(((PyArrayObject *)po_B)->data);
    A_dim = (npy_intp *)(((PyArrayObject *)po_A)->dimensions);

    N = A_dim[0];
    ev_dim = (npy_intp *)malloc(N *sizeof (int));
    ev_dim[0] = N;

    ev = ceigr_(N, A, B);

    pao_eV = (PyArrayObject *)PyArray_SimpleNewFromData(1, ev_dim, NPY_DOUBLE, ev);
    pao_A = (PyArrayObject *)PyArray_SimpleNewFromData(2, A_dim, NPY_DOUBLE, A);
    return Py_BuildValue("OO", pao_eV, pao_A);
}

double *ceigc_(int N, double complex *A, double complex *B) {
    //! complex version for generalized eigenvalue problem
    int lwork, liwork, lrwork;
    int info, i, j, m;
    int *iwork;
    int itype = 1;
    char jobz = 'V';
    char uplo = 'L';
    double *w, *rwork;
    double complex *work, tmp;

    liwork = 3 + 5*N;
    lwork  = 2*N + N*N;
    lrwork = 1 + 5*N + 2*N*N;

    w = (double*)malloc(N * sizeof (double));
    rwork = (double*)malloc(lrwork * sizeof (double));
    work = (double complex*)malloc(lwork * sizeof (double complex));
    iwork = (int*)malloc(liwork * sizeof (int));

    // We need to transpose the A and B matrices because zhegvd expects
    // fortran-ordered arrays. However, A and B are Hermitian, so
    // complex conjugate is enough.
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            A[i*N+j] = conj(A[i*N+j]);
        }
    }
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            B[i*N+j] = conj(B[i*N+j]);
        }
    }

    zhegvd_(&itype, &jobz, &uplo, &N, A, &N, B, &N, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);

    // The resulting matrix of eigenvectors must be transposed, because
    // we need the C-ordered matrix to return.
    for (i=0; i<N; i++) {
        for (j=i+1; j<N; j++) {
            tmp = A[i*N+j];
            A[i*N+j] = A[j*N+i];
            A[j*N+i] = tmp;
        }
    }
    free(work);
    free(iwork);
    free(rwork);

    if (info != 0) {
        m = min(N, 4);
        printf("***********************************************************");
        printf("Parameter INFO from dsygvd: %i\n",info);
        printf("matrix size N=%i\n",N);
        if ( info<0 ) {
            printf("The info-th argument for dsygvd had an illegal value.\n");
        }
        else if (info > 0 && info <= N) {
            printf("The algorithm failed to compute an eigenvalue while working\n");
            printf("on the submatrix lying in rows and columns INFO/(N+1) through mod(INFO,N+1)\n");
        }
        else if (info > N) {
            printf("If INFO = N + i, for 1 <= i <= N, then the leading minor of order i of \n");
            printf("S is not positive definite. The factorization of S could not be completed and\n");
            printf("no eigenvalues or eigenvectors were computed. (S=overlap matrix)\n");
        }
        printf("***********************************************************\n");
        printf("H = (complex numbers = number pairs) \n");
        for (i=0; i < m; i++) {
            for (j=0; j < m; j++) {
                printf("%0.4f %0.4fi", creal(A[i*N+j]), cimag(A[i*N+j]));
            }
        }
        printf("... and so on.\n");
        printf("S = (complex numbers = number pairs)\n");
        for (i=0; i < m; i++) {
            for (j=0; j < m; j++) {
                printf("%0.4f %0.4fi", creal(B[i*N+j]), cimag(B[i*N+j]));
            }
        }
        printf("... and so on.\n");
    }
    return w;
}

static PyObject *ceigc(PyObject *self, PyObject *args) {
    int N;
    npy_intp *A_dim, *ev_dim;
    double *ev;
    double complex *A, *B;
    PyObject *po_A, *po_B;
    PyArrayObject *pao_ev, *pao_A, *pao_Ar, *pao_Ai;

    if (!PyArg_ParseTuple(args, "OO", &po_A, &po_B))
        return NULL;

    A = (double complex*)(((PyArrayObject *)po_A)->data);
    B = (double complex*)(((PyArrayObject *)po_B)->data);
    A_dim = (npy_intp *)(((PyArrayObject *)po_A)->dimensions);

    N = A_dim[0];
    ev_dim = (npy_intp *)malloc(N *sizeof (int));
    ev_dim[0] = N;

    ev = ceigc_(N, A, B);

    pao_A = (PyArrayObject *)PyArray_SimpleNewFromData(2, A_dim, NPY_CDOUBLE, A);
    pao_ev = (PyArrayObject *)PyArray_SimpleNewFromData(1, ev_dim, NPY_DOUBLE, ev);
    return Py_BuildValue("OO", pao_ev, pao_A);
}

PyMethodDef methods[] = {
    {"ceigr", ceigr, METH_VARARGS, "Solve the generalized eigenvalue equation"},
    {"ceigc", ceigc, METH_VARARGS, "Solve the generalized (complex) eigenvalue equation"},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initceig() {
    (void) Py_InitModule("ceig", methods);
    import_array();
}

