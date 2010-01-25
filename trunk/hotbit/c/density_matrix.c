#include <Python.h>
#include <arrayobject.h>
#include <stdlib.h>
#include <complex.h>

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

double complex *c_rhoc_(double complex *wf, double *occ, int norb, int nk) {
    int i, j, k, n, mx;
    double complex *rho;

    rho = (double complex*)calloc(nk*norb*norb, sizeof(double complex));

    mx=-1;
    for (i=norb-1; i>=0; i--) {
        if (mx != -1) {
            break;
        }
        for (k=0; k<nk; k++) {
            if (occ[k*norb+i] > 1.0e-15) {
                mx = i;
                break;
            }
        }
    }
    // FIXME probably faster with different order of loops
    for (k=0; k<nk; k++) {
        for (i=0; i<norb; i++) {
            for (j=i; j<norb; j++) {
                for (n=0; n<mx+1; n++) {
                    rho[k*norb*norb+i*norb+j] += occ[k*norb+n]*wf[k*norb*norb+n*norb+i]*conj(wf[k*norb*norb+n*norb+j]);
                }
                if (i != j) {
                    rho[k*norb*norb+j*norb+i] = conj(rho[k*norb*norb+i*norb+j]);
                }
            }
        }
    }
    return rho;
}

static PyObject *c_rhoc(PyObject *self, PyObject *args) {
    int norb, nk;
    npy_intp *wf_dim;
    double *occ;
    double complex *wf, *rho;
    PyObject *po_wf, *po_occ;
    PyArrayObject *pao_rho;

    if (!PyArg_ParseTuple(args, "OO", &po_wf, &po_occ))
        return NULL;

    wf = (double complex*)(((PyArrayObject *)po_wf)->data);
    occ = (double *)(((PyArrayObject *)po_occ)->data);
    wf_dim = (npy_intp *)(((PyArrayObject *)po_wf)->dimensions);

    nk = wf_dim[0];
    norb = wf_dim[1];
    rho = c_rhoc_(wf, occ, norb, nk);

    pao_rho = (PyArrayObject *)PyArray_SimpleNewFromData(3, wf_dim, NPY_CDOUBLE, rho);
    return Py_BuildValue("O", pao_rho);
}

double complex *c_rhoec_(double complex *wf, double *occ, double *e, int norb, int nk) {
    int mx, i, j, k, n;
    double complex *rho;

    rho = (double complex*)calloc(nk*norb*norb, sizeof(double complex));

    mx=-1;
    for (i=norb-1; i>=0; i--) {
        if (mx != -1) {
            break;
        }
        for (k=0; k<nk; k++) {
            if (occ[k*norb+i] > 1.0e-15) {
                mx = i;
                break;
            }
        }
    }
    // FIXME probably faster with different order of loops
    for (k=0; k<nk; k++) {
        for (i=0; i<norb; i++) {
            for (j=i; j<norb; j++) {
                for (n=0; n<mx+1; n++) {
                    rho[k*norb*norb +i*norb +j] += e[k*norb +n]*occ[k*norb +n]*wf[k*norb*norb +n*norb +i]*conj(wf[k*norb*norb +n*norb +j]);
                    rho[k*norb*norb+j*norb+i] = conj(rho[k*norb*norb+i*norb+j]);
                }
            }
        }
    }
    return rho;
}

static PyObject *c_rhoec(PyObject *self, PyObject *args) {
    int norb, nk;
    npy_intp *wf_dim;
    double *occ, *e;
    double complex *wf, *rhoe;
    PyObject *po_wf, *po_occ, *po_e;
    PyArrayObject *pao_rhoe;

    if (!PyArg_ParseTuple(args, "OOO", &po_wf, &po_occ, &po_e))
        return NULL;

    wf = (double complex*)(((PyArrayObject *)po_wf)->data);
    occ = (double *)(((PyArrayObject *)po_occ)->data);
    e = (double *)(((PyArrayObject *)po_e)->data);
    wf_dim = (npy_intp *)(((PyArrayObject *)po_wf)->dimensions);

    nk = wf_dim[0];
    norb = wf_dim[1];
    rhoe = c_rhoec_(wf, occ, e, norb, nk);

    pao_rhoe = (PyArrayObject *)PyArray_SimpleNewFromData(3, wf_dim, NPY_CDOUBLE, rhoe);
    return Py_BuildValue("O", pao_rhoe);
}

PyMethodDef methods[] = {
    {"c_rhoc", c_rhoc, METH_VARARGS, "Calculate the complex density matrix."},
    {"c_rhoec", c_rhoec, METH_VARARGS, "Calculate the energy-weighted complex density matrix."},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initdensity_matrix() {
    (void) Py_InitModule("density_matrix", methods);
    import_array();
}

