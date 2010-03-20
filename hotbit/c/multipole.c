#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL HOTBIT_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include <complex.h>

#include "spherical.h"

#include "multipole.h"

#define min(a,b)        ((a)>(b)?(b):(a))
#define max(a,b)        ((a)<(b)?(b):(a))
#define lm2index(l, m)  ((l)*(l-1)/2 + (m) - 1)

/*
 * Helper module for the fast-multipole solvers.
 *
 * Multipole-to-multipole, multipole-to-local
 * and local-to-local transformations.
 */


/*
 * Shift the multipole moments of the child by dr and to the current multipole moments.
 */
void multipole_to_multipole(double *dr,                      /* shape [3] */
                            int l_max,
                            double *Ml0_of_child,            /* shape [0:l_max] */
                            double complex *Mlm_of_child,    /* shape [lm2index(l_max, l_max)] */
                            double *Ml0,                     /* shape [0:l_max] */
                            double complex *Mlm              /* shape [lm2index(l_max, l_max)] */
                            )
{
    double x, costh, phi;
    double Rl0[l_max+1];
    double complex Rlm[lm2index(l_max, l_max)+1];
    double complex cur_M, cur_R;

    int j, l, m, lambda, mu, sign;


    cartesian2spherical(dr, &x, &costh, &phi);

    solid_harmonic_R(x, costh, phi, l_max, Rl0, Rlm);

    j = 0;
    for (l = 0; l <= l_max; ++l) {

        /*
         * m = 0
         */

        for (lambda = 0; lambda <= l; ++lambda) {
            Ml0[l] = Ml0[l]
                + Ml0_of_child[l-lambda]*Rl0[lambda];

            sign = 1;
            for(mu = 1; mu <= min(lambda, -lambda+l); ++mu) {
                sign = sign*(-1);

                Ml0[l] = Ml0[l]
                    + 2*sign*creal(Mlm_of_child[lm2index(l-lambda, mu)]*Rlm[lm2index(lambda, mu)]);
            }
        }

        /*
         * m != 0
         */

        for (m = 1; m <= l; ++m) {
            for (lambda = 0; lambda <= l; ++lambda) {
                for (mu = max(-lambda, lambda-l+m); mu <= min(lambda, -lambda+l+m); ++mu) {

                    if (m-mu < 0) {
                        cur_M = pow(-1, mu-m) * conj(Mlm_of_child[lm2index(l-lambda, mu-m)]);
                    }
                    else if (m-mu == 0) {
                        cur_M = Ml0_of_child[l-lambda];
                    }
                    else {
                        cur_M = Mlm_of_child[lm2index(l-lambda, m-mu)];
                    }

                    if (mu < 0) {
                        cur_R = pow(-1, -mu) * conj(Rlm[lm2index(lambda, -mu)]);
                    }
                    else if (mu == 0) {
                        cur_R = Rl0[lambda];
                    }
                    else {
                        cur_R = Rlm[lm2index(lambda, mu)];
                    }

                    Mlm[j] = Mlm[j] + cur_M*conj(cur_R);

                }
            }

            j = j + 1;
        }

    }

}


/*
 * Computes the local expansion due to a multipole moment at a certain distance
 *
 * Note: This is the most expensive part of the fast-multipole expansion
 */
void multipole_to_local(double *dr,             /* shape [3] */
                        int l_max,
                        double *Ml0,            /* shape [0:l_max] */
                        double complex *Mlm,    /* shape [lm2index(l_max, l_max)] */
                        double *Ll0,            /* shape [0:l_max] */
                        double complex *Llm     /* shape [lm2index(l_max, l_max)] */
                        )
{
    double x, costh, phi;

    double Il0[2*l_max+1];
    double complex Ilm[lm2index(2*l_max, 2*l_max)+1];

    int l, lambda, m, mu, k;

    double complex cur_M, cur_I;


    cartesian2spherical(dr, &x, &costh, &phi);

    solid_harmonic_I(x, costh, phi, 2*l_max, Il0, Ilm);

    k = 0;
    /* l_loop */
    for (l = 0; l <= l_max; ++l) {

        /*
         * m = 0
         */

        /* lambda_loop1 */
        for (lambda = 0; lambda <= l_max; ++lambda) {

            /*
             * mu = 0
             */

            Ll0[l] = Ll0[l]
                + Ml0[lambda]*Il0[l+lambda];

            /*
             * mu != 0
             */

            for (mu = 1; mu <= lambda; ++mu) {
                Ll0[l] = Ll0[l]
                    + 2*creal(Mlm[lm2index(lambda, mu)]*Ilm[lm2index(l+lambda, mu)]);
            }

        }
        /* lambda_loop1 */

        /*
         * m != 0
         */

        /* m_loop */
        for (m = 1; m <= l; ++m) {
            /* lambda_loop2 */
            for (lambda = 0; lambda <= l_max; ++lambda) {

                /* mu_loop */
                for (mu = -lambda; mu <= lambda; ++mu) {

                    if (m+mu < 0) {
                        cur_I = pow(-1, -mu-m) * conj(Ilm[lm2index(l+lambda, -m-mu)]);
                    }
                    else if (m+mu == 0) {
                        cur_I = Il0[l+lambda];
                    }
                    else {
                        cur_I = Ilm[lm2index(l+lambda, m+mu)];
                    }

                    if (mu < 0) {
                        cur_M = pow(-1, -mu) * conj(Mlm[lm2index(lambda, -mu)]);
                    }
                    else if (mu == 0) {
                        cur_M = Ml0[lambda];
                    }
                    else {
                        cur_M = Mlm[lm2index(lambda, mu)];
                    }

                    Llm[k] = Llm[k]
                        + cur_M*cur_I;
                }
                /* mu_loop */
            }
            /* lambda_loop2 */

            k = k + 1;
        }
        /* m_loop */

    }
    /* l_loop */
}


/*
 * Shift the local expansion
 */
void local_to_local(double *dr,                /* shape [3] */
                    int l_max_in, 
                    double *Ll0,               /* shape [0:l_max_in] */
                    double complex *Llm,       /* shape [lm2index(l_max_in, l_max_in)] */
                    int l_max_out, 
                    double *Ll0_out,           /* shape [0:l_max_out] */
                    double complex *Llm_out    /* shape [lm2index(l_max_out, l_max_out)] */
                    )
{
    double x, costh, phi;

    double Rl0[l_max_in+1];
    double complex Rlm[lm2index(l_max_in, l_max_in)+1];

    int l, lambda, m, mu, k;

    double complex cur_L, cur_R;


    cartesian2spherical(dr, &x, &costh, &phi);

    solid_harmonic_R(x, costh, phi, l_max_in, Rl0, Rlm);
    for (l = 0; l <= lm2index(l_max_in, l_max_in); ++l) {
        Rlm[l]  = conj(Rlm[l]);
    }

    k = 0;
    /* l_loop2 */
    for (l = 0; l <= l_max_out; ++l) {

        /*
         * m = 0
         */

        /* lambda_loop3 */
        for (lambda = 0; lambda <= l_max_in - l; ++lambda) {

            /*
             * mu = 0
             */

            Ll0_out[l] = Ll0_out[l]
                + Ll0[l+lambda]*Rl0[lambda];

            /*
             * mu != 0
             */

            for (mu = 1; mu <= lambda; ++mu) {
                Ll0_out[l] = Ll0_out[l]
                    + 2*creal(Llm[lm2index(l+lambda, mu)]*Rlm[lm2index(lambda, mu)]);
            }

        }
        /* lambda_loop3 */

        /*
         * m != 0
         */

        /* m_loop2 */
        for (m = 1; m <= l; ++m) {
            /* lambda_loop4 */ 
            for (lambda = 0; lambda <= l_max_in - l; ++lambda) {

                /* mu_loop2 */
                for (mu = -lambda; mu <= lambda; ++mu) {
                    if (m+mu < 0) {
                        cur_L = pow(-1, -mu-m) * conj(Llm[lm2index(l+lambda, -m-mu)]);
                    }
                    else if (m+mu == 0) {
                        cur_L = Ll0[l+lambda];
                    }
                    else {
                        cur_L = Llm[lm2index(l+lambda, m+mu)];
                    }

                    if (mu < 0) {
                        cur_R = pow(-1, -mu) * conj(Rlm[lm2index(lambda, -mu)]);
                    }
                    else if (mu == 0) {
                        cur_R = Rl0[lambda];
                    }
                    else {
                        cur_R = Rlm[lm2index(lambda, mu)];
                    }

                    Llm_out[k] = Llm_out[k] 
                        + cur_L*cur_R;

                }
                /* mu_loop2 */

            }
            /* lambda_loop4 */

            k = k + 1;
        }
        /* m_loop2 */

    }
    /* l_loop2 */

}


/*
 * Python interface
 */


int check_size_and_type(int l_max, PyObject *l0, PyObject *lm)
{
    if (!PyArray_ISFLOAT(l0)) {
        PyErr_SetString(PyExc_TypeError, "l0 needs to be float.");
        return 0;
    }

    if (!PyArray_ISCOMPLEX(lm)) {
        PyErr_SetString(PyExc_TypeError, "lm needs to be complex.");
        return 0;
    }

    if (PyArray_DIM(l0, 0) != l_max+1) {
        PyErr_SetString(PyExc_ValueError, "l0 needs to have size l_max+1.");
        return 0;
    }

    if (PyArray_DIM(lm, 0) != lm2index(l_max, l_max)+1) {
        PyErr_SetString(PyExc_ValueError, "lm needs to have size lm2index(l_max, l_max)+1.");
        return 0;
    }

    return 1;
}


int check_size_and_type_or_allocate(int l_max, PyObject **l0, PyObject **lm)
{
    npy_intp dims[1];

    if (*l0) {
        if (!PyArray_ISFLOAT(*l0)) {
            PyErr_SetString(PyExc_TypeError, "l0 needs to be float.");
            return 0;
        }

        if (PyArray_DIM(*l0, 0) != l_max+1) {
            PyErr_SetString(PyExc_ValueError, "l0 needs to have size l_max+1.");
            return 0;
        }

        Py_INCREF(*l0);
    }
    else {
        dims[0] = l_max+1;
        *l0 = PyArray_ZEROS(1, dims, NPY_DOUBLE, NPY_FALSE);
    }

    if (*lm) {
        if (!PyArray_ISCOMPLEX(*lm)) {
            PyErr_SetString(PyExc_TypeError, "lm needs to be complex.");
            return 0;
        }

        if (PyArray_DIM(*lm, 0) != lm2index(l_max, l_max)+1) {
            PyErr_SetString(PyExc_ValueError, "lm needs to have size lm2index(l_max, l_max)+1.");
            return 0;
        }

        Py_INCREF(*lm);
    }
    else {
        dims[0] = lm2index(l_max, l_max)+1;
        *lm = PyArray_ZEROS(1, dims, NPY_CDOUBLE, NPY_FALSE);
    }

    return 1;
}


PyObject *py_multipole_to_multipole(PyObject *self, PyObject *args)
{
    PyObject *dr, *Ml0_of_child, *Mlm_of_child;
    int l_max;
    PyObject *Ml0 = NULL;
    PyObject *Mlm = NULL;

    if (!PyArg_ParseTuple(args, "O!iO!O!|O!O!",
                          &PyArray_Type, &dr, &l_max,
                          &PyArray_Type, &Ml0_of_child,
                          &PyArray_Type, &Mlm_of_child,
                          &PyArray_Type, &Ml0,
                          &PyArray_Type, &Mlm))
        return NULL;

    if (!PyArray_ISFLOAT(dr)) {
        PyErr_SetString(PyExc_TypeError, "Distance vector needs to be float.");
        return NULL;
    }

    if (PyArray_DIM(dr, 0) != 3) {
        PyErr_SetString(PyExc_TypeError, "Distance vector needs to be 3-dimensional.");
        return NULL;
    }

    if (!check_size_and_type(l_max, Ml0_of_child, Mlm_of_child))
        return NULL;

    if (!check_size_and_type_or_allocate(l_max, &Ml0, &Mlm))
        return NULL;

    multipole_to_multipole(PyArray_DATA(dr), l_max,
                           PyArray_DATA(Ml0_of_child), PyArray_DATA(Mlm_of_child), 
                           PyArray_DATA(Ml0), PyArray_DATA(Mlm));

    return Py_BuildValue("OO", Ml0, Mlm);
}


PyObject *py_multipole_to_local(PyObject *self, PyObject *args)
{
    PyObject *dr, *Ml0, *Mlm;
    int l_max;
    PyObject *Ll0 = NULL;
    PyObject *Llm = NULL;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "O!iO!O!|O!O!",
                          &PyArray_Type, &dr, &l_max,
                          &PyArray_Type, &Ml0,
                          &PyArray_Type, &Mlm,
                          &PyArray_Type, &Ll0,
                          &PyArray_Type, &Llm))
        return NULL;

    if (!PyArray_ISFLOAT(dr)) {
        PyErr_SetString(PyExc_TypeError, "Distance vector needs to be float.");
        return NULL;
    }

    if (PyArray_DIM(dr, 0) != 3) {
        PyErr_SetString(PyExc_TypeError, "Distance vector needs to be 3-dimensional.");
        return NULL;
    }

    if (!check_size_and_type(l_max, Ml0, Mlm))
        return NULL;

    if (!check_size_and_type_or_allocate(l_max, &Ll0, &Llm))
        return NULL;

    dims[0] = l_max+1;
    Ll0 = PyArray_ZEROS(1, dims, NPY_DOUBLE, NPY_FALSE);

    dims[0] = lm2index(l_max, l_max)+1;
    Llm = PyArray_ZEROS(1, dims, NPY_CDOUBLE, NPY_FALSE);

    multipole_to_local(PyArray_DATA(dr), l_max,
                           PyArray_DATA(Ml0), PyArray_DATA(Mlm),
                           PyArray_DATA(Ll0), PyArray_DATA(Llm));

    return Py_BuildValue("OO", Ll0, Llm);
}


PyObject *py_local_to_local(PyObject *self, PyObject *args)
{
    PyObject *dr, *Ll0, *Llm;
    int l_max_in, l_max_out;
    PyObject *Ll0_out = NULL;
    PyObject *Llm_out = NULL;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "O!iO!O!i|O!O!",
                          &PyArray_Type, &dr,
                          &l_max_in,
                          &PyArray_Type, &Ll0,
                          &PyArray_Type, &Llm,
                          &l_max_out,
                          &PyArray_Type, &Ll0_out,
                          &PyArray_Type, &Llm_out))
        return NULL;

    if (!PyArray_ISFLOAT(dr)) {
        PyErr_SetString(PyExc_TypeError, "Distance vector needs to be float.");
        return NULL;
    }

    if (PyArray_DIM(dr, 0) != 3) {
        PyErr_SetString(PyExc_TypeError, "Distance vector needs to be 3-dimensional.");
        return NULL;
    }

    if (!check_size_and_type(l_max_in, Ll0, Llm))
        return NULL;

    if (!check_size_and_type_or_allocate(l_max_out, &Ll0_out, &Llm_out))
        return NULL;

    dims[0] = l_max_out+1;
    Ll0_out = PyArray_ZEROS(1, dims, NPY_DOUBLE, NPY_FALSE);

    dims[0] = lm2index(l_max_out, l_max_out)+1;
    Llm_out = PyArray_ZEROS(1, dims, NPY_CDOUBLE, NPY_FALSE);

    local_to_local(PyArray_DATA(dr),
                   l_max_in, PyArray_DATA(Ll0), PyArray_DATA(Llm),
                   l_max_out, PyArray_DATA(Ll0_out), PyArray_DATA(Llm_out));

    return Py_BuildValue("OO", Ll0_out, Llm_out);
}


