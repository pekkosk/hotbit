#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL HOTBIT_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include <complex.h>

#include "spherical.h"


/*
 * Sperical symmetry module
 *
 * All kinds of stuff related to spherical symmetry, i.e.
 * - Legendre polynomials
 * - solid harmonics
 */



/*
 * Transform cartesian into spherical coordinates
 */
void cartesian_to_spherical(double *r,         /* shape [3], distance vector */
                         double *x,         /* absolute distance */
                         double *costh,     /* cos(theta) angle */
                         double *phi        /* phi angle */
                         )
{
    *x      = sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
    *costh  = r[2]/(*x);
    *phi    = atan2( r[1], r[0] );
}


/*
 * Compute the solid harmonics up to l=l_max.
 * Note that the normalization is non-standard.
 *
 *     
 *   m      1         l    i m phi    m
 * R   = -------    r    e          P   (cos theta)
 *   l    (l+m)!                      l
 */
void solid_harmonic_R(double x, 
                      double costh, 
                      double phi, 
                      int l_max, 
                      double *Rl0,            /* shape [0:l_max] */
                      double complex *Rlm     /* shape [lm2index(l_max, l_max)] */
                      )
{
    int l, m, j;
    double Pl0[l_max+1], Plm[lm2index(l_max, l_max)+1], xl, sinth, h;
    double complex ang[l_max];

    sinth = sqrt(1-costh*costh);

    /*
     * Compute the associated legendre polynomials divided by (l+m)!
     */

    Pl0[0]               = 1.0;
    Pl0[1]               = costh;
    Plm[lm2index(1, 1)]  = -sinth/2;

    for (l = 2; l <= l_max; ++l) {
        h  = sinth/(2*l);
        Plm[lm2index(l,   l)]  = -h*Plm[lm2index(l-1, l-1)];

        Plm[lm2index(l, l-1)]  =  costh*Plm[lm2index(l-1, l-1)];

        Pl0[l] = ( (2*l-1)*costh*Pl0[l-1] - Pl0[l-2] )/(l*l);

        j = lm2index(l, 0);
        for (m = 1; m <= l-2; ++m) {
            j = j+1;

            Plm[j] = ( (2*l-1)*costh*Plm[lm2index(l-1, m)] -
                       Plm[lm2index(l-2, m)] )/((l+m)*(l-m));
        }
    }

    /*
     * Compute angular (phi) terms
     */

    for (m = 1; m <= l_max; ++m) {
        ang[m-1] = cexp(I*m*phi);
    }

    /*
     * Compute solid harmonics
     */

    Rl0[0]               = 1.0;
    Rl0[1]               = x*Pl0[1];
    Rlm[lm2index(1, 1)]  = x*Plm[lm2index(1, 1)]*ang[0];

    xl    = x; /* x^l */
    for (l = 2; l <= l_max; ++l) {
        xl    = xl*x;

        Rl0[l] = xl*Pl0[l];

        j = lm2index(l, 0);
        for (m = 1; m <= l; ++m) {
            j       = j+1;

            Rlm[j]  = xl*Plm[j]*ang[m-1];
        }
    }
}


/*
 * Compute the solid harmonics up to l=l_max
 * Note that the normalization is non-standard.
 *
 *   m              -l+1    i m phi    m
 * I   = (l-m)!   r       e          P   (cos theta)
 *   l                                 l
 */
void solid_harmonic_I(double x,
                      double costh,
                      double phi,
                      int l_max,
                      double *Il0,            /* shape [0:l_max] */
                      double complex *Ilm     /* shape [lm2index(l_max, l_max)] */
                      )
{
    int l, m, j, h;

    double sinth, h_sinth, h_costh;
    double Pl0[l_max+1];
    double Plm[lm2index(l_max, l_max)+1];
    double xl, xlx, rec_x;
    double complex ang[l_max];

    sinth = sqrt(1-costh*costh);

    rec_x = 1.0/x;

    /*
     * Compute the associated legendre polynomials times by (l-m)!
     * and their derivatives
     *
     * Note: Plm is the associated Legendre polynomial
     * DIVIDED BY SIN(THETA). This way we don't run into
     * divide by zeros in the phi derivatives.
     */

    Pl0[0]               = 1.0;
    Pl0[1]               = costh;
    Plm[lm2index(1, 1)]  = -1.0;

    for (l = 2; l <= l_max; ++l) {
        h        = 2*l-1;
        h_sinth  = h*sinth;
        h_costh  = h*costh;
        Plm[lm2index(l,    l)]  = -h_sinth*Plm[lm2index(l-1,  l-1)];

        Plm[lm2index(l,  l-1)]  =  h_costh*Plm[lm2index(l-1,  l-1)];

        Pl0[l]                  = h_costh*Pl0[l-1] - (l-1)*(l-1)*Pl0[l-2];

        for (m = 1; m <= l-2; ++m) {
            Plm[lm2index(l, m)] = h_costh*Plm[lm2index(l-1, m)] - (l-m-1)*(l+m-1)*Plm[lm2index(l-2, m)];
        }
    }   

    /*
     * Compute angular (phi) terms
     */

    for (m = 1; m <= l_max; ++m) {
        ang[m-1] = cexp(I*m*phi);
    }

    /*
     * Compute solid harmonics
     */

    Il0[0]               = 1.0*rec_x;
    xl                   = rec_x*rec_x;
    Il0[1]               = Pl0[1]*xl;
    Ilm[lm2index(1, 1)]  = sinth*Plm[lm2index(1, 1)]*ang[0]*xl;

    xlx                  = xl*rec_x;

    for (l = 2; l <= l_max; ++l) {
        xl    = xlx;
        xlx   = xl*rec_x;

        Il0[l]     = Pl0[l]*xl;

        j          = lm2index(l, 0);
        for (m = 1; m <= l; ++m) {
            j = j+1;

            Ilm[j] = sinth*Plm[j]*ang[m-1]*xl;
        }
    }
}


/*
 * Python interface
 */


PyObject *py_cartesian_to_spherical(PyObject *self, PyObject *args)
{
    PyObject *r, *s;
    double *sd;
    npy_intp dims[1];

    if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &r))
        return NULL;

    if (!PyArray_ISFLOAT(r)) {
        PyErr_SetString(PyExc_TypeError, "Distance vector needs to be float.");
        return NULL;
    }

    if (PyArray_DIM(r, 0) != 3) {
        PyErr_SetString(PyExc_TypeError, "Distance vector needs to be 3-dimensional.");
        return NULL;
    }

    dims[0] = 3;
    s = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    sd = PyArray_DATA(s);

    cartesian_to_spherical(PyArray_DATA(r), &sd[0], &sd[1], &sd[2]);

    Py_DECREF(r);

    return s;
}


PyObject *py_solid_harmonic_R(PyObject *self, PyObject *args)
{
    PyObject *dr;
    double *drd;
    int l_max;

    double x, costh, phi;

    npy_intp dims[1];
    PyObject *Rl0 = NULL;
    PyObject *Rlm = NULL;


    if (!PyArg_ParseTuple(args, "O!i|O!O!",
                          &PyArray_Type, &dr,
                          &l_max,
                          &PyArray_Type, &Rl0,
                          &PyArray_Type, &Rlm
                          ))
        return NULL;

    if (!PyArray_ISFLOAT(dr)) {
        PyErr_SetString(PyExc_TypeError, "Distance vector needs to be float.");
        return NULL;
    }

    if (Rl0) {
        if (!PyArray_ISFLOAT(Rl0)) {
            PyErr_SetString(PyExc_TypeError, "Rl0 needs to be float.");
            return NULL;
        }

        if (PyArray_DIM(Rl0, 0) != l_max+1) {
            PyErr_SetString(PyExc_ValueError, "Rl0 needs to have size l_max+1.");
            return NULL;
        }
        
        Py_INCREF(Rl0);
    }
    else {
        dims[0] = l_max+1;
        Rl0 = PyArray_ZEROS(1, dims, NPY_DOUBLE, NPY_FALSE);
    }

    if (Rlm) {
        if (!PyArray_ISCOMPLEX(Rlm)) {
            PyErr_SetString(PyExc_TypeError, "Rlm needs to be complex.");
            return NULL;
        }

        if (PyArray_DIM(Rlm, 0) != lm2index(l_max, l_max)+1) {
            PyErr_SetString(PyExc_ValueError, "Rlm needs to have size lm2index(l_max, l_max)+1.");
            return NULL;
        }

        Py_INCREF(Rlm);
    }
    else {
        dims[0] = lm2index(l_max, l_max)+1;
        Rlm = PyArray_ZEROS(1, dims, NPY_CDOUBLE, NPY_FALSE);
    }

    drd = PyArray_DATA(dr);
    if (drd[0] == 0 && drd[1] == 0 && drd[2] == 0) {
        (*(double*) PyArray_DATA(Rl0)) += 1.0;
    } else {
        cartesian_to_spherical(drd, &x, &costh, &phi);
        solid_harmonic_R(x, costh, phi, l_max,
                         PyArray_DATA(Rl0), PyArray_DATA(Rlm));
    }

    return PyTuple_Pack(2, Rl0, Rlm);
}


PyObject *py_solid_harmonic_I(PyObject *self, PyObject *args)
{
    PyObject *dr;
    int l_max;

    double x, costh, phi;

    npy_intp dims[1];
    PyObject *Il0 = NULL;
    PyObject *Ilm = NULL;


    if (!PyArg_ParseTuple(args, "O!i",
                          &PyArray_Type, &dr,
                          &l_max,
                          &PyArray_Type, &Il0,
                          &PyArray_Type, &Ilm))
        return NULL;

    if (!PyArray_ISFLOAT(dr)) {
        PyErr_SetString(PyExc_TypeError, "Distance vector needs to be float.");
        return NULL;
    }

    cartesian_to_spherical(PyArray_DATA(dr), &x, &costh, &phi);

    if (Il0) {
        if (!PyArray_ISFLOAT(Il0)) {
            PyErr_SetString(PyExc_TypeError, "Il0 needs to be float.");
            return NULL;
        }

        if (PyArray_DIM(Il0, 0) != l_max+1) {
            PyErr_SetString(PyExc_ValueError, "Il0 needs to have size l_max+1.");
            return NULL;
        }
        
        Py_INCREF(Il0);
    }
    else {
        dims[0] = l_max+1;
        Il0 = PyArray_ZEROS(1, dims, NPY_DOUBLE, NPY_FALSE);
    }

    if (Ilm) {
        if (!PyArray_ISCOMPLEX(Ilm)) {
            PyErr_SetString(PyExc_TypeError, "Ilm needs to be complex.");
            return NULL;
        }

        if (PyArray_DIM(Ilm, 0) != lm2index(l_max, l_max)+1) {
            PyErr_SetString(PyExc_ValueError, "Ilm needs to have size lm2index(l_max, l_max)+1.");
            return NULL;
        }

        Py_INCREF(Ilm);
    }
    else {
        dims[0] = lm2index(l_max, l_max)+1;
        Ilm = PyArray_ZEROS(1, dims, NPY_CDOUBLE, NPY_FALSE);
    }

    solid_harmonic_I(x, costh, phi, l_max, PyArray_DATA(Il0), PyArray_DATA(Ilm));

    return PyTuple_Pack(2, Il0, Ilm);
}
