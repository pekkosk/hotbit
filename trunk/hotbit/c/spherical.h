#ifndef __SPHERICAL_H
#define __SPHERICAL_H

#include <Python.h>

#include <complex.h>

#define lm2index(l, m)  ((l)*(l-1)/2 + (m) - 1)


void cartesian_to_spherical(double *r, double *x, double *costh, double *phi);
void solid_harmonic_R(double x, double costh, double phi, int l_max, double *Rl0, double complex *Rlm);
void solid_harmonic_I(double x, double costh, double phi, int l_max, double *Il0, double complex *Ilm);

PyObject *py_cartesian_to_spherical(PyObject *self, PyObject *args);
PyObject *py_solid_harmonic_R(PyObject *self, PyObject *args);
PyObject *py_solid_harmonic_I(PyObject *self, PyObject *args);

#endif
