/*
 * Copyright (C) 2010 NSC Jyvaskyla, Fh-IWM
 * Please see the accompanying LICENSE file for further information.
 */

#ifndef __MULTIPOLE_H
#define __MUTLIPOLE_H

#include <complex.h>


void multipole_to_multipole(double *dr, int l_max,
                            double *Ml0_of_child, double complex *Mlm_of_child, 
                            double *Ml0, double complex *Mlm);
void multipole_to_local(double *dr, int l_max,
                        double *Ml0, double complex *Mlm,
                        double *Ll0, double complex *Llm);
void local_to_local(double *dr,
                    int l_max_in, 
                    double *Ll0, double complex *Llm,
                    int l_max_out, 
                    double *Ll0_out, double complex *Llm_out);
void transform_multipole(double *R, int l_max,
                      double *Rl0, double complex *Rlm,
                      double *Sl0, double complex *Slm);


PyObject *py_multipole_to_multipole(PyObject *self, PyObject *args);
PyObject *py_multipole_to_local(PyObject *self, PyObject *args);
PyObject *py_local_to_local(PyObject *self, PyObject *args);
PyObject *py_transform_multipole(PyObject *self, PyObject *args);

#endif
