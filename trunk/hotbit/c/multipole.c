/*
 * Copyright (C) 2010 NSC Jyvaskyla, Fh-IWM
 * Please see the accompanying LICENSE file for further information.
 */

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL HOTBIT_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include <complex.h>

#include "spherical.h"

#include "multipole.h"

#define min(a,b)        ((a)>(b)?(b):(a))
#define max(a,b)        ((a)<(b)?(b):(a))


/*
 * Helper module for the fast-multipole solver.
 *
 * Multipole-to-multipole, multipole-to-local
 * and local-to-local transformations.
 *
 * Multipole transformations (important for rotations).
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


    if (dr[0] == 0.0 && dr[1] == 0.0 && dr[2] == 0.0) {

        for (l = 0; l <= l_max; ++l) {
            Ml0[l] += Ml0_of_child[l];
        }

        for (l = 0; l <= lm2index(l_max, l_max); ++l) {
            Mlm[l] += Mlm_of_child[l];
        }

    } else {

		cartesian_to_spherical(dr, &x, &costh, &phi);

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
					for (mu = max(-lambda, lambda-l+m);
						 mu <= min(lambda, -lambda+l+m); ++mu) {

						if (m-mu < 0) {
							cur_M = pow(-1, mu-m)
								* conj(Mlm_of_child[lm2index(l-lambda, mu-m)]);
						}
						else if (m-mu == 0) {
							cur_M = Ml0_of_child[l-lambda];
						}
						else {
							cur_M = Mlm_of_child[lm2index(l-lambda, m-mu)];
						}

						if (mu < 0) {
							cur_R = pow(-1, -mu) 
								* conj(Rlm[lm2index(lambda, -mu)]);
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

    } /* dr == 0.0 */

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


    cartesian_to_spherical(dr, &x, &costh, &phi);

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


    if (dr[0] == 0.0 && dr[1] == 0.0 && dr[2] == 0.0) {

        m = min(l_max_in, l_max_out);

        for (l = 0; l <= m; ++l) {
            Ll0_out[l] += Ll0[l];
        }

        for (l = 0; l <= lm2index(m, m); ++l) {
            Llm_out[l] += Llm[l];
        }

    }
    else {

        /* Compute R for -dr */
        cartesian_to_spherical(dr, &x, &costh, &phi);

        solid_harmonic_R(x, -costh, M_PI+phi, l_max_in, Rl0, Rlm);

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
                        + 2*creal(Llm[lm2index(l+lambda, mu)]*
                                  conj(Rlm[lm2index(lambda, mu)]));
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
                            cur_L = pow(-1, -mu-m) *
                                conj(Llm[lm2index(l+lambda, -m-mu)]);
                        }
                        else if (m+mu == 0) {
                            cur_L = Ll0[l+lambda];
                        }
                        else {
                            cur_L = Llm[lm2index(l+lambda, m+mu)];
                        }
                        
                        if (mu < 0) {
                            cur_R = pow(-1, -mu) * Rlm[lm2index(lambda, -mu)];
                        }
                        else if (mu == 0) {
                            cur_R = Rl0[lambda];
                        }
                        else {
                            cur_R = conj(Rlm[lm2index(lambda, mu)]);
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

}


/*
 * All matrices have zero padding of width 2 since
 * in the sum i,j can actually become > l and < -l.
 */
#define idx(l,i,j)  ((2+l+j) + (2+l+i)*(2*(l+2)+1))
#define arrsize(l)  (idx(l, l+2, l+2)+1)

/*
 * Multiply a multipole by a transformation matrix D for
 * angular momentum l
 */
void matrix_dot_multipole(int l,
                          double complex *D,
                          double *Rl0,
                          double complex *Rlm,
                          double *Sl0,
                          double complex *Slm
                          )
{
    int m, n;

#ifdef DEBUG
    double complex c;

    printf("=== l = %i ===\n", l);

    for (m = -l; m <= l; ++m) {

        for (n = -l; n <= l; ++n) {
            printf("%7.3f+i%7.3f",
                   creal(D[idx(l, m, n)]), cimag(D[idx(l, m, n)]));

            c = pow(-1, m+n)*conj(D[idx(l, m, n)]);
            
            if (cabs(D[idx(l, -m, -n)] - c) > 1e-3) {
                printf("*   ");
            }
            else {
                printf("    ");
            }
        }

        printf("\n");
    }
#endif

    /* diagonal component */
    Sl0[l] = creal( D[idx(l, 0, 0)] ) * Rl0[l];
    for (m = 1; m <= l; ++m) {
        Sl0[l] +=
            creal( D[idx(l, 0,  m)]              *      Rlm[lm2index(l, m)]
                 + D[idx(l, 0, -m)] * pow(-1, m) * conj(Rlm[lm2index(l, m)]) );
    }

    /* off-diagonal components */
    for (n = 1; n <= l; ++n) {
        /* n = 0 contribution */
        Slm[lm2index(l, n)] = D[idx(l, n, 0)] * Rl0[l];

        for (m = 1; m <= l; ++m) {
            Slm[lm2index(l, n)] += 
                ( D[idx(l, n,  m)]              *      Rlm[lm2index(l, m)]
                + D[idx(l, n, -m)] * pow(-1, m) * conj(Rlm[lm2index(l, m)]) );
        }
    }
}


/*
 * Transform a multipole given the cartesian transformation matrix T
 *
 * The algorithm is derived and described in
 *    C.H. Choi, J. Ivanic, M.S. Gordon, K. Ruedenberg,
 *    J. Chem. Phys. 111, 8825 (1999)
 *
 * FIXME! This implementation currently computes D(l,m,n)
 * and D(l,-m,-n) which is unnessecary, since
 * D(l,-m,-n) = (-1)*(m+n) D(l,m,n). Bookkeeping becomes tedious
 * if one makes use of that symmetry though.
 */
void transform_multipole(double *T,                /* shape [3,3] */
                         int l_max,
                         double *Rl0,              /* shape [0:l_max] */
                         double complex *Rlm,      /* shape [lm2index(l_max, l_max)] */
                         double *Sl0,              /* shape [0:l_max] */
                         double complex *Slm       /* shape [lm2index(l_max, l_max)] */
                         )
{
    double complex D1[arrsize(1)];
    double complex D[arrsize(l_max)];
    double complex Dl[arrsize(l_max)];

    int m, n, l;


    /* FIXME!!! Tabulate integer sqrts */
    /* FIXME!!! Also tabulate (-1)**m? */
    /*
#define a(l,m,n)  sqrt( (l+m)*(l-m)/( (l+n)*(l-n) ) )
#define b(l,m,n)  sqrt( (l+m)*(l+m-1)/( 2*(l+n)*(l-n) ) )
#define c(l,m,n)  sqrt( 2*(l+m)*(l-m)/( (l+n)*(l+n-1) ) )
#define d(l,m,n)  sqrt( (l+m)*(l+m-1)/( (l+n)*(l+n-1) ) )
    */

#define a(l,m)  ( ( (double) ( (l)-(m)+1 ) * ( (l)-(m) ) ) / 2 )
#define b(l,m)  ( ( (double) ( (l)+(m)   ) * ( (l)-(m) ) )     )
#define c(l,m)  ( ( (double) ( (l)+(m)+1 ) * ( (l)+(m) ) ) / 2 )

    /*
#define aa(l,m,n)  ( ( (double) ( (l)+(n) ) * ( (l)+(n)-1 ) ) / ( 2*( ( (l)+(m) ) * ( (l)-(m) ) ) ) )
#define bb(l,m,n)  ( ( (double) ( (l)+(n) ) * ( (l)+(n)   ) ) / (   ( ( (l)+(m) ) * ( (l)-(m) ) ) ) )
#define cc(l,m,n)  ( ( (double) ( (l)+(n) ) * ( (l)+(n)+1 ) ) / ( 2*( ( (l)+(m) ) * ( (l)-(m) ) ) ) )
    */

    /* l = 0 does not transform under rotation */
    Sl0[0]  = Rl0[0];

#define idx2(i,j)  (j+3*i)
#define _XX_ idx2(0,0)
#define _YY_ idx2(1,1)
#define _ZZ_ idx2(2,2)
#define _XY_ idx2(0,1)
#define _YZ_ idx2(1,2)
#define _ZX_ idx2(2,0)
#define _YX_ idx2(1,0)
#define _ZY_ idx2(2,1)
#define _XZ_ idx2(0,2)

    bzero(D1, arrsize(1)*sizeof(double complex));
    bzero(D, arrsize(l_max)*sizeof(double complex));
    bzero(Dl, arrsize(l_max)*sizeof(double complex));

    /*
     * l = 1 starts the recurrence relation.
     * The conjugation is required because the multipole moments
     * are sums of R* rather than R.
     */
    D1[idx(1, -1, -1)]  = conj( T[_YY_] + T[_XX_]  + I*(T[_XY_] - T[_YX_]) )/2;
    D1[idx(1, -1,  0)]  = conj( T[_XZ_]            - I*T[_YZ_]             )/2;
    D1[idx(1, -1,  1)]  = conj( T[_YY_] - T[_XX_]  + I*(T[_XY_] + T[_YX_]) )/2;

    D1[idx(1,  0, -1)]  = conj( T[_ZX_]            + I*T[_ZY_]             );
    D1[idx(1,  0,  0)]  = conj( T[_ZZ_]                                    );
    D1[idx(1,  0,  1)]  = conj( -T[_ZX_]           + I*T[_ZY_]             );

    D1[idx(1,  1, -1)]  = conj( -T[_XX_] + T[_YY_] - I*(T[_XY_] + T[_YX_]) )/2;
    D1[idx(1,  1,  0)]  = conj( -T[_XZ_]           - I*T[_YZ_]             )/2;
    D1[idx(1,  1,  1)]  = conj( T[_XX_] + T[_YY_]  - I*(T[_XY_] - T[_YX_]) )/2;

    memcpy(D, D1, arrsize(1)*sizeof(double complex));

    matrix_dot_multipole(1, D, Rl0, Rlm, Sl0, Slm);

    /* l = 2...l_max */
    for (l = 2; l <= l_max; ++l) {
        /* copy last matrix to Dl */
        memcpy(Dl, D, arrsize(l-1)*sizeof(double complex));
        bzero(D, arrsize(l)*sizeof(double complex));

        /* transformation rules */
        for (n = -l; n <= l; ++n) {

            /* m = -l */
            D[idx(l, -l, n)] =
                (   a(l, n+1) * D1[idx(1, -1, -1)] * Dl[idx(l-1, -l+1, n+1)]
                  + b(l, n  ) * D1[idx(1, -1,  0)] * Dl[idx(l-1, -l+1, n  )]
                  + c(l, n-1) * D1[idx(1, -1,  1)] * Dl[idx(l-1, -l+1, n-1)] )
                /a(l, -l+1);

            /* -l+1 <= m <= l-1 */
            for (m = -l+1; m <= l-1; ++m) {

                D[idx(l, m, n)] =
                    (   a(l, n+1) * D1[idx(1,  0, -1)] * Dl[idx(l-1, m, n+1)]
                      + b(l, n  ) * D1[idx(1,  0,  0)] * Dl[idx(l-1, m, n  )]
                      + c(l, n-1) * D1[idx(1,  0,  1)] * Dl[idx(l-1, m, n-1)] )
                      /b(l, m);

            }

            /* m = l */
            D[idx(l,  l, n)] =
                (   a(l, n+1) * D1[idx(1,  1, -1)] * Dl[idx(l-1, l-1, n+1)]
                  + b(l, n  ) * D1[idx(1,  1,  0)] * Dl[idx(l-1, l-1, n  )]
                  + c(l, n-1) * D1[idx(1,  1,  1)] * Dl[idx(l-1, l-1, n-1)] )
                /c(l, l-1);

        }

        matrix_dot_multipole(l, D, Rl0, Rlm, Sl0, Slm);
    }

#undef a
#undef b
#undef c

}

#undef idx



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
    PyObject *ret;

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

    /* New tuple, SET_ITEM does not increase reference count */
    ret = PyTuple_New(2);
    PyTuple_SET_ITEM(ret, 0, Ml0);
    PyTuple_SET_ITEM(ret, 1, Mlm);

    return ret;
}


PyObject *py_multipole_to_local(PyObject *self, PyObject *args)
{
    PyObject *dr, *Ml0, *Mlm;
    int l_max;
    PyObject *Ll0 = NULL;
    PyObject *Llm = NULL;
    PyObject *ret;

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

    multipole_to_local(PyArray_DATA(dr), l_max,
		       PyArray_DATA(Ml0), PyArray_DATA(Mlm),
		       PyArray_DATA(Ll0), PyArray_DATA(Llm));

    /* New tuple, SET_ITEM does not increase reference count */
    ret = PyTuple_New(2);
    PyTuple_SET_ITEM(ret, 0, Ll0);
    PyTuple_SET_ITEM(ret, 1, Llm);

    return ret;
}


PyObject *py_local_to_local(PyObject *self, PyObject *args)
{
    PyObject *dr, *Ll0, *Llm;
    int l_max_in, l_max_out;
    PyObject *Ll0_out = NULL;
    PyObject *Llm_out = NULL;
    PyObject *ret;

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

    local_to_local(PyArray_DATA(dr),
                   l_max_in, PyArray_DATA(Ll0), PyArray_DATA(Llm),
                   l_max_out, PyArray_DATA(Ll0_out), PyArray_DATA(Llm_out));

    /* New tuple, SET_ITEM does not increase reference count */
    ret = PyTuple_New(2);
    PyTuple_SET_ITEM(ret, 0, Ll0_out);
    PyTuple_SET_ITEM(ret, 1, Llm_out);

    return ret;
}


PyObject *py_transform_multipole(PyObject *self, PyObject *args)
{
    PyObject *R, *Rl0, *Rlm;
    int l_max;
    PyObject *Sl0;
    PyObject *Slm;
    npy_intp dims[1];
    PyObject *ret;

    if (!PyArg_ParseTuple(args, "O!iO!O!",
                          &PyArray_Type, &R,
                          &l_max,
                          &PyArray_Type, &Rl0,
                          &PyArray_Type, &Rlm))
        return NULL;

    if (!PyArray_ISFLOAT(R)) {
        PyErr_SetString(PyExc_TypeError, "Rotation matrix needs to be float.");
        return NULL;
    }

    if (!PyArray_NDIM(R) == 2) {
        PyErr_SetString(PyExc_TypeError, "Rotation matrix needs to be 3x3.");
        return NULL;
    }

    if (PyArray_DIM(R, 0) != 3 || PyArray_DIM(R, 1) != 3) {
        PyErr_SetString(PyExc_TypeError, "Rotation matrix needs to be 3x3.");
        return NULL;
    }

    if (!check_size_and_type(l_max, Rl0, Rlm))
        return NULL;

    dims[0] = l_max+1;
    //    Sl0 = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    Sl0 = PyArray_ZEROS(1, dims, NPY_DOUBLE, NPY_FALSE);

    dims[0] = lm2index(l_max, l_max)+1;
    //    Slm = PyArray_SimpleNew(1, dims, NPY_CDOUBLE);
    Slm = PyArray_ZEROS(1, dims, NPY_CDOUBLE, NPY_FALSE);

    transform_multipole(PyArray_DATA(R), l_max,
                     PyArray_DATA(Rl0), PyArray_DATA(Rlm),
                     PyArray_DATA(Sl0), PyArray_DATA(Slm));


    /* New tuple, SET_ITEM does not increase reference count */
    ret = PyTuple_New(2);
    PyTuple_SET_ITEM(ret, 0, Sl0);
    PyTuple_SET_ITEM(ret, 1, Slm);

    return ret;
}



