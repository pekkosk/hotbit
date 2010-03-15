#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL HOTBIT_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include "slako.h"

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#define pow2(x)  ((x)*(x))


/*
 * Make Slater-Koster transformations
 *
 * Given direction cosine (rhat from atom i to atom j), distance (dist) and
 * interpolated matrix elements and their derivatives fos s and h,
 * apply transformation rules and return the transformed sub-matrices.
 *
 * noi and noj are the number of orbitals on atom i and j and
 * is 1,4 or 9 corresponding to s, sp and spd-valent elements.
 *
 * Orbital ordering: 
 * s, px, py, pz, dxy, dyz, dzx, dx2-y2, d3z2-r2    
 */
void 
fast_slako_transformations(
  double *rhat,    /* normal direction connecting the two atoms */
  double dist,     /* distance of the atoms */
  int noi,         /* number of orbitals of the first atom */
  int noj,         /* number of orbitals of the second atom */
  double *h,       /* shape [14] bond hamiltonian matrix elements, see orbital ordering above */
  double *s,       /* shape [14] bond overlap matrix elements, see orbital ordering above */
  double *dh,      /* shape [14] derivatives of the bond hamiltonian matrix elements */
  double *ds,      /* shape [14] derivatives of hte bond overlap matrix elements */
  double *ht,      /* shape [noi,noj] transformed hamiltonian matrix (output) */
  double *st,      /* shape [noi,noj] transformed overlap matrix (outpt) */
  double *dht,     /* shape [noi,noj,3] derivative of the transformed hamiltonian matrix (output) */
  double *dst      /* shape [noi,noj,3] derivative of the transformed overlap matrix (output) */
  )
{
    const double s3 = sqrt(3.0);

    int i, j, a, b, maxorb, nc;
    int cnt[9][9];
    int ind[9][9][3];

    double l, m, n, ll, mm, nn;
    double dl[3] = { 1.0, 0.0, 0.0 };
    double dm[3] = { 0.0, 1.0, 0.0 };
    double dn[3] = { 0.0, 0.0, 1.0 };
    double dll[3], dmm[3], dnn[3];
    double mat[9][9][3];
    double der[9][9][3][3];

    /* --- */

    l = rhat[0];
    m = rhat[1];
    n = rhat[2];
    
    ll = l*l;
    mm = m*m;
    nn = n*n;

    for (i = 0; i < 3; ++i) {
        dl[i] = (dl[i]-l*rhat[i])/dist;
        dm[i] = (dm[i]-m*rhat[i])/dist;
        dn[i] = (dn[i]-n*rhat[i])/dist;

        dll[i] = 2*l*dl[i];
        dmm[i] = 2*m*dm[i];
        dnn[i] = 2*n*dn[i];
    }

    for (i = 0; i < 9; ++i) {
        for (j = 0; j < 9; ++j) {

            if (i >= 4 && j >= 4) {
                cnt[i][j] = 3;
            } else if (i >= 1 && j >= 1) {
                cnt[i][j] = 2;
            } else {
                cnt[i][j] = 1;
            }

        }
    }

    bzero(ht, noi*noj*sizeof(double));
    bzero(st, noi*noj*sizeof(double));
    bzero(dht, noi*noj*3*sizeof(double));
    bzero(dst, noi*noj*3*sizeof(double));

    maxorb = max(noi, noj);

    mat[0][0][0] = 1.0;
    bzero(&der[0][0][0], 3*sizeof(double));
    ind[0][0][0] = 9;
    
    if (maxorb >= 2) {
        mat[0][1][0] = l;
        for (i = 0; i < 3; ++i) {
            der[0][1][0][i] = dl[i];
        }
        ind[0][1][0] = 8;
            
        mat[0][2][0] = m;
        for (i = 0; i < 3; ++i) {
            der[0][2][0][i] = dm[i];
        }
        ind[0][2][0] = 8;
                
        mat[0][3][0] = n;
        for (i = 0; i < 3; ++i) {
            der[0][3][0][i] = dn[i];
        }
        ind[0][3][0] = 8;
 
        mat[1][1][0] = ll;
        mat[1][1][1] = 1-ll;
        for (i = 0; i < 3; ++i) {
            der[1][1][0][i] = dll[i];
            der[1][1][1][i] = -dll[i];
        }
        ind[1][1][0] = 5;
        ind[1][1][1] = 6;
            
        mat[1][2][0] = l*m;
        mat[1][2][1] = -l*m;
        for (i = 0; i < 3; ++i) {
            der[1][2][0][i] = dl[i]*m+l*dm[i];
            der[1][2][1][i] = -(dl[i]*m+l*dm[i]);
        }
        ind[1][2][0] = 5.0;
        ind[1][2][1] = 6.0;
            
        mat[1][3][0] = l*n;
        mat[1][3][1] = -l*n;
        for (i = 0; i < 3; ++i) {
            der[1][3][0][i] = dl[i]*n+l*dn[i];
            der[1][3][1][i] = -(dl[i]*n+l*dn[i]);
        }
        ind[1][3][0] = 5;
        ind[1][3][1] = 6;
         
        mat[2][2][0] = mm;
        mat[2][2][1] = 1-mm;
        for (i = 0; i < 3; ++i) {
            der[2][2][0][i] = dmm[i];
            der[2][2][1][i] = -dmm[i];
        }
        ind[2][2][0] = 5;
        ind[2][2][1] = 6;
                
        mat[2][3][0] = m*n;
        mat[2][3][1] = -m*n;
        for (i = 0; i < 3; ++i) {
            der[2][3][0][i] = dm[i]*n+m*dn[i];
            der[2][3][1][i] = -(dm[i]*n+m*dn[i]);
        }
        ind[2][3][0] = 5;
        ind[2][3][1] = 6;

        mat[3][3][0] = nn;
        mat[3][3][1] = 1-nn;
        for (i = 0; i < 3; ++i) {
            der[3][3][0][i] = dnn[i];
            der[3][3][1][i] = -dnn[i];
        }
        ind[3][3][0] = 5;
        ind[3][3][1] = 6;
    }
        
    if (maxorb >=5) {
        mat[0][4][0] = s3*l*m;
        for (i = 0; i < 3; ++i) {
            der[0][4][0][i] = s3*(dl[i]*m+l*dm[i]);
        }
        ind[0][4][0] = 7;
                
        mat[0][5][0] = s3*m*n;
        for (i = 0; i < 3; ++i) {
            der[0][5][0][i] = s3*(dm[i]*n+m*dn[i]);
        }
        ind[0][5][0] = 7;
                
        mat[0][6][0] = s3*n*l;
        for (i = 0; i < 3; ++i) {
            der[0][6][0][i] = s3*(dn[i]*l+n*dl[i]);
        }
        ind[0][6][0] = 7;
                
        mat[0][7][0] = 0.5*s3*(ll-mm);
        for (i = 0; i < 3; ++i) {
            der[0][7][0][i] = 0.5*s3*(dll[i]-dmm[i]);
        }
        ind[0][7][0] = 7;
                
        mat[0][8][0] = nn-0.5*(ll+mm);
        for (i = 0; i < 3; ++i) {
            der[0][8][0][i] = dnn[i]-0.5*(dll[i]+dmm[i]);
        }
        ind[0][8][0] = 7;

        mat[1][4][0] = s3*ll*m;
        mat[1][4][1] = m*(1-2*ll);
        for (i = 0; i < 3; ++i) {
            der[1][4][0][i] = s3*(dll[i]*m+ll*dm[i]);
            der[1][4][1][i] = dm[i]*(1-2*ll)+m*(-2*dll[i]);
        }
        ind[1][4][0] = 3;
        ind[1][4][1] = 4;
                
        mat[1][5][0] = s3*l*m*n;
        mat[1][5][1] = -2*l*m*n;
        for (i = 0; i < 3; ++i) {
            der[1][5][0][i] = s3*(dl[i]*m*n+l*dm[i]*n+l*m*dn[i]);
            der[1][5][1][i] = -2*(dl[i]*m*n+l*dm[i]*n+l*m*dn[i]);
        }
        ind[1][5][0] = 3;
        ind[1][5][1] = 4;
                
        mat[1][6][0] = s3*ll*n;
        mat[1][6][1] = n*(1-2*ll);
        for (i = 0; i < 3; ++i) {
            der[1][6][0][i] = s3*(dll[i]*n+ll*dn[i]);
            der[1][6][1][i] = dn[i]*(1-2*ll)+n*(-2*dll[i]);
        }
        ind[1][6][0] = 3;
        ind[1][6][1] = 4;
                
        mat[1][7][0] = 0.5*s3*l*(ll-mm);
        mat[1][7][1] = l*(1-ll+mm);
        for (i = 0; i < 3; ++i) {
            der[1][7][0][i] = 0.5*s3*(dl[i]*(ll-mm)+l*(dll[i]-dmm[i]));
            der[1][7][1][i] = dl[i]*(1-ll+mm)+l*(-dll[i]+dmm[i]);
        }
        ind[1][7][0] = 3;
        ind[1][7][1] = 4;
                
        mat[1][8][0] = l*(nn-0.5*(ll+mm));
        mat[1][8][1] = -s3*l*nn;
        for (i = 0; i < 3; ++i) {
            der[1][8][0][i] = dl[i]*(nn-0.5*(ll+mm))+l*(dnn[i]-0.5*(dll[i]+dmm[i]));
            der[1][8][1][i] = -s3*(dl[i]*nn+l*dnn[i]);
        }
        ind[1][8][0] = 3;
        ind[1][8][1] = 4;

        mat[2][4][0] = s3*mm*l;
        mat[2][4][1] = l*(1-2*mm);
        for (i = 0; i < 3; ++i) {
            der[2][4][0][i] = s3*(dmm[i]*l+mm*dl[i]);
            der[2][4][1][i] = dl[i]*(1-2*mm)+l*(-2*dmm[i]);
        }
        ind[2][4][0] = 3;
        ind[2][4][1] = 4;
                
        mat[2][5][0] = s3*mm*n;
        mat[2][5][1] = n*(1-2*mm);
        for (i = 0; i < 3; ++i) {
            der[2][5][0][i] = s3*(dmm[i]*n+mm*dn[i]);
            der[2][5][1][i] = dn[i]*(1-2*mm)+n*(-2*dmm[i]);
        }
        ind[2][5][0] = 3;
        ind[2][5][1] = 4;
                
        mat[2][6][0] = s3*m*n*l;
        mat[2][6][1] = -2*m*n*l;
        for (i = 0; i < 3; ++i) {
            der[2][6][0][i] = s3*(dm[i]*n*l+m*dn[i]*l+m*n*dl[i]);
            der[2][6][1][i] = -2*(dm[i]*n*l+m*dn[i]*l+m*n*dl[i]);
        }
        ind[2][6][0] = 3;
        ind[2][6][1] = 4;
                
        mat[2][7][0] = 0.5*s3*m*(ll-mm);
        mat[2][7][1] = -m*(1+ll-mm);
        for (i = 0; i < 3; ++i) {
            der[2][7][0][i] = 0.5*s3*(dm[i]*(ll-mm)+m*(dll[i]-dmm[i]));
            der[2][7][1][i] = -(dm[i]*(1+ll-mm)+m*(dll[i]-dmm[i]));
        }
        ind[2][7][0] = 3;
        ind[2][7][1] = 4;
                
        mat[2][8][0] = m*(nn-0.5*(ll+mm));
        mat[2][8][1] = -s3*m*nn;
        for (i = 0; i < 3; ++i) {
            der[2][8][0][i] = dm[i]*(nn-0.5*(ll+mm))+m*(dnn[i]-0.5*(dll[i]+dmm[i]));
            der[2][8][1][i] = -s3*(dm[i]*nn+m*dnn[i]);
        }
        ind[2][8][0] = 3;
        ind[2][8][1] = 4;

        mat[3][4][0] = s3*l*m*n;
        mat[3][4][1] = -2*m*n*l;
        for (i = 0; i < 3; ++i) {
            der[3][4][0][i] = s3*(dl[i]*m*n+l*dm[i]*n+l*m*dn[i]);
            der[3][4][1][i] = -2*(dm[i]*n*l+m*dn[i]*l+m*n*dl[i]);
        }
        ind[3][4][0] = 3;
        ind[3][4][1] = 4;
                
        mat[3][5][0] = s3*nn*m;
        mat[3][5][1] = m*(1-2*nn);
        for (i = 0; i < 3; ++i) {
            der[3][5][0][i] = s3*(dnn[i]*m+nn*dm[i]);
            der[3][5][1][i] = dm[i]*(1-2*nn)+m*(-2*dnn[i]);
        }
        ind[3][5][0] = 3;
        ind[3][5][1] = 4;
                
        mat[3][6][0] = s3*nn*l;
        mat[3][6][1] = l*(1-2*nn);
        for (i = 0; i < 3; ++i) {
            der[3][6][0][i] = s3*(dnn[i]*l+nn*dl[i]);
            der[3][6][1][i] = dl[i]*(1-2*nn)+l*(-2*dnn[i]);
        }
        ind[3][6][0] = 3;
        ind[3][6][1] = 4;
                
        mat[3][7][0] = 0.5*s3*n*(ll-mm);
        mat[3][7][1] = -n*(ll-mm);
        for (i = 0; i < 3; ++i) {
            der[3][7][0][i] = 0.5*s3*(dn[i]*(ll-mm)+n*(dll[i]-dmm[i]));
            der[3][7][1][i] = -(dn[i]*(ll-mm)+n*(dll[i]-dmm[i]));
        }
        ind[3][7][0] = 3;
        ind[3][7][1] = 4;
                
        mat[3][8][0] = n*(nn-0.5*(ll+mm));
        mat[3][8][1] = s3*n*(ll+mm);
        for (i = 0; i < 3; ++i) {
            der[3][8][0][i] = dn[i]*(nn-0.5*(ll+mm))+n*(dnn[i]-0.5*(dll[i]+dmm[i]));
            der[3][8][1][i] = s3*(dn[i]*(ll+mm)+n*(dll[i]+dmm[i]));
        }
        ind[3][8][0] = 3;
        ind[3][8][1] = 4;

        mat[4][4][0] = 3*ll*mm;
        mat[4][4][1] = ll+mm-4*ll*mm;
        mat[4][4][2] = nn+ll*mm;
        for (i = 0; i < 3; ++i) {
            der[4][4][0][i] = 3*(dll[i]*mm+ll*dmm[i]);
            der[4][4][1][i] = dll[i]+dmm[i]-4*(dll[i]*mm+ll*dmm[i]);
            der[4][4][2][i] = dnn[i]+(dll[i]*mm+ll*dmm[i]);
        }
        ind[4][4][0] = 0;
        ind[4][4][1] = 1;
        ind[4][4][2] = 2;
            
        mat[4][5][0] = 3*l*mm*n;
        mat[4][5][1] = l*n*(1-4*mm);
        mat[4][5][2] = l*n*(mm-1);
        for (i = 0; i < 3; ++i) {
            der[4][5][0][i] = 3*(dl[i]*mm*n+l*dmm[i]*n+l*mm*dn[i]);
            der[4][5][1][i] = dl[i]*n*(1-4*mm)+l*dn[i]*(1-4*mm)+l*n*(-4*dmm[i]);
            der[4][5][2][i] = dl[i]*n*(mm-1)+l*dn[i]*(mm-1)+l*n*(dmm[i]);
        }
        ind[4][5][0] = 0;
        ind[4][5][1] = 1;
        ind[4][5][2] = 2;
            
        mat[4][6][0] = 3*ll*m*n;
        mat[4][6][1] = m*n*(1-4*ll);
        mat[4][6][2] = m*n*(ll-1);
        for (i = 0; i < 3; ++i) {
            der[4][6][0][i] = 3*(dll[i]*m*n+ll*dm[i]*n+ll*m*dn[i]);
            der[4][6][1][i] = dm[i]*n*(1-4*ll)+m*dn[i]*(1-4*ll)+m*n*(-4*dll[i]);
            der[4][6][2][i] = dm[i]*n*(ll-1)+m*dn[i]*(ll-1)+m*n*(dll[i]);
        }
        ind[4][6][0] = 0;
        ind[4][6][1] = 1;
        ind[4][6][2] = 2;
            
        mat[4][7][0] = 1.5*l*m*(ll-mm);
        mat[4][7][1] = 2*l*m*(mm-ll);
        mat[4][7][2] = 0.5*l*m*(ll-mm);
        for (i = 0; i < 3; ++i) {
            der[4][7][0][i] = 1.5*(dl[i]*m*(ll-mm)+l*dm[i]*(ll-mm)+l*m*(dll[i]-dmm[i]));
            der[4][7][1][i] = 2*(dl[i]*m*(mm-ll)+l*dm[i]*(mm-ll)+l*m*(dmm[i]-dll[i]));
            der[4][7][2][i] = 0.5*(dl[i]*m*(ll-mm)+l*dm[i]*(ll-mm)+l*m*(dll[i]-dmm[i]));
        }
        ind[4][7][0] = 0;
        ind[4][7][1] = 1;
        ind[4][7][2] = 2;
            
        mat[4][8][0] = s3*l*m*(nn-0.5*(ll+mm));
        mat[4][8][1] = - 2*s3*l*m*nn;
        mat[4][8][2] = 0.5*s3*l*m*(1+nn);
        for (i = 0; i < 3; ++i) {
            der[4][8][0][i] = s3*( dl[i]*m*(nn-0.5*(ll+mm))+l*dm[i]*(nn-0.5*(ll+mm))+l*m*(dnn[i]-0.5*(dll[i]+dmm[i])) );
            der[4][8][1][i] = -2*s3*(dl[i]*m*nn+l*dm[i]*nn+l*m*dnn[i]);
            der[4][8][2][i] = 0.5*s3*( dl[i]*m*(1+nn)+l*dm[i]*(1+nn)+l*m*(dnn[i]) );
        }
        ind[4][8][0] = 0;
        ind[4][8][1] = 1;
        ind[4][8][2] = 2;
            
        mat[5][5][0] = 3*mm*nn;
        mat[5][5][1] = (mm+nn-4*mm*nn);
        mat[5][5][2] = (ll+mm*nn);
        for (i = 0; i < 3; ++i) {
            der[5][5][0][i] = 3*(dmm[i]*nn+mm*dnn[i]);
            der[5][5][1][i] = (dmm[i]+dnn[i]-4*(dmm[i]*nn+mm*dnn[i]));
            der[5][5][2][i] = (dll[i]+dmm[i]*nn+mm*dnn[i]);
        }
        ind[5][5][0] = 0;
        ind[5][5][1] = 1;
        ind[5][5][2] = 2;
            
        mat[5][6][0] = 3*m*nn*l;
        mat[5][6][1] = m*l*(1-4*nn);
        mat[5][6][2] = m*l*(nn-1);
        for (i = 0; i < 3; ++i) {
            der[5][6][0][i] = 3*(dm[i]*nn*l+m*dnn[i]*l+m*nn*dl[i]);
            der[5][6][1][i] = dm[i]*l*(1-4*nn)+m*dl[i]*(1-4*nn)+m*l*(-4*dnn[i]);
            der[5][6][2][i] = dm[i]*l*(nn-1)+m*dl[i]*(nn-1)+m*l*(dnn[i]);
        }
        ind[5][6][0] = 0;
        ind[5][6][1] = 1;
        ind[5][6][2] = 2;
            
        mat[5][7][0] = 1.5*m*n*(ll-mm);
        mat[5][7][1] = - m*n*(1+2*(ll-mm));
        mat[5][7][2] = m*n*(1+0.5*(ll-mm));
        for (i = 0; i < 3; ++i) {
            der[5][7][0][i] = 1.5*( dm[i]*n*(ll-mm)+m*dn[i]*(ll-mm)+m*n*(dll[i]-dmm[i]) );
            der[5][7][1][i] = - ( dm[i]*n*(1+2*(ll-mm))+m*dn[i]*(1+2*(ll-mm))+m*n*(2*dll[i]-2*dmm[i]) );
            der[5][7][2][i] = dm[i]*n*(1+0.5*(ll-mm))+m*dn[i]*(1+0.5*(ll-mm))+m*n*(0.5*(dll[i]-dmm[i]));
        }
        ind[5][7][0] = 0;
        ind[5][7][1] = 1;
        ind[5][7][2] = 2;
            
        mat[5][8][0] = s3*m*n*(nn-0.5*(ll+mm));
        mat[5][8][1] = s3*m*n*(ll+mm-nn);
        mat[5][8][2] = -0.5*s3*m*n*(ll+mm);
        for (i = 0; i < 3; ++i) {
            der[5][8][0][i] = s3*( dm[i]*n*(nn-0.5*(ll+mm)) + m*dn[i]*(nn-0.5*(ll+mm))+m*n*(dnn[i]-0.5*(dll[i]+dmm[i])) );
            der[5][8][1][i] = s3*( dm[i]*n*(ll+mm-nn)+m*dn[i]*(ll+mm-nn)+m*n*(dll[i]+dmm[i]-dnn[i]) );
            der[5][8][2][i] = - 0.5*s3*( dm[i]*n*(ll+mm)+m*dn[i]*(ll+mm)+m*n*(dll[i]+dmm[i]) );
        }
        ind[5][8][0] = 0;
        ind[5][8][1] = 1;
        ind[5][8][2] = 2;
            
        mat[6][6][0] = 3*nn*ll;
        mat[6][6][1] = (nn+ll-4*nn*ll);
        mat[6][6][2] = (mm+nn*ll);
        for (i = 0; i < 3; ++i) {
            der[6][6][0][i] = 3*(dnn[i]*ll+nn*dll[i]);
            der[6][6][1][i] = dnn[i]+dll[i]-4*(dnn[i]*ll+nn*dll[i]);
            der[6][6][2][i] = (dmm[i]+dnn[i]*ll+nn*dll[i]);
        }
        ind[6][6][0] = 0;
        ind[6][6][1] = 1;
        ind[6][6][2] = 2;
            
        mat[6][7][0] = 1.5*n*l*(ll-mm);
        mat[6][7][1] = n*l*(1-2*(ll-mm));
        mat[6][7][2] = - n*l*(1-0.5*(ll-mm));
        for (i = 0; i < 3; ++i) {
            der[6][7][0][i] = 1.5*( dn[i]*l*(ll-mm)+n*dl[i]*(ll-mm)+n*l*(dll[i]-dmm[i]) );
            der[6][7][1][i] = dn[i]*l*(1-2*(ll-mm))+n*dl[i]*(1-2*(ll-mm))+n*l*(-2*(dll[i]-dmm[i]));
            der[6][7][2][i] = -( dn[i]*l*(1-0.5*(ll-mm))+n*dl[i]*(1-0.5*(ll-mm))+n*l*(-0.5*(dll[i]-dmm[i])) );
        }
        ind[6][7][0] = 0;
        ind[6][7][1] = 1;
        ind[6][7][2] = 2;
                
        mat[6][8][0] = s3*l*n*(nn-0.5*(ll+mm));
        mat[6][8][1] = s3*l*n*(ll+mm-nn);
        mat[6][8][2] = - 0.5*s3*l*n*(ll+mm);
        for (i = 0; i < 3; ++i) {
            der[6][8][0][i] = s3*( dl[i]*n*(nn-0.5*(ll+mm))+l*dn[i]*(nn-0.5*(ll+mm))+l*n*(dnn[i]-0.5*(dll[i]+dmm[i])) );
            der[6][8][1][i] = s3*( dl[i]*n*(ll+mm-nn)+l*dn[i]*(ll+mm-nn)+l*n*(dll[i]+dmm[i]-dnn[i]) );
            der[6][8][2][i] = - 0.5*s3*( dl[i]*n*(ll+mm)+l*dn[i]*(ll+mm)+l*n*(dll[i]+dmm[i]) );
        }
        ind[6][8][0] = 0;
        ind[6][8][1] = 1;
        ind[6][8][2] = 2;
    
        mat[7][7][0] = 0.75*pow2(ll-mm);
        mat[7][7][1] = (ll+mm-pow2(ll-mm));
        mat[7][7][2] = (nn+0.25*pow2(ll-mm));
        for (i = 0; i < 3; ++i) {
            der[7][7][0][i] = 0.75*2*(ll-mm)*(dll[i]-dmm[i]);
            der[7][7][1][i] = (dll[i]+dmm[i]-2*(ll-mm)*(dll[i]-dmm[i]));
            der[7][7][2][i] = (dnn[i]+0.25*2*(ll-mm)*(dll[i]-dmm[i]));
        }
        ind[7][7][0] = 0;
        ind[7][7][1] = 1;
        ind[7][7][2] = 2;
    
        mat[7][8][0] = 0.5*s3*(ll-mm)*(nn-0.5*(ll+mm));
        mat[7][8][1] = s3*nn*(mm-ll);
        mat[7][8][2] = 0.25*s3*(1+nn)*(ll-mm);
        for (i = 0; i < 3; ++i) {
            der[7][8][0][i] = 0.5*s3*( (dll[i]-dmm[i])*(nn-0.5*(ll+mm))+(ll-mm)*(dnn[i]-0.5*(dll[i]+dmm[i])) );
            der[7][8][1][i] = s3*( dnn[i]*(mm-ll)+nn*(dmm[i]-dll[i]) );
            der[7][8][2][i] = 0.25*s3*( dnn[i]*(ll-mm)+(1+nn)*(dll[i]-dmm[i]) );
        }
        ind[7][8][0] = 0;
        ind[7][8][1] = 1;
        ind[7][8][2] = 2;
                            
        mat[8][8][0] = pow2(nn-0.5*(ll+mm));
        mat[8][8][1] = 3*nn*(ll+mm);
        mat[8][8][2] = 0.75*pow2(ll+mm);
        for (i = 0; i < 3; ++i) {
            der[8][8][0][i] = 2*(nn-0.5*(ll+mm))*(dnn[i]-0.5*(dll[i]+dmm[i]));
            der[8][8][1][i] = 3*( dnn[i]*(ll+mm)+nn*(dll[i]+dmm[i]) );
            der[8][8][2][i] = 0.75*2*(ll+mm)*(dll[i]+dmm[i]);
        }
        ind[8][8][0] = 0;
        ind[8][8][1] = 1;
        ind[8][8][2] = 2;
    }
        
    /* use the same rules for orbitals when they are reversed (pd ->dp)... */
    for (a = 0; a < maxorb; ++a) {
        for (b = a+1; b < maxorb; ++b) {
            for (i = 0; i < 3; ++i) {
                mat[b][a][i] = mat[a][b][i];
                for (j = 0; j < 3; ++j) {
                    der[b][a][i][j] = der[a][b][i][j];
                }
                ind[b][a][i] = ind[a][b][i];
            }
        }
    }

    /* ...but use different indices from table            
       pd 3:5-->10:12 
       sd 7->12
       sp 8->13 */

    ind[1][0][0] = 13;
    ind[2][0][0] = 13;
    ind[3][0][0] = 13;

    ind[4][0][0] = 12;
    ind[5][0][0] = 12;
    ind[6][0][0] = 12;
    ind[7][0][0] = 12;
    ind[8][0][0] = 12;

    for (i = 0; i < 2; ++i) {
        ind[4][1][i] = 10+i;
        ind[5][1][i] = 10+i; 
        ind[6][1][i] = 10+i;
        ind[7][1][i] = 10+i;
        ind[8][1][i] = 10+i; 
        ind[4][2][i] = 10+i;
        ind[5][2][i] = 10+i; 
        ind[6][2][i] = 10+i;
        ind[7][2][i] = 10+i;
        ind[8][2][i] = 10+i;
        ind[4][3][i] = 10+i;
        ind[5][3][i] = 10+i;
        ind[6][3][i] = 10+i;
        ind[7][3][i] = 10+i;
        ind[8][3][i] = 10+i;
    }
        
    for (i = 0; i < noi; ++i) {
        for (j = 0; j < noj; ++j) {
            nc = cnt[i][j] - 1;

            ht[noj*i+j] = 0.0;
            st[noj*i+j] = 0.0;
            for (a = 0; a < 3; ++a) {
                dht[3*(noj*i+j)+a] = 0.0;
                dst[3*(noj*i+j)+a] = 0.0;
            }

            for (b = 0; b <= nc; ++b) {

                ht[noj*i+j] += mat[i][j][b]*h[ind[i][j][b]];
                st[noj*i+j] += mat[i][j][b]*s[ind[i][j][b]];

                for (a = 0; a < 3; ++a) {
                    dht[3*(noj*i+j)+a] += mat[i][j][b]*dh[3*ind[i][j][b]+a] + der[i][j][b][a]*h[ind[i][j][b]];
                    dst[3*(noj*i+j)+a] += mat[i][j][b]*ds[3*ind[i][j][b]+a] + der[i][j][b][a]*s[ind[i][j][b]];
                }

            }

        }
    }
}


PyObject *
py_fast_slako_transformations(PyObject *self, PyObject *args)
{
    int noi, noj;
    double dist;
    PyObject *rhat, *h, *s, *dh, *ds;
    PyObject *ht, *st, *dht, *dst;

    npy_intp dims[3];
    double *ht_data, *st_data, *dht_data, *dst_data;

    if (!PyArg_ParseTuple(args, "O!diiO!O!O!O!", 
                          &PyArray_Type, &rhat,
                          &dist, &noi, &noj,
                          &PyArray_Type, &h,
                          &PyArray_Type, &s,
                          &PyArray_Type, &dh,
                          &PyArray_Type, &ds))
        return NULL;


    /*
     * Error checking
     */

    if (noi != 1 && noi != 4 && noi != 9) {
        PyErr_SetString(PyExc_RuntimeError, "noi can only be 1, 4 or 9.");
        return NULL;
    }

    if (noj != 1 && noj != 4 && noj != 9) {
        PyErr_SetString(PyExc_RuntimeError, "noj can only be 1, 4 or 9.");
        return NULL;
    }

    if (PyArray_NDIM(rhat) != 1 || PyArray_DIM(rhat, 0) != 3) {
        PyErr_SetString(PyExc_TypeError, "rhat needs to be a 3-vector.");
        return NULL;
    }

    if (PyArray_NDIM(h) != 1 || PyArray_DIM(h, 0) != 14) {
        PyErr_SetString(PyExc_TypeError, "h needs to be a 14-vector.");
        return NULL;
    }

    if (PyArray_NDIM(s) != 1 || PyArray_DIM(s, 0) != 14) {
        PyErr_SetString(PyExc_TypeError, "s needs to be a 14-vector.");
        return NULL;
    }

    if (PyArray_NDIM(dh) != 2 || PyArray_DIM(dh, 0) != 14 || PyArray_DIM(dh, 1) != 3) {
        PyErr_SetString(PyExc_TypeError, "dh needs to be a 14 x 3 matrix.");
        return NULL;
    }

    if (PyArray_NDIM(ds) != 2 || PyArray_DIM(ds, 0) != 14 || PyArray_DIM(ds, 1) != 3) {
        PyErr_SetString(PyExc_TypeError, "ds needs to be a 14 x 3 matrix.");
        return NULL;
    }

    ht_data = (double *) malloc(noi*noj*sizeof(double));
    st_data = (double *) malloc(noi*noj*sizeof(double));
    dht_data = (double *) malloc(3*noi*noj*sizeof(double));
    dst_data = (double *) malloc(3*noi*noj*sizeof(double));

    fast_slako_transformations(PyArray_DATA(rhat),
                               dist, noi, noj,
                               PyArray_DATA(h),
                               PyArray_DATA(s),
                               PyArray_DATA(dh),
                               PyArray_DATA(ds),
                               ht_data,
                               st_data,
                               dht_data,
                               dst_data);

    dims[0] = noi;
    dims[1] = noj;
    dims[2] = 3;

    ht = PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, NULL, ht_data, NPY_OWNDATA, 0, NULL);
    st = PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, NULL, st_data, NPY_OWNDATA, 0, NULL);

    dht = PyArray_New(&PyArray_Type, 3, dims, NPY_DOUBLE, NULL, dht_data, NPY_OWNDATA, 0, NULL);
    dst = PyArray_New(&PyArray_Type, 3, dims, NPY_DOUBLE, NULL, dst_data, NPY_OWNDATA, 0, NULL);

    return Py_BuildValue("OOOO", ht, st, dht, dst);
}
