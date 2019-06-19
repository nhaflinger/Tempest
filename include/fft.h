////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This code is an implementation of a fast Fourier transform Cooley-Tukey algorithm
//
// Doug Creel
// Imagemovers Digital (2008)
//
// Derived from code by Paul Bourke 
// http://local.wasp.uwa.edu.au/~pbourke/other/dft
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <complex.h>

int FFT2D(complex<double>** c,int nx,int ny,int dir);
int FFT(int dir,int m,double *x,double *y);
int Powerof2(int n,int* m,int* twopm);
int DFT(int dir,int m,double *x1,double *y1);


using namespace std;

int FFT2D(complex<double>** c,int nx,int ny,int dir)
{
    int i,j;
    int m,twopm;
    double *real,*imag;

    /* Transform the rows */
    real = new double[nx];
    imag = new double[nx];
    if (real == NULL || imag == NULL)
        return(0);
    if (!Powerof2(nx,&m,&twopm) || twopm != nx) {
        return(0);
    }
    for (j=0;j<ny;j++) {
        for (i=0;i<nx;i++) {
            real[i] = c[i][j].real();
            imag[i] = c[i][j].imag();
        }
        FFT(dir,m,real,imag);
        for (i=0;i<nx;i++) {
            c[i][j].real() = real[i];
            c[i][j].imag() = imag[i];
        }
    }
    delete [] real;
    delete [] imag;

    /* Transform the columns */
    real = new double[ny];
    imag = new double[ny];
    if (real == NULL || imag == NULL)
        return(0);
    if (!Powerof2(ny,&m,&twopm) || twopm != ny) {
        return(0);
    }
    for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
            real[j] = c[i][j].real();
            imag[j] = c[i][j].imag();
        }
        FFT(dir,m,real,imag);
        for (j=0;j<ny;j++) {
            c[i][j].real() = real[j];
            c[i][j].imag() = imag[j];
        }
    }
    delete [] real;
    delete [] imag;

    return(1);
}

int FFT(int dir,int m,double *x,double *y)
{
    long nn,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;

    /* Calculate the number of points */
    nn = 1 << m;

    /* Do the bit reversal */
    i2 = nn >> 1;
    j = 0;
    for (i=0;i<nn-1;i++) {
        if (i < j) {
            tx = x[i];
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }

    /* Compute the FFT */
    c1 = -1.0;
    c2 = 0.0;
    l2 = 1;
    for (l=0;l<m;l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j=0;j<l1;j++) {
            for (i=j;i<nn;i+=l2) {
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            z =  u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }

        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }

    /* Scaling for forward transform */
    if (dir == 1) {
        for (i=0;i<nn;i++) {
            x[i] /= (double)nn;
            y[i] /= (double)nn;
        }
    }

    return(1);
}

int Powerof2(int n,int *m,int *twopm)
{
    if (n <= 1) {
        *m = 0;
        *twopm = 1;
        return(0);
    }

   *m = 1;
   *twopm = 2;
   do {
      (*m)++;
      (*twopm) *= 2;
   } while (2*(*twopm) <= n);

   if (*twopm != n)
       return(0);
   else
       return(1);
}

int DFT(int dir,int m,double *x1,double *y1)
{
    long i,k;
    double arg;
    double cosarg,sinarg;
    double *x2=NULL,*y2=NULL;
 
    x2 = (double*)malloc(m*sizeof(double));
    y2 = (double*)malloc(m*sizeof(double));
    if (x2 == NULL || y2 == NULL)
        return(0);

    for (i=0;i<m;i++) {
        x2[i] = 0;
        y2[i] = 0;
        arg = - dir * 2.0 * M_PI * (double)i / (double)m;
        for (k=0;k<m;k++) {
            cosarg = cos(k * arg);
            sinarg = sin(k * arg);
            x2[i] += (x1[k] * cosarg - y1[k] * sinarg);
            y2[i] += (x1[k] * sinarg + y1[k] * cosarg);
        }
    }

    /* Copy the data back */
    if (dir == 1) {
        for (i=0;i<m;i++) {
            x1[i] = x2[i] / (double)m;
            y1[i] = y2[i] / (double)m;
        }
    } 
    else {
        for (i=0;i<m;i++) {
            x1[i] = x2[i];
            y1[i] = y2[i];
        }
    }

    free(x2);
    free(y2);

    return(1);
}


