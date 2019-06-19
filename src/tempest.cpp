////////////////////////////////////////////////////////////////////////////////////////////////////
//
// tempest 
//
// This code is an implementation of "Simulating Ocean Water" by Jerry Tessendorf
//
// Doug Creel
// Imagemovers Digital (2008)
//
////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <stdlib.h>
#include <fft.h>
#include <Magick++.h>


using namespace std;
using namespace Magick;

#define GRAVITY 9.81
#define ACONST 0.0008
#define INVSQRT2 0.70710678118654746
#define BUFFER 2048
#define Complex complex<double>


// this is the Phillips spectrum term
double phillips(double a, double k[2], double wind[2])
{
    double k2, km, v2, l, ret;
    k2 = k[0]*k[0]+k[1]*k[1];
    if (k2==0)
        return 0;
    km = sqrt(k2);
    v2 = wind[0]*wind[0]+wind[1]*wind[1];
    l = v2 / GRAVITY;
    // the factor exp(-sqrt(k2)*1.0) is to get rid of small waves
    //ret = a*(exp(-1/(k2*l*l))/(k2*k2))*((k[0]*wind[0]+k[1]*wind[1])*(k[0]*wind[0]+k[1]*wind[1]) / (k2*v2))*exp(-sqrt(k2)*1.0);
    ret = a*(exp(-1/(k2*l*l))/(k2*k2))*((k[0]*wind[0]+k[1]*wind[1])*(k[0]*wind[0]+k[1]*wind[1]) / (k2*v2));
    return ret;
}

int main(int argc, char* argv[])
{
    int dir, nx, ny, lx, ly;
    double ctime;
    double wind[2];
    char outfile[BUFFER];

    if (argc != 14) {
        cout << "Usage:  tempest -r <xres yres>  -l <xsize ysize> -w <xspeed yspeed> -t <time> -o <output file>" << endl;
        return(0);
    }

    for (int i=1; i<argc; i++) {
        if (!strcmp(argv[i], "-r")) {nx = atoi(argv[++i]); ny = atoi(argv[++i]);}
        if (!strcmp(argv[i], "-l")) {lx = atoi(argv[++i]); ly = atoi(argv[++i]);}
        if (!strcmp(argv[i], "-w")) {wind[0] = atof(argv[++i]); wind[1] = atof(argv[++i]);}
        if (!strcmp(argv[i], "-t")) ctime = atof(argv[++i]);
        if (!strcmp(argv[i], "-o")) sprintf(outfile, "%s", argv[++i]);
    }

    // define new image
    Image image(Geometry(nx,ny),Color(0, 0, 0, 0));

    // assign memory
    Complex** ht0;
    Complex** fftdata;
    fftdata = new Complex*[nx];
    ht0 = new Complex*[nx];
    for (int i=0; i<nx; i++) {
        fftdata[i] = new Complex[ny];
        ht0[i] = new Complex[ny];
    }

    double kvec[2], gauss[2]; 
    unsigned int seed = 0;
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            kvec[0] = 2.0*M_PI*((double)i - 0.5*nx)/lx;
            kvec[1] = 2.0*M_PI*((double)j - 0.5*ny)/ly;
            Complex zeta(((double)rand_r(&seed)/RAND_MAX,  (double)rand_r(&seed)/RAND_MAX));
            ht0[i][j] = INVSQRT2*sqrt(phillips(ACONST, kvec, wind))*zeta;
        }
    }

    double km, wkt;
    double speed = 1./24.;
    // we only need to loop over half the values since we can calculate the rest using conjugate
    int xhalf = nx/2 + 1;
    for (int i=0; i<xhalf; ++i) {
        for (int j=0; j<ny; ++j) {
            kvec[0] = 2.0*M_PI*((double)i - 0.5*nx)/lx;
            kvec[1] = 2.0*M_PI*((double)j - 0.5*ny)/ly;
            km = sqrt(kvec[0]*kvec[0] + kvec[1]*kvec[1]);
            wkt = sqrt(km*GRAVITY)*ctime*speed;

            Complex ep(cos(wkt), sin(wkt)); 
            Complex em = std::conj(ep); 
            Complex ht0star = std::conj(ht0[nx-i-1][j]);
            fftdata[i][j] = ht0[i][j]*ep + ht0star*em;

            // h~(-K) = conj(h~(K))
            if (i != xhalf-1) {
                fftdata[nx-i-1][j] = -std::conj(fftdata[i][j]);
            }
        }
    }

    // run fft 
    dir = -1;
    if (!FFT2D(fftdata, nx, ny, dir)) {
        cout << "FFT2D returned error" << endl;
        exit(0);
    }

    // write out data
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            double value = fftdata[i][j].real() * powf(-1,i+j) / lx;
            value = (1 + value)*0.5;
            image.pixelColor(i,j,Color((int)(MaxRGB*value+0.5), (int)(MaxRGB*value+0.5), (int)(MaxRGB*value+0.5), MaxRGB));
        }
    }

    image.write(outfile);

    return 0;
}


