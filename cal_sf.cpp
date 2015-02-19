#include <stdio.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <fftw3.h>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <string>
using namespace std;

//#define OUPUT3DSF
//#define DUMPXYZ // yes, shrinking works as intended

double inline sqr( double a ) { return a*a; }
const int nh = 50;
const double drm = double(RAND_MAX)+1.0;

int main( int argc, char *argv[] )
{
  if( argc < 2 || argc > 4 ) {
    cerr << "./structurefactor file.xyz [ss [shrinkfactor]]" << endl;
    return 1;
  }
  int shrink = (argc==4) ? atoi(argv[3]) : 1;
  int hb, sm, ss = (argc>=3) ?  atoi(argv[2]) : 64;
  int nhisto[nh], i, j, x[3], x2[3];
  
  //if shrink, test if ss is power of 2 
  if( shrink > 1 ) {
    for( i = 0; i < 64 && ss>>i != 1; i++ ); 
    cout << i << ' ' << (1<<i) << endl;
    if( (1<<i) != ss ) {
      cerr << "illegal cell size (must be power of 2)" << endl;
      return 1;
    }
    hb = 1 << (i-1); // high bit carryover
  }

  int na = ss*ss*ss;
  double hbin = 0.6*double(ss)/double(nh);
  double histo[nh], x3[3];
  int skip = 0; // skip this many frames
  bool sfaverage = true;

  if( sfaverage )
  {
    memset( histo, 0, nh * sizeof(double) );
    memset( nhisto, 0, nh * sizeof(int) );
  }

  fftw_complex* f3dmap = (fftw_complex*)fftw_malloc( ss * ss * ss * sizeof(fftw_complex) );
  fftw_plan     plan3d = fftw_plan_dft_3d( ss, ss, ss, f3dmap, f3dmap, FFTW_FORWARD, FFTW_ESTIMATE );
  unsigned char *rbin = (unsigned char*)malloc( sizeof(unsigned char) * na );

  // initialize bin lookup table
  for( x[0] = 0; x[0] < ss; x[0]++ )
    for( x[1] = 0; x[1] < ss; x[1]++ )
      for( x[2] = 0; x[2] < ss; x[2]++ )
      {
        double mind2 = -1, d2, x3[3];
        // iterate over neighboring periodic copies
        for( x2[0] = -1; x2[0] < 1; x2[0]++ )
          for( x2[1] = -1; x2[1] < 1; x2[1]++ )
            for( x2[2] = -1; x2[2] < 1; x2[2]++ )
            {
              d2 = 0;
              for( i = 0; i < 3; ++i )
                x3[i] = x2[i]*ss + x[i];
              for( i = 0; i < 3; ++i )
                d2 += 0.25 * sqr( -x3[i] + x3[(i+1)%3] + x3[(i+2)%3] );

              if( mind2 < 0 || d2 < mind2 ) mind2 = d2;
            }
        rbin[x[0]+ss*(x[1]+ss*x[2])] = (unsigned char)(sqrt(mind2)/hbin);
      }

  FILE *in = fopen( argv[1], "rt" );
#ifdef DUMPXYZ
  FILE *out = fopen( "dump.xyz", "wt" );
#endif
  int s;
  for( s = 0; ; ++s )
  {
    int nb, a,b,c;
    char buffer[1000], el;

    fgets( buffer, 1000, in);
    sscanf( buffer, "%d\n", &nb );
    if( feof( in ) ) break;

    // read b atoms
    fgets( buffer, 1000, in);

    double cb = double(nb-1)/double(na); // minus 1 vacancy
    //fprintf( stderr, "Reading frame %d (cb=%f, nb=%d)\n", s+1, float(cb),nb );
    if( s >= skip )
    {
      // initialize
      for( x[0] = 0; x[0] < ss; x[0]++ )
        for( x[1] = 0; x[1] < ss; x[1]++ )
          for( x[2] = 0; x[2] < ss; x[2]++ )
          {
            f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] = -cb;
            f3dmap[x[0]+ss*(x[1]+ss*x[2])][1] = 0.0;
          }
    }

    if( !sfaverage )
    {
      memset( histo, 0, nh * sizeof(double) );
      memset( nhisto, 0, nh * sizeof(int) );
    }

    for( i = 0; i < nb; ++i )
    {
      fgets( buffer, 1000, in);
      sscanf( buffer, "%c %d %d %d %*d %*d\n", &el, &a, &b, &c );
//      if( el == 'B' )
//      {
        // project from realspace back into fcc-basis coordinates
//        x[0] = ( b + c - a ) / 2;
//        x[1] = ( c + a - b ) / 2;
//        x[2] = ( a + b - c ) / 2;

        // now shrink coordinates
        for( j = 1; j < shrink; ++j ) {
          for( int k = 0; k < 3; ++k ) {
            x[k] = (x[k]>>1) + (x[k]&1)*hb; // rotate with wraping
          }
        }
        f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] = 1.0-cb;
//      }
    }

    #ifdef DUMPXYZ
      fprintf( out, "%d\n Frame number %lld %f fs boxsize %f %f %f\n", nb-1, s, float(s), float(ss), float(ss), float(ss) ); 
      for( x[0] = 0; x[0] < ss; x[0]++ )
        for( x[1] = 0; x[1] < ss; x[1]++ )
          for( x[2] = 0; x[2] < ss; x[2]++ )
            if( f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] == 1.0-cb ) {
              fprintf( out, "B %d %d %d\n", x[1]+x[2], x[2]+x[0], x[0]+x[1] );
            }
    #endif

    // insert sphere for calibration
    //#define CALIBRATE
    #ifdef CALIBRATE
    const int particle_radius = 16/2;
    nb = 0;
    memset( f3dmap, 0, na * sizeof(fftw_complex) );
    for( x[0] = 0; x[0] < ss; x[0]++ )
      for( x[1] = 0; x[1] < ss; x[1]++ )
        for( x[2] = 0; x[2] < ss; x[2]++ )
        {
          a = x[1]+x[2];
          b = x[0]+x[2];
          c = x[1]+x[0];
          if( sqr(a-32)+sqr(b-32)+sqr(c-32) <= particle_radius*particle_radius )
          {
            f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] = 1.0;
            nb++;
          }
          else 
            f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] = 0.0;
        }
    cb = double(nb)/double(na);
    for( x[0] = 0; x[0] < ss; x[0]++ )
      for( x[1] = 0; x[1] < ss; x[1]++ )
        for( x[2] = 0; x[2] < ss; x[2]++ )
          f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] -= cb;
    #endif

    //#define RANDOMCALIBRATE
    #ifdef RANDOMCALIBRATE
    // generate a random permutation set and initialize id table
    cb = 0.15;
    int *perm = (int*)calloc( na, sizeof(int) );
    int swap;
    nb = int(double(na)*cb);
    for( i = 0; i < na; i++ ) perm[i] = i;
    for( i = 0; i < nb; i++ )
    {
      j = int( (double(rand())/drm) * (na-i) ) + i;
      swap = perm[j]; perm[j] = perm[i]; perm[i] = swap;
    }
    for( i = 0; i < na; i++ )
      f3dmap[perm[i]][0] = i<nb ? 1-cb : -cb;
    free(perm);
    #endif

    // insert checkerboard pattern for calibration
    //#define CHECKCALIBRATE
    #ifdef CHECKCALIBRATE
    const int checker_size = 8;
    nb = na/2;
    memset( f3dmap, 0, na * sizeof(fftw_complex) );
    for( x[0] = 0; x[0] < ss; x[0]++ )
      for( x[1] = 0; x[1] < ss; x[1]++ )
        for( x[2] = 0; x[2] < ss; x[2]++ )
          if( (x[0]/checker_size + x[1]/checker_size + x[2]/checker_size )%2 )
            f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] = -0.5;
          else 
            f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] = 0.5;
    printf( "%d\n\n", nb );
    for( x[0] = 0; x[0] < ss; x[0]++ )
      for( x[1] = 0; x[1] < ss; x[1]++ )
        for( x[2] = 0; x[2] < ss; x[2]++ )
        {
          if( f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] > 0.0 )
          {
            printf( "B %d %d %d\n", x[1]+x[2], x[0]+x[2], x[1]+x[0] );
          }
        }
    return 0;
    #endif


    // sum
    /*
    double sum = 0.0;
    for( x[0] = 0; x[0] < ss; x[0]++ )
      for( x[1] = 0; x[1] < ss; x[1]++ )
        for( x[2] = 0; x[2] < ss; x[2]++ )
          sum += f3dmap[x[0]+ss*(x[1]+ss*x[2])][0];
    fprintf( stderr, "  sum = %e\n", sum );
    */

    if( s >= skip )
    {
      // execute transform
      fftw_execute( plan3d );

      // process the result
      double k2, k[3];
      #ifdef OUPUT3DSF
        printf( "%d\n\n", na );
      #endif
      for( x[0] = 0; x[0] < ss; x[0]++ )
        for( x[1] = 0; x[1] < ss; x[1]++ )
          for( x[2] = 0; x[2] < ss; x[2]++ )
          {
            for( i = 0; i < 3; ++i )
              k[i] = -x[i] + x[(i+1)%3] + x[(i+2)%3];

            i = rbin[x[0]+ss*(x[1]+ss*x[2])];
            if( i < nh )
            {
              histo[i] += sqr( f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] ) + sqr( f3dmap[x[0]+ss*(x[1]+ss*x[2])][1] );
              nhisto[i]++;
            }
            #ifdef OUPUT3DSF
            printf( "B %f %f %f  %f %d\n", k[0],k[1],k[2],
                    log( sqr( f3dmap[x[0]+ss*(x[1]+ss*x[2])][0] ) +
                    sqr( f3dmap[x[0]+ss*(x[1]+ss*x[2])][1] ) ), i );
            #endif
          }

      // output histogram
      #ifndef OUPUT3DSF
      if( !sfaverage )
      {
        for( i = 0; i < nh; ++i )
          printf( "%d %f %e\n", s, double(i)*hbin, 
                  nhisto[i]>0 ? double(i*i*hbin*hbin)*histo[i]/double(nhisto[i]) : 0.0 );
        printf( "\n" );
      }
      #endif
    }
  }

  // output averaged structure factor and first three moments
  double mom[3] = { 0.0 }, kint[3] = { 0.0 }; 
  int ii;
  if( sfaverage )
  {
    for( i = 0; i < nh; ++i )
    {
      printf( "%f %e\n", double(i)*hbin, nhisto[i]>0 ? histo[i]/double(nhisto[i]) : 0.0 );
      ii = 1;
      for( j = 0; j < 3; j++ )
      {
        mom[j]  += hbin * double(i*ii) * histo[i]/double(nhisto[i]==0?1:nhisto[i]);
        kint[j] +=                 ii  * histo[i]/double(nhisto[i]==0?1:nhisto[i]);
        ii *= i;
      }
    }
  }
  cout << "#MM";
  for( j = 0; j < 3; j++ )
    cout << ' ' << double(ss)*kint[j]/mom[j]; // == ss/(mom/kint)
  cout << endl;

  /*
  fftw_complex in1d[nh], out1d[nh];
  fftw_plan plan1d = fftw_plan_dft_1d(nh, in1d, out1d, FFTW_BACKWARD, FFTW_ESTIMATE); 
  for( i = 0; i < nh; ++i )
  {
    in1d[i][0] = nhisto[i]>0 ? histo[i]/double(nhisto[i]) : 0.0;
    in1d[i][1] = 0.0;
  }
  fftw_execute(plan1d);
  for( i = 0; i < nh/2; ++i )
    printf( "%e %e\n", i*hbin, out1d[i][0]*out1d[i][0]+out1d[i][1]*out1d[i][1] );
  fclose(in);
  */
}
