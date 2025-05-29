
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <unistd.h>
#include <algorithm>
using namespace std;


/* random number generator: Mersenne Twister **/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
gsl_rng_type *rndT = (gsl_rng_type *)gsl_rng_default;
gsl_rng *mt;


#define GNUPLOT_SHOW (1)      // make output into GNUPLOT or not
#define GNUPLOT_SAVE_GIF (0)  // write output onto gif_filename or not
char  gif_filename[] = "testanime.gif";
#define GNUPLOT_OUT_INTERVAL (10000) // make output every GNUPLOT_OUT_INTERVAL*Dt time



// Parameters
#define INITIAL_CONDITION_TYPE (1) //1: clustered condition, 0: uniformly-distributed condition
#define N  5000   // Number of particles
#define LX 100  // x-axis length
#define LY 100  // y-axis length

#define DOF (5*N)  // Degree of freedome to be solved
#define Dt 0.001    // disctetized time step
const int SIMULATION_TIME = 40000; //  = Dt*SIMULATION_STEP

// parameters of equatios
double gm = 10.0; // drag coeff.
double Temp = 0.0001;  // kBT
double AngTemp = 3*Temp; //  often used criteria, varidated for spherical colloid (?)
double RDD, AngS; // noise strength,
double f0 = 0.0; // active force, argv[1], indicated as 'f' in the paper.
double u0 = 0.0; // adhesion, argv[2], indicated as 'a' in the paper.
double tau = 1/(3*gm*Temp);// tau_theta, time-scale of change in theta

const double sigma0 = 1; // unit of length, interaction
const double sigma0sq = sigma0*sigma0;
const double Particle_size = pow(2.0,1.0/6.0)*sigma0;
const double Particle_sizeSq = Particle_size*Particle_size;


/** parametes used in lib/particle_neighbours.cpp **/
int h= 5;   // hx >= verlet_radius, divisible by LX,LY
int cell_size_x = LX/h, cell_size_y = LY/h;
double MARGIN = 0.5*sigma0;
double verlet_radius_sq = (2.5*sigma0+MARGIN)*(2.5*sigma0+MARGIN); // > 1.0 = Lb, interaction length

#include "lib/particle_neighbours.cpp"
#include "lib/adhABPlib.c"


void Accelelation_BAOAB( double *X, double *kF ){
  double *x  = X,  *y  = X+N, *theta = X+4*N;
  double *kfx = kF, *kfy = kF+N;

  for(int i=0; i<N; ++i )kfx[i] = f0*cos(theta[i]); // drag force
  for(int i=0; i<N; ++i )kfy[i] = f0*sin(theta[i]);

  for (const auto& itr : neighbours) {
    int i = itr.i, j = itr.j;
    double dx = x[j]-x[i] + itr.bx*LX;
    double dy = y[j]-y[i] + itr.by*LY;
    double dl = dx*dx+dy*dy;
    if( dl< Particle_sizeSq ){
      if(dl<0.3)dl=0.3;  // avoid to high interaction force
      double r6 = dl*dl*dl;
      double f = (r6-2.0)/(r6*r6*dl);
      kfx[i]+=f*dx;
      kfy[i]+=f*dy;
      kfx[j]-=f*dx;
      kfy[j]-=f*dy;
    }
    else{
      double r6 = dl*dl*dl;
      double f = u0*(r6-2.0)/(r6*r6*dl);
      if(dl>6.25*sigma0sq)f=0;  // avoid to high interaction force
      kfx[i]+=f*dx;
      kfy[i]+=f*dy;
      kfx[j]-=f*dx;
      kfy[j]-=f*dy;
    }
  }
}
    
/** BAOAB (improved Stochastic Velocity Verlet) **/
// https://pubs.rsc.org/en/content/articlepdf/2017/sm/c7sm01526g
const double HDt = Dt/2.0;
const double EDt = exp(-gm*Dt);
inline void SVV_BAOAV(double *X, double *kF){
  int i;
  double lns[2*N];
  //for( i=0; i<2*N; ++i )lns[i]= RDD_BAOAB*normal(mt);
  for( i=0; i<2*N; ++i )lns[i]= RDD*gsl_ran_gaussian_ziggurat(mt, 1.);
  for( i=0; i<2*N; ++i )X[i+2*N] += kF[i]*HDt;
  for( i=0; i<2*N; ++i )X[i] += HDt*X[i+2*N]; // update vel.
  for( i=0; i<2*N; ++i )X[i+2*N] = X[i+2*N]*EDt + lns[i];
  for( i=0; i<2*N; ++i )X[i] += HDt*X[i+2*N];
  Accelelation_BAOAB( X, kF );
  for( i=0; i<2*N; ++i )X[i+2*N] += kF[i]*HDt;
  for( i=0; i<N; ++i )X[4*N+i]+=AngS*gsl_ran_gaussian_ziggurat(mt, 1.);
}




////////////////////////////////////////////////
int main(int argc, char *argv[]){
  double X[DOF],kF[2*N];
  double *x = X, *y = X+N, *vx = X+2*N, *vy = X+3*N, *theta = X+4*N;
  double *kfx = kF,*kfy = kF+N;
 
  f0 = atof(argv[1]);
  u0 = atof(argv[2]);
 
  //parameter record
  cout << "#" << endl;
  cout << "#" << " " << "particle_num N =" << N << endl;
  cout << "#" << " " << "system_size LX = " << LX << endl;
  cout << "#" << " " << "system_size LY = " << LY << endl;
  cout << "#" << " " << "vol_fraction = " << (double)( (M_PI*Particle_sizeSq*N)/(4.0*LX*LY) )<< endl;
  cout << "#" << " " << "Dt = " <<  Dt << endl;
  cout << "#" << " " << "simulation time = " << SIMULATION_TIME << endl;
  cout << "#" << " " << "gamma = " << gm << endl;
  cout << "#" << " " << "kBT = " << Temp << endl;
  cout << "#" << " " << "tau = " << tau <<endl;
  cout << "#" << " " << "f0 = " << f0 << endl;
  cout << "#" << " " << "u0 = " << u0 << endl;
   

  //
  RDD = sqrt(Temp*(1.0-exp(-2.0*gm*Dt)));
  AngS = sqrt(2*AngTemp*gm*Dt);

  //omp_set_num_threads(NUMBER_THREADS);



  // random number initilization
  {
    gsl_rng_env_setup();
    mt = gsl_rng_alloc(rndT); // rand num generator
    //int seed = (unsigned long int) time(NULL);
    int seed = atoi(argv[3]);
    cout << "# seed = " << seed << endl;
    gsl_rng_set(mt, seed); // seed of rand
  }

  // initial conditions
    
  //for(int i=0; i<N; ++i)x[i] = (L)*gsl_rng_uniform(mt);
  //for(int i=0; i<N; ++i)y[i] = (L)*gsl_rng_uniform(mt);

  if(INITIAL_CONDITION_TYPE){
    double dcx = (double)LX/(double)(sqrt(N)+.0);
    double dcy = (double)LY/(double)(sqrt(N)+.0);
    for( int i=0; i<(int)sqrt(N)+1; i++){
      for(int j=0; j<(int)sqrt(N)+1; j++){
        if(i*((int)sqrt(N)+1)+j<N){
          x[i*((int)sqrt(N)+1)+j] = LX/4+(i+0.5)*dcx/2;
          y[i*((int)sqrt(N)+1)+j] = LY/4+(j+0.5)*dcy/2;
        }
      }
    }
  }
  else{
    double dcx = (double)LX/(double)(sqrt(N)+.0);
    double dcy = (double)LY/(double)(sqrt(N)+.0);
    for( int i=0; i<(int)(5*sqrt(N/5)+1); i++){
      for(int j=0; j<(int)(sqrt(N/5)+1); j++){
        if(i*((int)(sqrt(N)+1))+j<N){
          x[i*((int)(sqrt(N)+1))+j] = (i+0.5)*dcx;
          y[i*((int)(sqrt(N)+1))+j] = (j+0.5)*dcy;
        }
      }
    }
  }
    
    
    
  for(int i=0; i<N; ++i)vx[i] = 0.2*gsl_ran_gaussian_ziggurat(mt, 1.);
  for(int i=0; i<N; ++i)vy[i] = 0.2*gsl_ran_gaussian_ziggurat(mt, 1.);
  for(int i=0; i<N; ++i)theta[i] = 2*M_PI*gsl_rng_uniform(mt);
  makeneighbours( x, y );
  Accelelation_BAOAB( X, kF );

  // For output images
  FILE *gp = popen("gnuplot -persist","w");
  fprintf(gp, "unset key\n");
  fprintf(gp, "set xrange [0:%d]\n",LX/2);
  fprintf(gp, "set yrange [0:%d]\n",LY);
  fprintf(gp, "set size 1.0, 0.5\n");
  //fprintf(gp, "set size square\n");
  fprintf(gp, "set style fill solid 0.2\n");
  #if GNUPLOT_SAVE_GIF
    fprintf(gp, "set term gif animate\n");	
    fprintf(gp, "set output '%s'\n",gif_filename);
  #endif

  for(int tc=0; tc*Dt < SIMULATION_TIME; tc++ ){
    #if GNUPLOT_SHOW
    if(tc%GNUPLOT_OUT_INTERVAL==0)Output_gnuplot(gp,X,Particle_size,0.5,Dt*tc,f0,u0);
    #endif

      // Output x,y,vx,vy,theta to file
    if(tc%GNUPLOT_OUT_INTERVAL==0){
      for(int i=0; i<N; i++){
        cout << "#time=" << tc*Dt << " " << i << " " << x[i] << " " << y[i] << " " << vx[i] << " " << vy[i]
             << " " << theta[i] << " " << kfx[i] << " " << kfy[i] <<endl;
      }
    }

    // update by BAOAB SVV
    SVV_BAOAV(X, kF);
      
    check_pairlist(vx,vy,x,y);


  }
  return 0;
}


