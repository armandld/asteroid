#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <valarray>
#include <cmath> 
#include <numeric>
using namespace std;

class Exercice3
{

private:
  double t, dt, tFin;
  double a;
  double GM=6.674e-11;
  double msol, mjup, msat;
  double mtot;
  double Om;
  double xs,xj;
  int nsteps;
  int last;
  bool adapt;
  double alpha = 0e0;
  double beta = 0e0;
  double tol= 30;
  double f=0.99;
  int count_steps=1;
  valarray<double> y0 = std::valarray<double>(0.e0, 4); // Correctly initialized
  valarray<double> y  = std::valarray<double>(0.e0, 4); // Correctly initialized
  ofstream *outputFile;
  
  int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics

  void printOut(bool write)
  {
    if((!write && last>=sampling) || (write && last!=1))
    {
      double Energy = compute_energy(y[0],y[1],y[2],y[3]);
      *outputFile << count_steps << " " << t << " " << y[0] << " " << y[1] << " "<< y[2] << " " << y[3] << " " \
      << Energy<< " "<< endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }
  
  double dist_s(double x) const {return sqrt(pow(y[0]-x,2)+pow(y[1],2));}
  
  std::valarray<double> get_f(const std::valarray<double>& x, double t) {
    std::valarray<double> xbis(0.0, 4);
    //TO DO
    if (mjup != 0) {
			
		double dist_s_s = dist_s(xs);
		double dist_s_j = dist_s(xj);
		
		xbis[2] = - pow(Om, 2) * (pow(a,3)*beta*(x[0]+alpha*a)/pow(dist_s_s,3)
								 + pow(a,3)*alpha*(x[0]-beta*a)/pow(dist_s_j,3)
								 - x[0])
				   + 2 * Om * x[3];

		xbis[3] = - pow(Om, 2) * (pow(a,3)*beta*x[1]/pow(dist_s_s,3)
								 + pow(a,3)*alpha*x[1]/pow(dist_s_j,3)
								 - x[1])
				   - 2 * Om * x[2];
	} else {
		xbis[2] = - GM * msol  * x[0] / pow(dist_s(xs),3);
		xbis[3] = - GM * msol * x[1] / pow(dist_s(xs), 3);
	}
//AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    // Velocities
    xbis[0] = x[2];
    xbis[1] = x[3];

    return xbis;
  }  

// Function to compute mechanical energy per mass in R'
double compute_energy(double xx, double yy, double vx, double vy) {
	
	double dist_s_s = dist_s(xs);
    double dist_s_j = (mjup != 0) ? dist_s(xj) : 0;
	double energy = 0.5*msat*(vx*vx + vy*vy) // Energie cinétique
					- GM * msat * msol / dist_s_s // Énergie potentielle gravitationnelle Soleil;
					- 0.5 * msat * Om * Om * (xx*xx + yy*yy); // Énergie potentielle centrifuge
	
    if (mjup!=0) {energy += - GM * msat * mjup / dist_s_j;} // Énergie potentielle gravitationnelle Jupiter}
    return energy;
}

std::valarray<double> RK4_do_onestep(const std::valarray<double>& yold, double ti, double dt) {
    std::valarray<double> k1, k2, k3, k4, ynew;
	k1=dt*get_f(yold,ti);
	k2=dt*get_f(yold+0.5*k1,ti+0.5*dt);
	k3=dt*get_f(yold+0.5*k2,ti+0.5*dt);
	k4=dt*get_f(yold+k3,ti+dt);
    ynew=yold+(k1+2.0*k2+2.0*k3+k4)/6.0;
    return ynew;
}


public:
  Exercice3(int argc, char* argv[])
  {
    string inputPath("configuration.in.example"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tFin         = configFile.get<double>("tFin");            // t final (overwritten if N_excit >0)
    msol         = configFile.get<double>("msol");              // mass of the sun
    mjup         = configFile.get<double>("mjup");              // mass of the earth
    msat         = configFile.get<double>("msat");              // mass of the satellite
    a            = configFile.get<double>("a");               // distance sol jup
    adapt        = configFile.get<bool>("adapt");             //if 1=true -> adaptive dt, if 0=false -> fixed dt
    tol          = configFile.get<double>("tol");             //tolerance of the adaptive scheme
    nsteps       = configFile.get<int>("nsteps");        // number of time step per period
    sampling = configFile.get<unsigned int>("sampling",sampling);
    mtot=msol+mjup;
    alpha = mjup / mtot;
    beta = msol / mtot;
    //TO DO
    
    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
  }

  ~Exercice3()
  {
    outputFile->close();
    delete outputFile;
  };

  void run()
  {
    t = 0.e0;

    y[0]=2.0*a;
    y[1]=0.0;
    y[2]=-11000; 
    y[3]=2000;
    
    if (mjup!=0){
		Om = sqrt(GM *  (msol+mjup) / pow(a,3));
		xs = a*alpha;
		xj = -a*beta;
	}else{
		Om=0;
		xs=0;
	}
    
    last = 0;
    printOut(true);
    
    std::valarray<double> y1, y2;
    dt=tFin/nsteps;
    if (!adapt){ //Pas d'adaptation de dt
      for (int i = 0; i < nsteps; i++) {
            y = RK4_do_onestep(y, t, dt);
            t += dt;
			++count_steps;
            printOut(false);
        }
    }else{ 				// adaptation de dt
      while (t<=tFin) {
		  dt=min(dt,tFin-t);
		  ++count_steps;
		  y1=RK4_do_onestep(y,t,dt);
		  
		  std::valarray<double> y_inter=RK4_do_onestep(y,t,0.5*dt);
		  y2=RK4_do_onestep(y_inter,t+0.5*dt,0.5*dt);
		  
		  double d=sqrt(pow(y1[0] - y2[0], 2) + pow(y1[1] - y2[1], 2));
		  if (d<=tol){
			  y=y2;
			  t+=dt;
			  dt*=pow(tol/d,0.2); // if n=4, 1/(n+1) = 0.2
			  printOut(false);
			  }
		  else {
			  while (d>tol) {
				  dt*=f*pow(d/tol,0.2);			  
				  std::valarray<double>y_a=RK4_do_onestep(y,t,dt);
				  std::valarray<double>y_temp=RK4_do_onestep(y,t,0.5*dt);
			      std::valarray<double>y_b=RK4_do_onestep(y_temp,t+0.5*dt,0.5*dt);
			      d = sqrt(pow(y_a[0] - y_b[0], 2) + pow(y_a[1] - y_b[1], 2));			  
			  }
			  y=y2;
			  t+=dt;
			  printOut(false);
		  }
	  }
    };
    
    printOut(true); // ecrire le dernier pas de temps
  }

};

int main(int argc, char* argv[])
{
  Exercice3 engine(argc, argv);
  engine.run();

  return 0;
}
