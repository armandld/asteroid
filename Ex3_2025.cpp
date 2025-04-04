#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h"
#include <valarray>
#include <cmath>
#include <numeric>

using namespace std;

class Exercice3 {
private:
    double t, dt, tFin;
    double r0, r1, v0, a;
    double GM = 6.674e-11;
    valarray<double> m = std::valarray<double>(0.e0, 2);
    int N_excit, nsteps;
    int sampling;
    int last;
    int nsel_physics;
    bool adapt;
    double tol = 0e0;
    valarray<double> x0 = std::valarray<double>(0.e0, 4);
    valarray<double> x  = std::valarray<double>(0.e0, 4);
    ofstream *outputFile;

    void printOut(bool write) {
        if ((!write && last >= sampling) || (write && last != 1)) {
            double Energy = compute_energy(x[0], x[1], x[2], x[3]);
                *outputFile << t << " " << x[0] << " " << x[1] << " "
                        << x[2] << " " << x[3] << " " << Energy;
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
	
	if (mjup==0){
		Om=0;
		xs=0;
		y[0]=2.0*a;
		y[1]=0.0;
		y[2]=-11000; 
		y[3]=2000;
	}else{
		Om = sqrt(GM *  (msol+mjup) / pow(a,3));
		xs = -a*alpha;
		xj = a*beta;
		y[0]=2.0*a+xs;
		y[1]=0.0;
		y[2]=-11000 - Om * y[1]; //??????????????
		y[3]=2000 + Om * y[0];	//?????????????
		
	
	}    
            if (adapt) {
                *outputFile << " " << dt;  // ajouter dt en dernière colonne
            }
    
            *outputFile << endl;
    
            last = 1;
        } else {
            last++;
        }
    }


    std::valarray<double> get_f(double t, const std::valarray<double>& x) {
        std::valarray<double> xdot(0.0, 4);

        if (nsel_physics == 1) {
            double r = sqrt(x[0] * x[0] + x[1] * x[1]);
            double acc = -GM * m[0] / (r * r * r);
            xdot[0] = x[2];
            xdot[1] = x[3];
            xdot[2] = acc * x[0];
            xdot[3] = acc * x[1];
        } else if (nsel_physics == 2) {
            double omega = sqrt(GM * (m[0] + m[1]) / (a * a * a));
            double x_sun = -a * m[1] / (m[0] + m[1]);
            double x_jup =  a * m[0] / (m[0] + m[1]);

            double dx_s = x[0] - x_sun;
            double dy_s = x[1];
            double r_s = sqrt(dx_s * dx_s + dy_s * dy_s);

            double dx_j = x[0] - x_jup;
            double dy_j = x[1];
            double r_j = sqrt(dx_j * dx_j + dy_j * dy_j);

			double acc_sun_x = -GM * m[0] * dx_s / pow(r_s, 3);
            double acc_sun_y = -GM * m[0] * dy_s / pow(r_s, 3);

            double acc_jup_x = -GM * m[1] * dx_j / pow(r_j, 3);
            double acc_jup_y = -GM * m[1] * dy_j / pow(r_j, 3);

            xdot[0] = x[2];
            xdot[1] = x[3];
            xdot[2] = acc_sun_x + acc_jup_x + 2 * omega * x[3] + omega * omega * x[0];
            xdot[3] = acc_sun_y + acc_jup_y - 2 * omega * x[2] + omega * omega * x[1];
        } else {
            cerr << "No dynamics corresponds to this index" << endl;
        }

        return xdot;
    }

		double get_Epot(double xx, double yy) {
			if (nsel_physics == 1) {
				double r = sqrt(xx * xx + yy * yy);
				return -GM * m[0] / r;
			}
			else if (nsel_physics == 2) {
				// Cas avec Jupiter : référentiel du centre de masse
				double mS = m[0];
				double mJ = m[1];
				double mT = mS + mJ;

				double x_sun = -a * mJ / mT;
				double x_jup =  a * mS / mT;

				double dx_s = xx - x_sun;
				double dy_s = yy;

				double dx_j = xx - x_jup;
				double dy_j = yy;

				double r_s = sqrt(dx_s * dx_s + dy_s * dy_s);
				double r_j = sqrt(dx_j * dx_j + dy_j * dy_j);

				double Ep_sun = -GM * mS / r_s;
				double Ep_jup = -GM * mJ / r_j;

				return Ep_sun + Ep_jup;
			}
			else {
				cerr << "Erreur dans get_Epot : nsel_physics invalide" << endl;
				return 0.0;
			}
		}


double compute_energy(double xx, double yy, double vx, double vy) {
    double Ekin = 0.5 * (vx * vx + vy * vy);
    double Epot = get_Epot(xx, yy);

    if (nsel_physics == 2) {
        // Potentiel centrifuge dans R'
        double omega = sqrt(GM * (m[0] + m[1]) / pow(a, 3));
        double r2 = xx * xx + yy * yy;
        double Ecentrifuge = -0.5 * omega * omega * r2;
        return Ekin + Epot + Ecentrifuge;
    }

    return Ekin + Epot;
}


    void initial_condition() {
        double x_sun = -a * (m[1] / (m[0] + m[1]));
        double x_jupiter = a * (m[0] / (m[0] + m[1]));

        x0[0] = 2 * a;
        x0[1] = 0;
        x0[2] = -11000;
        x0[3] = 2000;

        if (nsel_physics == 2) {
            double omega = sqrt(GM * (m[0] + m[1]) / (a * a * a));
            double x_rel = x0[0];
            double y_rel = x0[1];
            double vx_rel = x0[2];
            double vy_rel = x0[3];
            x0[2] = vx_rel + omega * y_rel;
            x0[3] = vy_rel - omega * x_rel;
        }
    }

    std::valarray<double> RK4_do_onestep(const std::valarray<double>& yold, double ti, double dt) {
        std::valarray<double> k1 = dt * get_f(ti, yold);
        std::valarray<double> k2 = dt * get_f(ti + dt / 2.0, yold + 0.5 * k1);
        std::valarray<double> k3 = dt * get_f(ti + dt / 2.0, yold + 0.5 * k2);
        std::valarray<double> k4 = dt * get_f(ti + dt, yold + k3);
        return yold + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    }

    double compute_adaptive_dt(std::valarray<double>& y1, std::valarray<double>& y2, double dt) {
        double error = sqrt(pow(y1[0] - y2[0], 2) + pow(y1[1] - y2[1], 2));
        double factor;
        double f = 0.999;

        if (error <= tol + 1e-12) {
            factor = pow(tol / error, 1.0 / 5.0);
            dt = f * dt * factor;
        } else {
            while (error > tol) {
                factor = pow(tol / error, 1.0 / 5.0);
                dt = f * dt * factor;

                y1 = RK4_do_onestep(x, t, dt);
                std::valarray<double> y_intermediate = RK4_do_onestep(x, t, dt / 2.0);
                y2 = RK4_do_onestep(y_intermediate, t + dt / 2.0, dt / 2.0);
                error = sqrt(pow(y1[0] - y2[0], 2) + pow(y1[1] - y2[1], 2));
            }
        }

        if (t + dt > tFin) {//pour que tfin corresponde au temps final demandé
            dt = tFin - t;
        }

        return dt;
    }

public:
    Exercice3(int argc, char* argv[]) {
        string inputPath("configuration.in.example");
        if (argc > 1)
            inputPath = argv[1];

        ConfigFile configFile(inputPath);

        for (int i = 2; i < argc; ++i)
            configFile.process(argv[i]);

        tFin = configFile.get<double>("tFin");
        m[0] = configFile.get<double>("m1");
        m[1] = configFile.get<double>("m2");
        a = configFile.get<double>("a");
        nsel_physics = configFile.get<int>("nsel_physics");
        adapt = configFile.get<bool>("adapt");
        tol = configFile.get<double>("tol");
        sampling = configFile.get<int>("sampling");
        nsteps = configFile.get<int>("nsteps");

        outputFile = new ofstream(configFile.get<string>("output").c_str());
        outputFile->precision(15);

        dt = tFin / nsteps;
        //std::cout << "[DEBUG] adapt = " << adapt << ", nsteps = " << nsteps << ", tol = " << tol << std::endl;
    }

    ~Exercice3() {
        outputFile->close();
        delete outputFile;
    }

    void run() {
        t = 0.0;
        initial_condition();
        x = x0;
        last = 0;
        printOut(true);

        std::valarray<double> y1, y2;
        if (!adapt) {
            for (int i = 0; i < nsteps; i++) {
                x = RK4_do_onestep(x, t, dt);
                t += dt;
                printOut(false);
            }
        } else {
            while (t < tFin) {
                y1 = RK4_do_onestep(x, t, dt);
                std::valarray<double> y_intermediate = RK4_do_onestep(x, t, dt / 2.0);
                y2 = RK4_do_onestep(y_intermediate, t + dt / 2.0, dt / 2.0);
                dt = compute_adaptive_dt(y1, y2, dt);
                x = y2;
                t += dt;
                printOut(false);
            }
        }

        printOut(true);
    }
};

int main(int argc, char* argv[]) {
    
    Exercice3 engine(argc, argv);
    engine.run();
    return 0;
}

