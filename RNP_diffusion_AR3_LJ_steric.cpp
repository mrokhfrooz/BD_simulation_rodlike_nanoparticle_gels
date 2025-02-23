#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <unistd.h>  // Included this for the "getpid" function
#pragma warning(disable : 4996) // https://stackoverflow.com/questions/13550864/error-c4996-ctime-this-function-or-variable-may-be-unsafe
using namespace std;


double dotProduct(vector<double> a, vector<double> b) {  // input must be 2 3-elements vectors  
	double result_dot = 0;
	result_dot = (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
	return result_dot;
}

vector<double> crossProduct(vector<double> a, vector<double> b) { 
	vector<double> result(3);

	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];

	return result;
}


void parallelSepCon( vector <double> rjri, const vector<double> P_i, const vector<double> P_j, double dot_switzer, double L_switzer, double L_mucin, double& xmin, double& ymin, const vector<double> R_i, const vector<double> R_j) {
	
	
	double aaa= dotProduct(rjri, P_i); 
	double posneg;
	double epslion= 1e-14; // Small epsilon to account for round-off errors in floating-point comparisons


	if (abs(aaa) < 1e-12) {     
	
		xmin=0;
		ymin=0;
	}
	else if (abs (aaa)>= (L_mucin+L_switzer) ) {   
	    
		double sijp = dotProduct(rjri, P_i) + L_mucin * dot_switzer;
		double sijm = dotProduct(rjri, P_i) - L_mucin * dot_switzer;
		double sjip = -dotProduct(rjri, P_j) + L_switzer * dot_switzer;
		double sjim = -dotProduct(rjri, P_j) - L_switzer * dot_switzer;

		// Finding xmin for object i: RNP

		if (abs(sijp) < abs(sijm)) { 
			xmin = sijp;
		}
		else if (abs(sijp) > abs(sijm)) {
			xmin = sijm;
		}

		if (xmin >= L_switzer) { 
		xmin = L_switzer;
		}
		if (xmin <= -L_switzer) {
			xmin = -L_switzer;
		}

		// Finding ymin for Object j: Fiber
		if (abs(sjip) < abs(sjim)) {
			ymin = sjip;
		}
		else if (abs(sjip) > abs(sjim)) {
			ymin = sjim;
		}

		if (ymin >= L_mucin) {
			ymin = L_mucin;
		}
		if (ymin <= -L_mucin) {
			ymin = -L_mucin;
		}

	}
	else if (abs (aaa) <= epslion+L_mucin-L_switzer )  { 

		xmin=0;
		ymin= -dotProduct(rjri, P_j);

	}

	else if ( (abs (aaa)>=L_mucin-L_switzer) && (dot_switzer>0) && (abs (aaa) < L_mucin+L_switzer) ) {
	
		double sijp = dotProduct(rjri, P_i) + L_mucin * dot_switzer;
		double sijm = dotProduct(rjri, P_i) - L_mucin * dot_switzer;

		if (abs(sijp) < abs(sijm)) {
			xmin = sijp;
			posneg = 1;
		}
		else if (abs(sijp) > abs(sijm)) {
			xmin = sijm;
			posneg = -1;
		}

		xmin= (xmin-(posneg*L_switzer))/2;
		vector <double> midpoint = {R_i[0]+xmin*P_i[0], R_i[1]+xmin*P_i[1], R_i[2]+xmin*P_i[2]};
		vector <double> rjmidpoint = { R_j[0]-midpoint[0], R_j[1]-midpoint[1], R_j[2]-midpoint[2]};

		ymin= -dotProduct(rjmidpoint, P_j);
		
	}

	else if ( (abs (aaa)>=L_mucin-L_switzer) && (dot_switzer < 0) && (abs (aaa) < L_mucin+L_switzer) ) {
	
		double sijp = dotProduct(rjri, P_i) + L_mucin * dot_switzer;
		double sijm = dotProduct(rjri, P_i) - L_mucin * dot_switzer;


		if (abs(sijp) < abs(sijm)) {
			xmin = sijp;
			posneg = -1;
		}
		else if (abs(sijp) > abs(sijm)) {
			xmin = sijm;
			posneg = 1;
		}

		xmin= (xmin-(posneg*L_switzer))/2;
		vector <double> midpoint = {R_i[0]+xmin*P_i[0], R_i[1]+xmin*P_i[1], R_i[2]+xmin*P_i[2]};
		vector <double> rjmidpoint = { R_j[0]-midpoint[0], R_j[1]-midpoint[1], R_j[2]-midpoint[2]};

		ymin= -dotProduct(rjmidpoint, P_j);
	}

	else {
	xmin=0;  
	ymin=0;
	}

}







void RNP_fiber_distance (const vector<double> R_i, const vector<double> P_i, const vector<double> R_j, const vector<double> P_j, double L_switzer, double L_mucin, double d_mucin, double p, int uu , double& hij) {
	double d_effective= (p+d_mucin)/2;

	double dot_switzer = dotProduct(P_i, P_j);
	vector <double> rirj = { R_i[0] - R_j[0],R_i[1] - R_j[1], R_i[2] - R_j[2] };  // Ri-Rj
	vector <double> rjri = { -rirj[0], -rirj[1], -rirj[2] };   // Rj-Ri
	


	if (abs(dot_switzer * dot_switzer - 1) <= 1e-6) {  // Handling the case where RNP and fiber are in parallel

 
		double xmin; // xmin : distance from COM of RNP to a point along the RNP axis at which shortest distance between RNP and fiber occurs
		double ymin; // ymin : same concept as xmin, but is applicable to the fiber

		
		parallelSepCon(rjri, P_i, P_j, dot_switzer, L_switzer, L_mucin,xmin, ymin, R_i, R_j); // Calling parallelSepCon function

		double	sepahan = sqrt(pow(rjri[0] + ymin * P_j[0] - xmin * P_i[0], 2) + pow(rjri[1] + ymin * P_j[1] - xmin * P_i[1], 2) + pow(rjri[2] + ymin * P_j[2] - xmin * P_i[2], 2));



		double gij[] = { rjri[0] + ymin * P_j[0] - xmin * P_i[0], rjri[1] + ymin * P_j[1] - xmin * P_i[1], rjri[2] + ymin * P_j[2] - xmin * P_i[2] };

		double gij_norm = sepahan;
		double hij = gij_norm - d_effective;  // Surface-to-surface distance in obtained by having the center-to-center distance and effective diameter of the objects



	}

	
	else {  // Finding distance between RNP and fiber when they are not in parallel 


		double term1 = dotProduct(rirj, P_j); 
		double term2 = dotProduct(rirj, P_i); 
		double term3 = -1 * term2; 
		double term4 = -1 * term1; 
		vector <double> sij(9);  
		vector <double> sji(9);
		int it = 0;

		sij[0] = (term1 * dot_switzer - term2) / (1 - pow(dot_switzer, 2));  
		sji[0] = (term3 * dot_switzer - term4) / (1 - pow(dot_switzer, 2));


		double sep = sqrt(pow(rjri[0] + sji[0] * P_j[0] - sij[0] * P_i[0], 2) + pow(rjri[1] + sji[0] * P_j[1] - sij[0] * P_i[1], 2) + pow(rjri[2] + sji[0] * P_j[2] - sij[0] * P_i[2], 2));
		

		if (abs(sij[0]) < L_switzer && abs(sji[0]) < L_mucin) {   // Scenario #1

			it = 0;
		}
		else {
			vector <double> termA = { rjri[0] + L_mucin * P_j[0], rjri[1] + L_mucin * P_j[1] , rjri[2] + L_mucin * P_j[2] };
			sij[1] = dotProduct(termA, P_i);
			sji[1] = L_mucin;
			// Scenario #2 



			vector <double> termB = { rjri[0] - L_mucin * P_j[0],rjri[1] - L_mucin * P_j[1] , rjri[2] - L_mucin * P_j[2] };
			sij[2] = dotProduct(termB, P_i);
			sji[2] = -L_mucin;
			// Scenario #3


			sij[3] = L_switzer;
			vector <double> termC = { rirj[0] + L_switzer * P_i[0],rirj[1] + L_switzer * P_i[1],rirj[2] + L_switzer * P_i[2] };
			sji[3] = dotProduct(termC, P_j);
			// Scenario #4


			sij[4] = -L_switzer;
			vector <double> termD = { rirj[0] - L_switzer * P_i[0],rirj[1] - L_switzer * P_i[1],rirj[2] - L_switzer * P_i[2] };
			sji[4] = dotProduct(termD, P_j);
			// Scenario #5


			sij[5] = L_switzer;
			sji[5] = L_mucin;
			// Scenario #6


			sij[6] = L_switzer;
			sji[6] = -L_mucin;
			// Scenario #7



			sij[7] = -L_switzer;
			sji[7] = L_mucin;
			// Scenario #8


			sij[8] = -L_switzer;
			sji[8] = -L_mucin;
			// Scenario #9

			sep = 100 * 1e-9; // A large distance is initially considered as the distance between the objects, if their distance is < sep, then it is updated

			for (int ii = 1; ii < 9; ii++) {  // ii starts from 1, because Scenario #1 was already considered in the "if" block
				double 	sep_tmp = sqrt(pow(rjri[0] + sji[ii] * P_j[0] - sij[ii] * P_i[0], 2) + pow(rjri[1] + sji[ii] * P_j[1] - sij[ii] * P_i[1], 2) + pow(rjri[2] + sji[ii] * P_j[2] - sij[ii] * P_i[2], 2));
				if (sep_tmp < sep && abs(sij[ii]) <= L_switzer && abs(sji[ii]) <= L_mucin) {
					sep = sep_tmp;
					it = ii;
				}

			}
		}


		double gij[] = { rjri[0] + sji[it] * P_j[0] - sij[it] * P_i[0], rjri[1] + sji[it] * P_j[1] - sij[it] * P_i[1], rjri[2] + sji[it] * P_j[2] - sij[it] * P_i[2] };
		double gij_norm = sep;
		
		 hij = gij_norm - d_effective;  // hij is the SURFACE_to_SURFACE distance, gij_norm is CENTER_TO_CENTER distance


	}



} 


void center_to_center(const vector<double> R_i, const vector<double> P_i, const vector<double> R_j, const vector<double> P_j, double L_switzer, double L_mucin,  double& gij_norm, vector<double>& Gij ) {

	double dot_switzer = dotProduct(P_i, P_j);
	vector <double> rirj = { R_i[0] - R_j[0],R_i[1] - R_j[1], R_i[2] - R_j[2] }; // Ri-Rj
	vector <double> rjri = { -rirj[0], -rirj[1], -rirj[2] };  // Rj-Ri
	

	if (abs(dot_switzer * dot_switzer - 1) <= 1e-6) { 


		double xmin;
		double ymin;


		parallelSepCon(rjri, P_i, P_j, dot_switzer, L_switzer, L_mucin,xmin, ymin, R_i, R_j);



		double	sepahan = sqrt(pow(rjri[0] + ymin * P_j[0] - xmin * P_i[0], 2) + pow(rjri[1] + ymin * P_j[1] - xmin * P_i[1], 2) + pow(rjri[2] + ymin * P_j[2] - xmin * P_i[2], 2));



		double gij[] = { rjri[0] + ymin * P_j[0] - xmin * P_i[0], rjri[1] + ymin * P_j[1] - xmin * P_i[1], rjri[2] + ymin * P_j[2] - xmin * P_i[2] };
		 gij_norm = sepahan;


	    Gij = { xmin * P_i[0] ,  xmin * P_i[1] ,  xmin * P_i[2] };



	}


	else {


		double term1 = dotProduct(rirj, P_j); 
		double term2 = dotProduct(rirj, P_i); 
		double term3 = -1 * term2; 
		double term4 = -1 * term1; 
		vector <double> sij(9);  
		vector <double> sji(9);
		int it = 0;

		sij[0] = (term1 * dot_switzer - term2) / (1 - pow(dot_switzer, 2));  
		sji[0] = (term3 * dot_switzer - term4) / (1 - pow(dot_switzer, 2));



		double sep = sqrt(pow(rjri[0] + sji[0] * P_j[0] - sij[0] * P_i[0], 2) + pow(rjri[1] + sji[0] * P_j[1] - sij[0] * P_i[1], 2) + pow(rjri[2] + sji[0] * P_j[2] - sij[0] * P_i[2], 2));
		

		if (abs(sij[0]) < L_switzer && abs(sji[0]) < L_mucin) {   // Scenario #1

			it = 0;
		}
		else {
			vector <double> termA = { rjri[0] + L_mucin * P_j[0], rjri[1] + L_mucin * P_j[1] , rjri[2] + L_mucin * P_j[2] };
			sij[1] = dotProduct(termA, P_i);
			sji[1] = L_mucin;
			



			vector <double> termB = { rjri[0] - L_mucin * P_j[0],rjri[1] - L_mucin * P_j[1] , rjri[2] - L_mucin * P_j[2] };
			sij[2] = dotProduct(termB, P_i);
			sji[2] = -L_mucin;
			


			sij[3] = L_switzer;
			vector <double> termC = { rirj[0] + L_switzer * P_i[0],rirj[1] + L_switzer * P_i[1],rirj[2] + L_switzer * P_i[2] };
			sji[3] = dotProduct(termC, P_j);
			


			sij[4] = -L_switzer;
			vector <double> termD = { rirj[0] - L_switzer * P_i[0],rirj[1] - L_switzer * P_i[1],rirj[2] - L_switzer * P_i[2] };
			sji[4] = dotProduct(termD, P_j);
			


			sij[5] = L_switzer;
			sji[5] = L_mucin;
			



			sij[6] = L_switzer;
			sji[6] = -L_mucin;
			



			sij[7] = -L_switzer;
			sji[7] = L_mucin;
			


			sij[8] = -L_switzer;
			sji[8] = -L_mucin;
		

			sep = 100 * 1e-9;

			for (int ii = 1; ii < 9; ii++) {  

				double 	sep_tmp = sqrt(pow(rjri[0] + sji[ii] * P_j[0] - sij[ii] * P_i[0], 2) + pow(rjri[1] + sji[ii] * P_j[1] - sij[ii] * P_i[1], 2) + pow(rjri[2] + sji[ii] * P_j[2] - sij[ii] * P_i[2], 2));
				if (sep_tmp < sep && abs(sij[ii]) <= L_switzer && abs(sji[ii]) <= L_mucin) {
					sep = sep_tmp;
					it = ii;
				}

			}
		}


		double gij[] = { rjri[0] + sji[it] * P_j[0] - sij[it] * P_i[0], rjri[1] + sji[it] * P_j[1] - sij[it] * P_i[1], rjri[2] + sji[it] * P_j[2] - sij[it] * P_i[2] };
		gij_norm = sep;

		Gij = { sij[it] * P_i[0] ,  sij[it] * P_i[1] ,  sij[it] * P_i[2]  };




	}



} 



void LJ_force_torque_calc (const vector<double> R_i, const vector<double> P_i, const vector<double> R_j, const vector<double> P_j, vector<double>& Force_temp, vector<double>& Tork_temp, double epsilon, double delta, double cutoff_LJ, double L_switzer, double L_mucin, double d_mucin, double p, int uu ,int j, int k) {


	double gij;  // center to center distance
	vector <double> arm (3); // Arm of Torque
	center_to_center(R_i, P_i, R_j, P_j,  L_switzer,L_mucin ,gij, arm);

	if (gij > cutoff_LJ)	{
			Force_temp = { 0,0,0 };
			Tork_temp = { 0,0,0 };
		}
		else {

			double s=p+d_mucin; // Steric diameter
			double  disx,disy,disz;
			double  dis= gij;
			vector <double> R_ix (3);
			vector <double> R_iy (3);
			vector <double> R_iz (3);
			R_ix = { R_i[0]+delta , R_i[1], 		 R_i[2] };
			R_iy = { R_i[0] , 		R_i[1]+delta, 	 R_i[2] };
			R_iz = { R_i[0] , 		R_i[1],			 R_i[2]+delta };

			double FSx, FSy, FSz;
			vector <double> arm_fake (3); // I do not use this arm in calculation of the torque, that's why called fake

			center_to_center(R_ix, P_i, R_j, P_j,  L_switzer, L_mucin , disx, arm_fake);
			center_to_center(R_iy, P_i, R_j, P_j,  L_switzer, L_mucin , disy, arm_fake);
			center_to_center(R_iz, P_i, R_j, P_j,  L_switzer, L_mucin , disz, arm_fake);
			
			// Used this video to simplify LJ calculation : https://www.youtube.com/watch?v=CIbKCTflpG8
			double s_r= s/(2*dis);
			double A6= pow (s_r,6); 
			double U_basis= A6*(A6-1);
		  	FSx=(epsilon/delta)*((pow ((s/(2*disx)),12)-pow ((s/(2*disx)),6))-U_basis);
  	  		FSy=(epsilon/delta)*((pow ((s/(2*disy)),12)-pow ((s/(2*disy)),6))-U_basis);
  	  		FSz=(epsilon/delta)*((pow ((s/(2*disz)),12)-pow ((s/(2*disz)),6))-U_basis);

			Force_temp = { FSx,FSy,FSz };


			Tork_temp = crossProduct(arm, Force_temp);

		}


		if (Force_temp[0] > 1e-8 || Force_temp[1] > 1e-8 || Force_temp[2] > 1e-8) {
		
			// In some cases, overlap between RNP and fiber might occur (large time step, very concentrated system, use this block of code for debugging and printing the parameters)	
			cout << "Large forces are generated" << endl;

		}


} 




int main() {


	const int N = 10; // 10; //number of rod-like NPs
	const int step = 2000000; //2000000 ; // Number of steps
	const int Saveeverystep = 1000; //  the output of the code is saved @ every "Saveeverystep" step. Larger storage is needed if you use smaller "Saveeverystep"
	

	int SaveN = 0;
	for (int kk = 0; kk < step + 1; kk++) {
		if (fmod(kk, Saveeverystep) == 0) {
			SaveN = SaveN + 1;
		}
	}
	

	long pid = static_cast<long>(getpid()); // I used computer's process ID for seeding random number generator
	srand(pid);  // Seed the random number generator 
	
	float pidpid=pid;
	
	
	auto start = chrono::high_resolution_clock::now(); // chat gpb this part
	// get the current time
	time_t currentTime = time(0);
	// convert the current time to a string
	char* timeString = ctime(&currentTime);
	// print the current time
	cout << "Simulation started at: " << timeString << endl;


	const double pi = 3.14159265358979323846;
	double k_B = 1.3806e-23; // unit J / K
	double T = 300;  // temperature
	double mu = 0.001; // Pa·s% viscosity
	

	double p =  85.958e-9; // RNP diameter
	double radius= p/2.0; // RNP radius


	double L_lowen = 257.875151e-9; // Total length of RNP with end-caps
    double L_chen = L_lowen-p; // Length of RNP without the end-caps
	double L_switzer = L_chen / 2; // Switzer (Switzer III, Leonard H. Simulating systems of flexible fibers. The University of Wisconsin-Madison, 2002.) wrote the formulation based on L_Chen/2
	

	double b = 1000e-9; // Periodic box size

	double d_mucin = 0 * 1e-9;// mucin diameter[meter]
	double L_mucin_total= 1000e-9; //
	double L_mucin= (L_mucin_total)/2; // Half length of fiber



	double pratio = L_lowen / p; // Aspect ratio, See Fig 1's caption in : Löwen, Hartmut. "Brownian dynamics of hard spherocylinders." Physical Review E 50.2 (1994): 1232.
	double r_p= pratio-1; // Chen et al. defined the aspect ratio without considering the end-caps: Chen, Jing-Yao, et al. "Rheology and structure of suspensions of spherocylinders via Brownian dynamics simulations." Journal of Rheology 65.2 (2021): 273-288.
	

	double D_not = (k_B * T) / (L_lowen * mu); // equation 4 in lowen's paper
	
	// equations 1 - 3 in Lowen's paper, three short time diffusivitiy
	double D_amodi = (D_not / (4 * pi)) * (log(pratio) + 0.839 + (0.185 / pratio) + (0.233 / pow(pratio, 2))); //% equation 1
	double D_movazi = (D_not / (2 * pi)) * (log(pratio) - 0.207 + (0.98 / pratio) - (0.133 / pow(pratio, 2))); //% equation 2
	double D_r = (3 * D_not / (pi * L_lowen * L_lowen)) * (log(pratio) - 0.662 + (0.917 / pratio) - (0.05 / pow(pratio, 2))); //% equation 3
	double D0 = (1.0 / 3.0) * (2 * D_amodi + D_movazi);



	double zamanscale_lowen = p * p / D_not;
	double dt_lowen = 0.0005 * zamanscale_lowen;  
//	double endzaman_lowen = 35 * zamanscale_lowen;  


	double dt = 1.0 * dt_lowen; // time step of simulation
	double zamanscale1 = zamanscale_lowen;
	double endzaman = step * dt;  




	// LJ potentialll
	double epsilon= -4*(1*k_B*T); //See Eq. (7) in: Hansing, Johann, et al. "Nanoparticle filtering in charged hydrogels: Effects of particle size, charge asymmetry and salt concentration." The European Physical Journal E 39 (2016): 1-13.
	double cutoff_LJ=(pow(2,(-5.0/6.0))*(p+d_mucin)); // LJ cutoff



	vector <double>  zaman(step + 1); // zaman means time, (used this word as some of C++ libraries already used the keyword of "time")
	vector <double>  zamanD(step + 1); // D: Dimensionless

	vector <double>  zaman_internal(SaveN);  // Number of data points for internal sampling can be different
	vector <double>  zamanD_internal(SaveN); // Capital D: dimensionlesss time

	zaman[0] = 0;

	for (int i = 1; i < step + 1; i++) {
		zaman[i] = zaman[i - 1] + dt;
	}

	for (int i = 1; i < SaveN; i++) {
		zaman_internal[i] = zaman_internal[i - 1] + Saveeverystep * dt;
		zamanD_internal[i] = zaman_internal[i] / zamanscale1;
	}


	for (int i = 0; i < step + 1; i++) {
		zamanD[i] = zaman[i] / zamanscale1;
	}

	vector <int>     counter(step + 1);
	vector <double>  msd(step + 1);
    vector <double>  Ati(step + 1); // Equ 10 in "Molecular simulation of the diffusion mechanism of nanorods in cross-linked networks"†
    vector <double>  msd_movazi(step + 1);
	vector <double>  msd_internal(step + 1);
	vector <double>  msd_theory(step + 1);
	vector <double>  mdr(step + 1);
	vector <double>  mdr_internal(step + 1);
	vector <double>  slope(step + 1);
	vector <double>  slope_internal(step + 1);
	vector <double>  D_rotation(step + 1);
	vector <double>  D_rotation_internal(step + 1);
	vector <double>  Drratio(step + 1);
	vector <double>  Drratio_internal(step + 1);
	vector <double>  Dtransratio(step + 1);
	vector <double>  Dtransratio_internal(step + 1);
	vector <double>  error(step + 1);
	
	for (int i = 0; i <= step; i++) {
		counter[i] = 0;
		msd_theory[i] = 0;
		msd[i] = 0;
		msd_movazi[i]=0;
		Ati[i]=0;
		msd_internal[i] = 0;
		mdr[i] = 0;
		mdr_internal[i] = 0;
		error[i] = 0;
		slope[i] = 0;
		slope_internal[i] = 0;
		D_rotation[i] = 0;
		D_rotation_internal[i] = 0;
		Drratio[i] = 0;
		Drratio_internal[i] = 0;
		Dtransratio[i] = 0;
		Dtransratio_internal[i] = 0;
	}
	


	double omega[3][N];
	double omega0[3][N];
	double e1[3][N];
	double e2[3][N];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < N; j++) {
			omega[i][j] = 0;
		}
	}
	
	for (int j = 0; j < N; j++) {
		omega[2][j] = 1;
	}
	
	double initial_cord = b/2; // Place all RNPs at the center of the box initially, later release them randomly in the free space of the gel, Othwewise, C++ might use garbage values

	static  double qx[2][N];
	static	double qx_total[2][N];
	static	double qy[2][N];
	static	double qy_total[2][N];
	static	double qz[2][N];
	static	double qz_total[2][N];
	static	double qx0[1][N];
	static	double qy0[1][N];
	static	double qz0[1][N];

	vector <vector<double >> save_qx(SaveN, vector<double>(N));
	vector <vector<double >> save_qx_total(SaveN, vector<double>(N));
	vector <vector<double >> save_qy(SaveN, vector<double>(N));
	vector <vector<double >> save_qy_total(SaveN, vector<double>(N));
	vector <vector<double >> save_qz(SaveN, vector<double>(N));
	vector <vector<double >> save_qz_total(SaveN, vector<double>(N));

	vector <vector<double >> save_omega(3 * SaveN, vector<double>(N));


	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < N; j++) {
			qx[i][j] = initial_cord;
			qy[i][j] = initial_cord;
			qz[i][j] = initial_cord;
		}
	}

	
	// Bulding the random gel

	const int Nmucin=160; // Number of fibers  
	double Nu= (Nmucin)/(b*b*b);  //   Number of fibers per unit volume

	double TwoNuL= Nu*(L_mucin_total) ; // Total fiber length per unit volume, Table 2 of Wang et al. https://doi.org/10.1016/j.jmps.2017.12.014

	double mucin_omega[3][Nmucin];
	double mucin_location[3][Nmucin];


	double aa = 0* b * 1e9;  // create random numbers inside the box for x,y,z coordinates
	double bb = 1* b * 1e9;  // create random numbers inside the box for x,y,z coordinates
	int a_int = static_cast<int>(aa);
	int b_int = static_cast<int>(bb);

	for (int jj=0;jj<Nmucin;jj++) {
		
		int random_num = rand() % (b_int - a_int + 1) + a_int;
		double x_candidate = random_num * 1e-9;
		random_num = rand() % (b_int - a_int + 1) + a_int;
		double y_candidate = random_num * 1e-9;
		random_num = rand() % (b_int - a_int + 1) + a_int;
		double z_candidate = random_num * 1e-9;

		int random_orient_0 = rand();
		double pi0 = (double)random_orient_0 / RAND_MAX * 2 - 1;
		int random_orient_1 = rand();
		double pi1 = (double)random_orient_1 / RAND_MAX * 2 - 1;
		int random_orient_2 = rand();
		double pi2 = (double)random_orient_2 / RAND_MAX * 2 - 1;
		double norm_pi_vec = sqrt(pi0 * pi0 + pi1 * pi1 + pi2 * pi2);  // To normalize unit vectors 
		pi0 = pi0 / norm_pi_vec;
		pi1 = pi1 / norm_pi_vec;
		pi2 = pi2 / norm_pi_vec;

		mucin_location[0][jj] = x_candidate;
		mucin_location[1][jj] = y_candidate;
		mucin_location[2][jj] = z_candidate;
		mucin_omega[0][jj] = pi0;
		mucin_omega[1][jj] = pi1;
		mucin_omega[2][jj] = pi2;
	}

	ofstream myfile88("save_boxomega.txt");  // save Fiber orientaions!
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < Nmucin; j++) {
			myfile88 << mucin_omega[k][j] << ";";
		}
		myfile88 << endl;
	}
	myfile88.close();

	ofstream kyfile90("save_boxLocation.txt"); // save fiber's COM in a single txt file
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < Nmucin; j++) {

			kyfile90 << mucin_location[k][j] << ";";

		}
		kyfile90 << endl;
	}
	kyfile90.close();
	





// Finding random position and orientation for RNPs in the space of gel with no contact with fibers. "Initialization"

	int shomarande = 1; // counter
	aa = 0.00* b * 1e9;  // create random numbers inside the box for x,y,z coordinates
	bb = 1.00* b * 1e9;  // create random numbers inside the box for x,y,z coordinates
	a_int = static_cast<int>(aa);
	b_int = static_cast<int>(bb);

	cout << "Initialization started " << endl;

	vector <double> R_j_periodic= {0,0,0};


	
	while (shomarande < N) {  
	
		int random_num = rand() % (b_int - a_int + 1) + a_int;
		double x_candidate = random_num * 1e-9;
		random_num = rand() % (b_int - a_int + 1) + a_int;
		double y_candidate = random_num * 1e-9;
		random_num = rand() % (b_int - a_int + 1) + a_int;
		double z_candidate = random_num * 1e-9;

		int random_orient_0 = rand();
		double pi0 = (double)random_orient_0 / RAND_MAX * 2 - 1;
		int random_orient_1 = rand();
		double pi1 = (double)random_orient_1 / RAND_MAX * 2 - 1;
		int random_orient_2 = rand();
		double pi2 = (double)random_orient_2 / RAND_MAX * 2 - 1;
		double norm_pi_vec = sqrt(pi0 * pi0 + pi1 * pi1 + pi2 * pi2);  // to normalize the vector
		pi0 = pi0 / norm_pi_vec;
		pi1 = pi1 / norm_pi_vec;
		pi2 = pi2 / norm_pi_vec;

		vector <double> R_K = { x_candidate ,y_candidate,  z_candidate }; // New particle placed randomly
		vector <double> P_K = { pi0 ,pi1,  pi2 };

		
		int switzer_counter = 0; // Counter 

		
		for (int uu = 0; uu < Nmucin ; uu++) 
		{ // Checking if RNP would contact with any fiber in the randomly selected location
		
				double hij = 0.0;
				vector <double> R_j = { mucin_location[0][uu] , mucin_location[1][uu], mucin_location[2][uu] }; 
				vector <double> P_j = { mucin_omega[0][uu] ,mucin_omega[1][uu],mucin_omega[2][uu] }; 

				// We need to check not only with the fiber in the main cell, but with all its periodic images
				for (int w1=-1;w1<2;w1++) {
					for (int w2=-1;w2<2;w2++) {
						for (int w3=-1;w3<2;w3++) {

							R_j_periodic[0] = R_j[0] + (b * w1);     
							R_j_periodic[1] = R_j[1] + (b * w2);
							R_j_periodic[2] = R_j[2] + (b * w3);							
							
							RNP_fiber_distance (R_K, P_K, R_j_periodic, P_j, L_switzer,L_mucin,d_mucin , p, uu, hij);
							
							if (hij > 10e-9) {
								switzer_counter = switzer_counter + 1;
							}
							else {
								break;
							}

						}
					}
				}			
        }

		if (switzer_counter == (Nmucin*27)) { // When this condition met, RNP is successfully initialized
			qx[0][shomarande] = x_candidate;
			qy[0][shomarande] = y_candidate;
			qz[0][shomarande] = z_candidate;
			omega[0][shomarande] = pi0;
			omega[1][shomarande] = pi1;
			omega[2][shomarande] = pi2;
			shomarande = shomarande + 1;
		}
	}


	qx[0][0] = qx[0][1];
	qy[0][0] = qy[0][1];
	qz[0][0] = qz[0][1];
	omega[0][0] =omega[0][1];
	omega[1][0] = omega[1][1];
	omega[2][0] = omega[2][1];
	// RNP 1 &2 have the same initial location and orientation




	// Printing the current time in the console
	// get the current time
	currentTime = time(0);
	// convert the current time to a string
	timeString = ctime(&currentTime);
	// print the current time
	cout << "Initialization succusfully done at: " << timeString << endl;




	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < N; j++) {
			omega0[i][j] = omega[i][j];
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < N; j++) {
			save_omega[i][j] = omega[i][j];
		}
	}


	// e1 , e2 calculation: they are perpendicular to RNP's orientation vector
	
	double W[] = { 1 / sqrt(3),1 / sqrt(3), 1 / sqrt(3) };

	for (int j = 0; j < N; j++) {
		e1[0][j] = (omega[1][j] * W[2]) - (omega[2][j] * W[1]);
		e1[1][j] = (omega[2][j] * W[0]) - (omega[0][j] * W[2]);
		e1[2][j] = (omega[0][j] * W[1]) - (omega[1][j] * W[0]);


		if ((e1[0][j] == 0) && (e1[1][j] == 0) && (e1[2][j] == 0))
		{
			W[0] = -1 / sqrt(3);
			W[1] = -1 / sqrt(3);
			W[2] = -1 / sqrt(3);
			e1[0][j] = (omega[1][j] * W[2]) - (omega[2][j] * W[1]);
			e1[1][j] = (omega[2][j] * W[0]) - (omega[0][j] * W[2]);
			e1[2][j] = (omega[0][j] * W[1]) - (omega[1][j] * W[0]);
		}
		
		// e1 Must be normalzied
		double norm_e1 = sqrt((e1[0][j] * e1[0][j]) + (e1[1][j] * e1[1][j]) + (e1[2][j] * e1[2][j]));
		e1[0][j] = e1[0][j] / norm_e1;
		e1[1][j] = e1[1][j] / norm_e1;
		e1[2][j] = e1[2][j] / norm_e1;



		e2[0][j] = (omega[1][j] * e1[2][j]) - (omega[2][j] * e1[1][j]);
		e2[1][j] = (omega[2][j] * e1[0][j]) - (omega[0][j] * e1[2][j]);
		e2[2][j] = (omega[0][j] * e1[1][j]) - (omega[1][j] * e1[0][j]);

		// e2 Must be normalzied
		double norm_e2 = sqrt((e2[0][j] * e2[0][j]) + (e2[1][j] * e2[1][j]) + (e2[2][j] * e2[2][j]));
		e2[0][j] = e2[0][j] / norm_e2;
		e2[1][j] = e2[1][j] / norm_e2;
		e2[2][j] = e2[2][j] / norm_e2;
	}


	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < N; j++) {
			qx_total[i][j] = qx[i][j];
			qy_total[i][j] = qy[i][j];
			qz_total[i][j] = qz[i][j];
		}
	}
	
	for (int i = 0; i < 1; i++) {
		for (int j = 0; j < N; j++) {
			qx0[i][j] = qx[i][j];
			qy0[i][j] = qy[i][j];
			qz0[i][j] = qz[i][j];
			save_qx[i][j] = qx[i][j];
			save_qy[i][j] = qy[i][j];
			save_qz[i][j] = qz[i][j];
			save_qx_total[i][j] = qx[i][j];
			save_qy_total[i][j] = qy[i][j];
			save_qz_total[i][j] = qz[i][j];
		}
	}

	// Qx, Qy, Qz amodi, movazi, movazi= parallel, qx_amodi = perpendicular
	static	double qx_movazi[2][N];
	static	double qx_amodi[2][N];
	static	double qy_movazi[2][N];
	static	double qy_amodi[2][N];
	static	double qz_movazi[2][N];
	static	double qz_amodi[2][N];


	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < N; j++) {
			qx_movazi[i][j] = 0;
			qx_amodi[i][j] = 0;
			qy_movazi[i][j] = 0;
			qy_amodi[i][j] = 0;
			qz_movazi[i][j] = 0;
			qz_amodi[i][j] = 0;
		}
	}
	
	
	for (int j = 0; j < N; j++) {  
		double dot_product = (qx0[0][j] * omega0[0][j]) + (qy0[0][j] * omega0[1][j]) + (qz0[0][j] * omega0[2][j]);
		qx_movazi[0][j] = dot_product * omega0[0][j];
		qx_amodi[0][j] = qx[0][j] - qx_movazi[0][j];
		qy_movazi[0][j] = dot_product * omega0[1][j];
		qy_amodi[0][j] = qy[0][j] - qy_movazi[0][j];
		qz_movazi[0][j] = dot_product * omega0[2][j];
		qz_amodi[0][j] = qz[0][j] - qz_movazi[0][j];
	}

	double progress[20]; 
	progress[0] = round(0.05 * step);
	for (int i = 1; i < 20; i++) {
		progress[i] = round(progress[i - 1] + (0.05 * step));
	}


   // Print all parameter of the program in a single txt file

	ofstream ourfile("parameters.txt");
	ourfile << "Process ID-rand seed number" << ";" << pidpid << endl;
	ourfile << "N" << ";" << N << endl;

	ourfile << "step" << ";" << step << endl;
	ourfile << "Saveeverystep" << ";" << Saveeverystep << endl;
	ourfile << "#ofSavedDataPoints" << ";" << SaveN << endl;
	ourfile << "b-mesh-size" << ";" << b << endl;
	ourfile << "b-mesh-size[nm]" << ";" << b*1e9 << endl;
	
	ourfile << "k_B" << ";" << k_B << endl;
	ourfile << "T" << ";" << T << endl;
	ourfile << "mu" << ";" << mu << endl;

	ourfile << "p-partcile-diameter" << ";" << p << endl;
	ourfile << "p-partcile-diameter-[nm]" << ";" << p*1e9 << endl;
	ourfile << "L_lowen" << ";" << L_lowen << endl;
	ourfile << "L_lowen[nm]" << ";" << L_lowen*1e9 << endl;
	ourfile << "L_chen" << ";" << L_chen << endl;
	ourfile << "L_chen[nm]" << ";" << L_chen*1e9 << endl;
	ourfile << "L_switzer" << ";" << L_switzer << endl;
	ourfile << "L_switzer[nm]" << ";" << L_switzer*1e9 << endl;
	ourfile << "Aspect-pratio_lowen" << ";" << pratio << endl;
	ourfile << "a-mucin-diameter" << ";" << d_mucin << endl;
 	ourfile << "a-mucin-diameter[nm]" << ";" << d_mucin*1e9 << endl;
 	ourfile << "Half L-mucin-diameter" << ";" << L_mucin_total/2 << endl;
 	ourfile << "Half L-mucin-diameter[nm]" << ";" << (L_mucin_total/2)*1e9 << endl;
	ourfile << "NumverOfFibers" << ";" << Nmucin << endl;
	ourfile << "2*nu*L_fiber" << ";" << TwoNuL << endl;
	ourfile << "Nu" << ";" << Nu << endl;

	ourfile << "**************" << ";" << endl;
	ourfile << "zamanscale1_lowen" << ";" << zamanscale_lowen << endl;
	ourfile << "dt_lowen" << ";" << dt_lowen << endl;

	ourfile << "**************" << ";" << endl;


	ourfile << "**************" << ";" << endl;
	ourfile << "dt_selected" << ";" << dt << endl;
	ourfile << "dt_selected/tav" << ";" << dt/zamanscale_lowen << endl;
	ourfile << "endzaman_Given" << ";" << endzaman << endl;
	ourfile << "**************" << ";" << endl;
	ourfile << "D_amodi" << ";" << D_amodi << endl;
	ourfile << "D_movazi" << ";" << D_movazi << endl;
	ourfile << "D_r" << ";" << D_r << endl;
	ourfile << "D0_translation_avg" << ";" << D0 << endl;
	ourfile << "D_not_lowen" << ";" << D_not << endl;
	ourfile << "**************" << ";" << endl;

	ourfile << "*** Steric interaction parameters***" << ";" << endl;
	ourfile << "Cutoff-LJ [m]" << ";" << cutoff_LJ << endl;
	ourfile << "Cutoff-LJ [nm]" << ";" << cutoff_LJ*1e9 << endl;
	ourfile << "epsilon [J]" << ";" << epsilon << endl;
	
	ourfile.close();

	
	
	// Print initial Omega of RNP before start of simulation
	ofstream yourfile1("RNP_omega_initial_simulation.txt");
	for (int k = 0; k < 3; k++) {
		for (int j = 0; j < N; j++) {
			yourfile1 << omega[k][j] << ";";
		}
		yourfile1 << endl;
	}
	yourfile1.close();





	ofstream yourfile20("RNP_coordinates_initial_simulation.txt");
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile20 << qx0[k][j] << ";";
		}
		yourfile20 << endl;
	}
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile20 << qy0[k][j] << ";";
		}
		yourfile20 << endl;
	}
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile20 << qz0[k][j] << ";";
		}
		yourfile20 << endl;
	}
	yourfile20.close();


	/////////////// Stric force and torque vectors

	vector <double> NewtonForceX(N); 
	vector <double> NewtonForceY(N);
	vector <double> NewtonForceZ(N);

	vector <double> NewtonTorkX(N); 
	vector <double> NewtonTorkY(N);
	vector <double> NewtonTorkZ(N);

	vector <vector <double> > Force_total_movazi(3, vector<double>(N, 0.0));
	vector <vector <double> > Force_total_amodi(3, vector<double>(N, 0.0));


	// initialization
	vector <double > Force_total = { 0,0,0 };
	vector <double > Force_temp = { 0,0,0 };
	vector <double > Tork_total = { 0,0,0 };
	vector <double > Tork_temp = { 0,0,0 };
	vector <double > Tork_temp2 = { 0,0,0 };
	


	double norm_omega;
	default_random_engine generator(pid); // To generate random forces, torques due to Brownian motion, 
	normal_distribution<double> qq1(0, sqrt(2 * D_movazi * dt));
	normal_distribution<double> qq2(0, sqrt(2 * D_amodi * dt));
	normal_distribution<double> qq3(0, sqrt(2 * D_amodi * dt));

	normal_distribution<double> qq5(0, sqrt(2 * D_r * dt));
	normal_distribution<double> qq6(0, sqrt(2 * D_r * dt));

	
	W[0] = 1 / sqrt(3);
	W[1] = 1 / sqrt(3);
	W[2] = 1 / sqrt(3);

	

	int k_save = 0;
	shomarande = 0;
	int i = 0;   // key index for the main loop below, always should remain 0



     ///// You may print the results without internal sampling as well
     
//	ofstream myfile0("Results_NonInternal_rotational.txt");
//	ofstream myfile("Results_NonInternal_translational.txt");
//	myfile0 << "zaman[k]" << ";" << "zamanD" << ";" << "mdr" << ";" << "D_rotation[k]" << ";" << "Drratio[k]" << ";" << "Dtransratio[k]" << ";" << endl;
// myfile << "zaman[k]" << ";" << "zamanD" << ";" << "msd_theory[k]" << ";" << "msd[k]" << ";" << "msd_internal[k]" << ";" << "slope[k]" << ";" << "slope_internal[k]" << ";" << "error% " << "; " << "Dtransratio " << "; " << "msd_movazi" <<"; " << "Ati" << ";" << endl;




	double delta=1e-11; // Very small number used for numerical integration
	double square_length_threshold = pow (L_mucin+ (L_lowen/2)+10e-9,2); // A early rejection if RNP and Fiber are too far away
	R_j_periodic= {0,0,0};



////////////////////////////////////////////////////////////////////////////// Heart of the code
cout << "************* Main loop started " << endl;
for (int k = 0; k < step; k++) {

	for (int j = 0; j < N; j++) { // RNP - Fiber interaction

            vector <double> R_i = { qx[i][j] , qy[i][j], qz[i][j] };  // center of mass of i-th particle
			vector <double> P_i = { omega[0][j] ,omega[1][j],omega[2][j] };  // oritentaion of i-th particle

			Force_temp[0] = 0;
			Force_temp[1] = 0;
			Force_temp[2] = 0;

			Tork_temp[0] = 0;
			Tork_temp[1] = 0;
			Tork_temp[2] = 0;

			for (int uu = 0; uu < Nmucin ; uu++) 
			{
				vector <double> R_j = { mucin_location[0][uu] , mucin_location[1][uu], mucin_location[2][uu] }; // center of mass of j th particle,i.e, FIBER
				vector <double> P_j = { mucin_omega[0][uu] ,mucin_omega[1][uu],mucin_omega[2][uu] }; // oritentaion of i-th particle

				for (int w1=-1;w1<2;w1++) 
				{
					for (int w2=-1;w2<2;w2++) 
					{
						for (int w3=-1;w3<2;w3++) 
						{
							
							R_j_periodic[0] = R_j[0] + (b * w1);     
							R_j_periodic[1] = R_j[1] + (b * w2);
							R_j_periodic[2] = R_j[2] + (b * w3);
														
							double square_distance= pow (R_j_periodic[0]- R_i[0],2)+ pow (R_j_periodic[1]- R_i[1] ,2) +pow (R_j_periodic[2]- R_i[2],2);  
							
							if (square_distance > square_length_threshold ) { // RNP-Fiber are far from each other, no overlap is plaussible 
								continue;
							}
							else 
							{ 
							
								LJ_force_torque_calc (R_i, P_i, R_j_periodic, P_j, Force_temp, Tork_temp, epsilon, delta, cutoff_LJ, L_switzer,L_mucin,d_mucin , p, uu,j, k);
	
								NewtonForceX[j] = NewtonForceX[j] + Force_temp[0];
								NewtonForceY[j] = NewtonForceY[j] + Force_temp[1];
								NewtonForceZ[j] = NewtonForceZ[j] + Force_temp[2];
		
		
								NewtonTorkX[j] = NewtonTorkX[j] + Tork_temp[0];
								NewtonTorkY[j] = NewtonTorkY[j] + Tork_temp[1];
								NewtonTorkZ[j] = NewtonTorkZ[j] + Tork_temp[2];
							}
						}
					}
				}
			}


			double TermQ = (NewtonForceX[j] * omega[0][j]) + (NewtonForceY[j] * omega[1][j]) + (NewtonForceZ[j] * omega[2][j]);
			Force_total_movazi[0][j] = TermQ * omega[0][j];
			Force_total_movazi[1][j] = TermQ * omega[1][j];
			Force_total_movazi[2][j] = TermQ * omega[2][j];

			Force_total_amodi[0][j] = NewtonForceX[j] - Force_total_movazi[0][j];
			Force_total_amodi[1][j] = NewtonForceY[j] - Force_total_movazi[1][j];
			Force_total_amodi[2][j] = NewtonForceZ[j] - Force_total_movazi[2][j];

			double rand_movazi = qq1(generator);

			// Updating RNP position and orientation based on Lowen's paper algorithm 
			qx_movazi[i + 1][j] = qx_movazi[i][j] + (rand_movazi * omega[0][j]) + ((D_movazi * dt / (k_B * T)) * Force_total_movazi[0][j]);
			qy_movazi[i + 1][j] = qy_movazi[i][j] + (rand_movazi * omega[1][j]) + ((D_movazi * dt / (k_B * T)) * Force_total_movazi[1][j]);
			qz_movazi[i + 1][j] = qz_movazi[i][j] + (rand_movazi * omega[2][j]) + ((D_movazi * dt / (k_B * T)) * Force_total_movazi[2][j]);

			double rand_amodi = qq2(generator);
			double rand_amodi2 = qq3(generator);


			qx_amodi[i + 1][j] = qx_amodi[i][j] + (rand_amodi * e1[0][j]) + (rand_amodi2 * e2[0][j]) + ((D_amodi * dt / (k_B * T)) * Force_total_amodi[0][j]);
			qy_amodi[i + 1][j] = qy_amodi[i][j] + (rand_amodi * e1[1][j]) + (rand_amodi2 * e2[1][j]) + ((D_amodi * dt / (k_B * T)) * Force_total_amodi[1][j]);
			qz_amodi[i + 1][j] = qz_amodi[i][j] + (rand_amodi * e1[2][j]) + (rand_amodi2 * e2[2][j]) + ((D_amodi * dt / (k_B * T)) * Force_total_amodi[2][j]);

			double rand_rotation = qq5(generator);
			double rand_rotation2 = qq6(generator);
			
			// Calculate external Torque
			vector<double> Ext_tork = crossProduct({ NewtonTorkX[j],NewtonTorkY[j],NewtonTorkZ[j] }, { omega[0][j] ,omega[1][j],omega[2][j] });

			omega[0][j] = omega[0][j] + (rand_rotation * e1[0][j]) + (rand_rotation2 * e2[0][j]) + ((D_r * dt / (k_B * T)) * Ext_tork[0]); // equation 10 of Lowen 1994
			omega[1][j] = omega[1][j] + (rand_rotation * e1[1][j]) + (rand_rotation2 * e2[1][j]) + ((D_r * dt / (k_B * T)) * Ext_tork[1]); // equation 10
			omega[2][j] = omega[2][j] + (rand_rotation * e1[2][j]) + (rand_rotation2 * e2[2][j]) + ((D_r * dt / (k_B * T)) * Ext_tork[2]); // equation 10


			norm_omega = sqrt((omega[0][j] * omega[0][j]) + (omega[1][j] * omega[1][j]) + (omega[2][j] * omega[2][j]));
			omega[0][j] = omega[0][j] / norm_omega;
			omega[1][j] = omega[1][j] / norm_omega;
			omega[2][j] = omega[2][j] / norm_omega;

			e1[0][j] = (omega[1][j] * W[2]) - (omega[2][j] * W[1]);
			e1[1][j] = (omega[2][j] * W[0]) - (omega[0][j] * W[2]);
			e1[2][j] = (omega[0][j] * W[1]) - (omega[1][j] * W[0]);


			if ((e1[0][j] == 0) && (e1[1][j] == 0) && (e1[2][j] == 0))
			{
			
				W[0] = -1 / sqrt(3);
				W[1] = -1 / sqrt(3);
				W[2] = -1 / sqrt(3);
				e1[0][j] = (omega[1][j] * W[2]) - (omega[2][j] * W[1]);
				e1[1][j] = (omega[2][j] * W[0]) - (omega[0][j] * W[2]);
				e1[2][j] = (omega[0][j] * W[1]) - (omega[1][j] * W[0]);
			}
			
			// e1 Must be normalzied
			double norm_e1 = sqrt((e1[0][j] * e1[0][j]) + (e1[1][j] * e1[1][j]) + (e1[2][j] * e1[2][j]));
			e1[0][j] = e1[0][j] / norm_e1;
			e1[1][j] = e1[1][j] / norm_e1;
			e1[2][j] = e1[2][j] / norm_e1;

			e2[0][j] = (omega[1][j] * e1[2][j]) - (omega[2][j] * e1[1][j]);
			e2[1][j] = (omega[2][j] * e1[0][j]) - (omega[0][j] * e1[2][j]);
			e2[2][j] = (omega[0][j] * e1[1][j]) - (omega[1][j] * e1[0][j]);


			// e2 Must be normalzied
			double norm_e2 = sqrt((e2[0][j] * e2[0][j]) + (e2[1][j] * e2[1][j]) + (e2[2][j] * e2[2][j]));
			e2[0][j] = e2[0][j] / norm_e2;
			e2[1][j] = e2[1][j] / norm_e2;
			e2[2][j] = e2[2][j] / norm_e2;

			qx[i + 1][j] = qx_movazi[i + 1][j] + qx_amodi[i + 1][j];
			qy[i + 1][j] = qy_movazi[i + 1][j] + qy_amodi[i + 1][j];
			qz[i + 1][j] = qz_movazi[i + 1][j] + qz_amodi[i + 1][j];

			
			
			// Peirodic boundary condition
			// x-axis
			if (qx[i + 1][j] > b)
			{
				qx[i + 1][j] = qx[i + 1][j] - b;
				qx_total[i + 1][j] = qx_total[i][j] + (qx[i + 1][j] - qx[i][j] + b);
				counter[k] = counter[k] + 1;
			}
			else if (qx[i + 1][j] < 0) {
				qx[i + 1][j] = qx[i + 1][j] + b;
				qx_total[i + 1][j] = qx_total[i][j] + (qx[i + 1][j] - qx[i][j] - b);
				counter[k] = counter[k] + 1;
			}
			else {
				qx_total[i + 1][j] = qx_total[i][j] + (qx[i + 1][j] - qx[i][j]);
			}


			// y-axis
			if (qy[i + 1][j] > b)
			{
				qy[i + 1][j] = qy[i + 1][j] - b;
				qy_total[i + 1][j] = qy_total[i][j] + (qy[i + 1][j] - qy[i][j] + b);
				counter[k] = counter[k] + 1;
			}
			else if (qy[i + 1][j] < 0) {
				qy[i + 1][j] = qy[i + 1][j] + b;
				qy_total[i + 1][j] = qy_total[i][j] + (qy[i + 1][j] - qy[i][j] - b);
				counter[k] = counter[k] + 1;
			}
			else {
				qy_total[i + 1][j] = qy_total[i][j] + (qy[i + 1][j] - qy[i][j]);
			}

			// z-axis
			if (qz[i + 1][j] > b)
			{
				qz[i + 1][j] = qz[i + 1][j] - b;
				qz_total[i + 1][j] = qz_total[i][j] + (qz[i + 1][j] - qz[i][j] + b);
				counter[k] = counter[k] + 1;
			}
			else if (qz[i + 1][j] < 0) {
				qz[i + 1][j] = qz[i + 1][j] + b;
				qz_total[i + 1][j] = qz_total[i][j] + (qz[i + 1][j] - qz[i][j] - b);
				counter[k] = counter[k] + 1;
			}
			else {
				qz_total[i + 1][j] = qz_total[i][j] + (qz[i + 1][j] - qz[i][j]);
			}



			// End of periodic BC
			msd[k + 1] = msd[k + 1] + pow(qx_total[i + 1][j] - qx0[0][j], 2) + pow(qy_total[i + 1][j] - qy0[0][j], 2) + pow(qz_total[i + 1][j] - qz0[0][j], 2);
			msd_movazi[k + 1] = msd_movazi[k + 1] + pow (((qx_total[i + 1][j] - qx0[0][j])*omega0[0][j])+((qy_total[i + 1][j] - qy0[0][j])*omega0[1][j])+((qz_total[i + 1][j] - qz0[0][j])*omega0[2][j]),2);
			mdr[k + 1] = mdr[k + 1] + ((omega[0][j] * omega0[0][j]) + (omega[1][j] * omega0[1][j]) + (omega[2][j] * omega0[2][j]));

			qx[i][j] = qx[i + 1][j];
			qx_total[i][j] = qx_total[i + 1][j];

			qx_movazi[i][j] = ((omega[0][j] * qx[i][j]) + (omega[1][j] * qy[i][j]) + (omega[2][j] * qz[i][j])) * omega[0][j];
			qx_amodi[i][j] = qx[i][j] - qx_movazi[i][j];

			qy[i][j] = qy[i + 1][j];
			qy_total[i][j] = qy_total[i + 1][j];
			qy_movazi[i][j] = ((omega[0][j] * qx[i][j]) + (omega[1][j] * qy[i][j]) + (omega[2][j] * qz[i][j])) * omega[1][j];
			qy_amodi[i][j] = qy[i][j] - qy_movazi[i][j];

			qz[i][j] = qz[i + 1][j];
			qz_total[i][j] = qz_total[i + 1][j];
			qz_movazi[i][j] = ((omega[0][j] * qx[i][j]) + (omega[1][j] * qy[i][j]) + (omega[2][j] * qz[i][j])) * omega[2][j];
			qz_amodi[i][j] = qz[i][j] - qz_movazi[i][j];

			if (fmod(k + 1, Saveeverystep) == 0) {

				save_qx[k_save + 1][j] = qx[i][j];
				save_qy[k_save + 1][j] = qy[i][j];
				save_qz[k_save + 1][j] = qz[i][j];

				save_qx_total[k_save + 1][j] = qx_total[i][j];
				save_qy_total[k_save + 1][j] = qy_total[i][j];
				save_qz_total[k_save + 1][j] = qz_total[i][j];

				save_omega[3 * (k_save + 1) + 0][j] = omega[0][j];
				save_omega[3 * (k_save + 1) + 1][j] = omega[1][j];
				save_omega[3 * (k_save + 1) + 2][j] = omega[2][j];


				if (j == N - 1) {  
					k_save++;
				}
			}



			if (k == progress[shomarande])
			{
				cout << "progress is: %" << progress[shomarande] * 100 / step << endl;
				shomarande = shomarande + 1;
			}





		} // End of loop for N Rod-like nanoparticles



		for (int j = 0; j < N; j++) {  // Force and torque must be set to ZERO  for the next iteraction
		
			NewtonForceX[j] = 0;
			NewtonForceY[j] = 0;
			NewtonForceZ[j] = 0;
			NewtonTorkX[j] = 0;
			NewtonTorkY[j] = 0;
			NewtonTorkZ[j] = 0;
		}
		
		// Ensemble averaging 
		msd[k] = msd[k] / N;
		msd_movazi[k] = msd_movazi[k] / N;
		mdr[k] = log(mdr[k] / N);
		msd_theory[k] = 6 * D0 * zaman[k];
		slope[k] = msd[k] / (6 * zaman[k]);
		Ati[k]= (3*msd_movazi[k]/msd[k] )-1;
		D_rotation[k] = mdr[k] / (-2 * zaman[k]);
		Drratio[k] = D_rotation[k] / D_r;
		Dtransratio[k] = slope[k] / D0;
		error[k] = ((slope[k] - D0) / D0) * 100; 

		// Print above results in txt file:
		//if (k%20==0) {  // Saving the results at every 1dt create very large txt files...
			//	myfile0 << zaman[k] << ";" << zamanD[k] << ";" << mdr[k] << ";" << D_rotation[k] << ";" << Drratio[k] << ";" << Dtransratio[k] << ";";
			//	myfile0 << endl;
	
			//	myfile << zaman[k] << ";" << zamanD[k] << ";" << msd_theory[k] << ";" << msd[k] << ";" << msd_internal[k] << ";" << slope[k] << ";" << slope_internal[k] << ";" << error[k] << ";" << Dtransratio[k] << ";" << msd_movazi[k] << ";" << Ati[k] << ";"  ;
			//	myfile << endl;		
		//}



} // Main loop is done




	// get the current time
	currentTime = time(0);
	// convert the current time to a string
	timeString = ctime(&currentTime);
	// print the current time
	cout << "End of Main loop : " << timeString << endl;
	
	
	///////////////////////////////////////////////////////////////////////////////////// Print results

	ofstream yourfile21("RNP_coordinates_end_simulation.txt");
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile21 << qx[k][j] << ";";
		}
		yourfile21 << endl;
	}
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile21 << qy[k][j] << ";";
		}
		yourfile21 << endl;
	}
	for (int k = 0; k < 1; k++) {
		for (int j = 0; j < N; j++) {
			yourfile21 << qz[k][j] << ";";
		}
		yourfile21 << endl;
	}
	yourfile21.close();



	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////// Internal sampling

	// get the current time
	currentTime = time(0);
	// convert the current time to a string
	timeString = ctime(&currentTime);
	// print the current time
	cout << "Internal sampling  started at: " << timeString << endl;
	// Print current time  End of block
	
	


	int numberOfdeltaT = save_qx_total.size();
	int left = numberOfdeltaT - floor(numberOfdeltaT / 2);

	for (int dt = 1; dt < numberOfdeltaT; dt++) { 
		double msd_dt = 0.0;
		double mdr_dt = 0.0;
		int numElements = 0;

		if (dt < (numberOfdeltaT / 2)) {

			for (int i = dt; i < left + dt; i++) {  
				numElements = numElements + 1;
				for (int j = 0; j < N; j++) {

					mdr_dt = mdr_dt + ((save_omega[3 * i][j] * save_omega[3 * (i - dt)][j]) + (save_omega[3 * i + 1][j] * save_omega[3 * (i - dt) + 1][j]) + (save_omega[3 * i + 2][j] * save_omega[3 * (i - dt) + 2][j]));
					msd_dt = msd_dt + (pow(save_qx_total[i][j] - save_qx_total[i - dt][j], 2) + pow(save_qy_total[i][j] - save_qy_total[i - dt][j], 2) + pow(save_qz_total[i][j] - save_qz_total[i - dt][j], 2));
					
				}
			
			}
		}
		
		else {
			numElements = numberOfdeltaT - dt;
			for (int i = dt; i < numberOfdeltaT; i++) {
				for (int j = 0; j < N; j++) {
					
					mdr_dt = mdr_dt + ((save_omega[3 * i][j] * save_omega[3 * (i - dt)][j]) + (save_omega[3 * i + 1][j] * save_omega[3 * (i - dt) + 1][j]) + (save_omega[3 * i + 2][j] * save_omega[3 * (i - dt) + 2][j]));
					msd_dt = msd_dt + (pow(save_qx_total[i][j] - save_qx_total[i - dt][j], 2) + pow(save_qy_total[i][j] - save_qy_total[i - dt][j], 2) + pow(save_qz_total[i][j] - save_qz_total[i - dt][j], 2));
				
				}

			}
		}
		msd_internal[dt] = msd_dt / numElements;
		mdr_internal[dt] = mdr_dt / numElements;
	}



// Compute results under internal sampling
	for (int k = 1; k < SaveN; k++) {
		msd_internal[k] = msd_internal[k] / N;
		mdr_internal[k] = log(mdr_internal[k] / N);
		slope_internal[k] = msd_internal[k] / (6 * zaman_internal[k]);
		D_rotation_internal[k] = mdr_internal[k] / (-2 * zaman_internal[k]);
		Drratio_internal[k] = D_rotation_internal[k] / D_r;
		Dtransratio_internal[k] = slope_internal[k] / D0;
	}


	ofstream myfile99("Internal_TranslationRotation.txt");
	myfile99 << "zamanD_internal" << ";" << "Drratio_internal " << "; " << " Dtransratio_internal" << "; " << endl;
	for (int k = 0; k < SaveN; k++) {
		myfile99 << zamanD_internal[k] << ";" << Drratio_internal[k] << ";" << Dtransratio_internal[k] << ";";
		myfile99 << endl;
	}
	myfile99.close();



	// get the current time
	currentTime = time(0);
	// convert the current time to a string
	timeString = ctime(&currentTime);
	cout << "Internal sampling finished at: " << timeString << endl;
	


	// Print qX-save
	ofstream myfile90("save_qx.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile90 << save_qx[k][j] << ";";
		}
		myfile90 << endl;
	}
	myfile90.close();

	// Print qY-save
	ofstream myfile91("save_qy.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile91 << save_qy[k][j] << ";";
		}
		myfile91 << endl;
	}
	myfile91.close();

	// Print qZ-save
	ofstream myfile92("save_qz.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile92 << save_qz[k][j] << ";";
		}
		myfile92 << endl;
	}
	myfile92.close();

	// Print save_qx_total
	ofstream myfile93("save_qx_total.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile93 << save_qx_total[k][j] << ";";
		}
		myfile93 << endl;
	}
	myfile93.close();

	// Print save_qy_total
	ofstream myfile94("save_qy_total.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile94 << save_qy_total[k][j] << ";";
		}
		myfile94 << endl;
	}
	myfile94.close();

	// Print save_qz_total
	ofstream myfile95("save_qz_total.txt");
	for (int k = 0; k < SaveN; k++) {
		for (int j = 0; j < N; j++) {
			myfile95 << save_qz_total[k][j] << ";";
		}
		myfile95 << endl;
	}
	myfile95.close();

	// Print Saveomega
	ofstream myfile96("save_omega.txt");
	for (int k = 0; k < 3 * SaveN; k++) {
		for (int j = 0; j < N; j++) {

			myfile96 << save_omega[k][j] << ";";


		}
		myfile96 << endl;
	}
	myfile96.close();


	// get the current time
	currentTime = time(0);
	// convert the current time to a string
	timeString = ctime(&currentTime);
	// print the current time
	cout << "Simulation finished at: " << timeString << endl;


	auto end = chrono::high_resolution_clock::now();  // // chat gpb this part
	auto duration = chrono::duration_cast<chrono::seconds>(end - start);
	cout << "Execution time: " << duration.count() << " seconds\n";

	return 0;

}








