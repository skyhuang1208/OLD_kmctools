#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>

using namespace std;

void error(int nexit, string errinfo, int nnum=0, double num1=0, double num2=0){
	cout << "\nError: ";
	cout << errinfo << " ";
	switch(nnum){
		case 0:	 cout << endl;			      break;
		case 1:  cout << num1 << endl; 		      break;
		case 2:  cout << num1 << " " << num2 << endl; break;
		default: cout << "!!!ERROR FUNCTION MALFUNCTION!!! WRONG NNUM!!!" << endl;
	}
	exit(1); 
}


// parameters //
#define PI 3.14159265358979323846
const double vbra[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
const bool is_images= true;
const int  N_images= 2; // -N_images ~ N_images in each dimension
const int  N_totbox= pow((2*N_images+1), 3); // total boxes

const int nx=  10;
const int ny=  10;
const int nz=  10;

#define SUB 0.0000000001 // a small number to avoid error
const double r0= 0.1 - SUB; // smallest sphere radius of rdf
double dr; // width of the shell
const int nr= 1000; // number of sphere shells

const double sigma= 0.1; // standard deviation of the gaussian function: vibration distance 

const char name_in_t0[20]=  "./t0.ltcp";
// parameters //

// global variables
const int nall= nx*ny*nz;
int states[nall];
int x[nall];
int y[nall];
int z[nall];

FILE * out;
FILE * out_nor;
// global variables

//######################### functions #########################//
// Read files
void read_t0();
// Calculation
double cal_rho();
void cal_rdf(double rho);
// Sub-functions
double cal_dis(int a1, int b1, int c1, int a2, int b2, int c2);
int pbc(int x_, int nx_); // Periodic Boundary Condition
void cal_dis_i(double dis[], int a1, int b1, int c1, int a2, int b2, int c2);
double gaussian(double x_);
//######################### functions #########################//

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	cout << "Calculation is starting..." << endl;
	cout << "The system size is " << nx << " x " << ny << " x " << nz << endl;
	if(nArg != 2) error(1, "nArg != 2, please enter the dr by ./a.out [dr]");
	dr= atof(Arg[1]);
	cout << "dr imported is: " << dr << endl;
	double rho= cal_rho();
	cout << "The density is: " << rho << endl;

	// OPEN OUTPUT FILES  
	char name[7]  = "t0.rdf";
	out= fopen(name, "w");
	if(NULL==out) error(1, "out file was not open");

	char name_nor[11]  = "t0_nor.rdf";
	out_nor= fopen(name_nor, "w");
	if(NULL==out_nor) error(1, "nor out file was not open");
	// CALCULATIONS
	read_t0();
	cal_rdf(rho);

	int tfcpu= time(0);
	cout << "\nCalculation is completed!! Total CPU time: " << tfcpu - t0cpu << " secs" << endl;
}

double cal_rho(){ // calculate the density of the whole system
	double v1[3]= {nx*vbra[0][0], nx*vbra[0][1], nx*vbra[0][2]};
	double v2[3]= {ny*vbra[1][0], ny*vbra[1][1], ny*vbra[1][2]};
	double v3[3]= {nz*vbra[2][0], nz*vbra[2][1], nz*vbra[2][2]};

	double cross[3]= {v2[1]*v3[2]-v2[2]*v3[1], v2[2]*v3[0]-v2[0]*v3[2], v2[0]*v3[1]-v2[1]*v3[0]};

	return nall / (v1[0]*cross[0] + v1[1]*cross[1] + v2[2]*cross[2]);
}

void read_t0(){ // Reading t0.ltcp
	ifstream in_t0(name_in_t0, ios::in);
	if(! in_t0.is_open()) error(1, "in_t0 was not open"); // check

	int ntotal; in_t0 >> ntotal; in_t0.ignore(); 
	if(ntotal != nx*ny*nz) error(1, "ntotal != nx*ny*nz", 2, ntotal, nx*ny*nz); // check
	
	string line2; getline(in_t0, line2);

	for(int a=0; a<ntotal; a ++){
		if(in_t0.eof()) error(1, "reach end of file before finish reading all data");
		
		int state_in;
		int i, j, k;
		in_t0 >> state_in >> i >> j >> k;

		states[a]= state_in;
		x[a]= i;
		y[a]= j;
		z[a]= k;
	}
	cout << "t0.ltcp file reading completed" << endl;
	in_t0.close();
}

void cal_rdf(double rho){
	cout << "Now calculating rdf ..." << endl;

	double sum[nr]= {0}; // sum of all ltcp included within a certain shell
	for(int i=0; i<nall; i ++){
		for(int j=0; j<nall; j ++){
			double dis[N_totbox];
			if(is_images)
				cal_dis_i(dis, x[i], y[i], z[i], x[j], y[j], z[j]);
			else
				dis[0]= cal_dis(x[i], y[i], z[i], x[j], y[j], z[j]);

			for(int a=0; a < nr; a++){
//				double r1= r0 + a*dr;
				double r_shell= r0 + 0.5*dr + a*dr;
				
				if(is_images){
					for(int b=0; b<N_totbox; b ++){
						if(0==dis[b]) continue;
						sum[a] += gaussian(dis[b]-r_shell);
//						if((dis[b] >= r1) && (dis[b] < r1+dr)) sum[a] ++;
					}
				}
				else
						if(0==dis[0]) continue;
						sum[a] += gaussian(dis[0]-r_shell);
//						if((dis[0] >= r1) && (dis[0] < r1+dr)) sum[a] ++;

			}
		
			printf("\r%d", i);
		}

	}
	cout << "\nRDF calculations completed!! Now output to the file..." << endl; 

	for(int a=0; a < nr; a++){
		double r_print= r0 + 0.5*dr + a*dr;
		double r1= r0 + dr* a;
		double r2= r0 + dr*(a+1);
		double volume = 4*PI/3*(pow(r2,3) - pow(r1,3));
		fprintf(out, "%f %f\n", r_print, sum[a] / nall);
		fprintf(out_nor, "%f %f\n", r_print, sum[a] / nall / volume / rho);
	}
}

double cal_dis(int a1, int b1, int c1, int a2, int b2, int c2){
	// index pbc distance
	int da= pbc((a1-a2), nx);
	int db= pbc((b1-b2), ny);
	int dc= pbc((c1-c2), nz);

	// distance in Cartesian Coordinate
	double dx= da*vbra[0][0] + db*vbra[1][0] + dc*vbra[2][0];
	double dy= da*vbra[0][1] + db*vbra[1][1] + dc*vbra[2][1];
	double dz= da*vbra[0][2] + db*vbra[1][2] + dc*vbra[2][2];

	return sqrt(dx*dx + dy*dy + dz*dz);
}

int pbc(int dx_, int nx_){ // Periodic Boundary Condition
	if (dx_<=-nx_ || dx_>=nx_) error(1, "(pbc) diff value out of bound", 2, dx_, nx_);

	if (abs(dx_) > nx_/2){ 
		if(dx_>0) return (dx_-nx_);
		else	  return (dx_+nx_);
	}
	else		  return  dx_;
}

// includes periodic image atoms
void cal_dis_i(double dis[], int a1, int b1, int c1, int a2, int b2, int c2){ // index pbc distance
	int index= -1;
	for(int a=-N_images; a<=N_images; a++){
		for(int b=-N_images; b<=N_images; b++){
			for(int c=-N_images; c<=N_images; c++){
				index ++;
				
				int a2i= a2 + a*nx;
				int b2i= b2 + b*ny;
				int c2i= c2 + c*nz;
				double dx= (a1-a2i)*vbra[0][0] + (b1-b2i)*vbra[1][0] + (c1-c2i)*vbra[2][0];
				double dy= (a1-a2i)*vbra[0][1] + (b1-b2i)*vbra[1][1] + (c1-c2i)*vbra[2][1];
				double dz= (a1-a2i)*vbra[0][2] + (b1-b2i)*vbra[1][2] + (c1-c2i)*vbra[2][2];

				dis[index]= sqrt(dx*dx + dy*dy + dz*dz);
	}}}
}

double gaussian(double x_){
	if(abs(x_) > 3*sigma)	return 0;
	else 			return dr * ( exp(-x_*x_/2/sigma/sigma) / sqrt(2*PI*sigma*sigma) );
}
