#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>

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
const int nx= 256;
const int ny= 64;
const int nz= 64;

int Ttype= -1;
// parameters //

// global variables
long long int timestep; 
double realtime;
int states[nx][ny][nz];

int N_check= 0;
int den[nx]= {0};

FILE * out_den; 
// global variables

//######################### functions #########################//
// Read files
void read_xyz(char ifname[]);
// Calculations
void cal_den();
//######################### functions #########################//

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	cout << "Calculation is starting..." << endl;
	cout << "The system size is " << nx << " x " << ny << " x " << nz << endl;

	// READING FILES
	if(nArg != 2) error(0, "nArg must be 2\nUse den.exe <.xyz>\n");
	read_xyz(Arg[1]);
	// READING FILES
	
	// OPEN OUTPUT FILES 
	char name[50]="out.", s_ts[50]; 
	sprintf(s_ts, "%lld", timestep);
	strcat(name, s_ts);
	
	out_den= fopen(name, "w");
	if(NULL==out_den) error(1, "out_den was not open");
	// OPEN OUTPUT FILES  
	
	// CALCULATIONS
	cal_den();
	// CALCULATIONS

	int tfcpu= time(0);
	cout << "\nCalculation is completed!! Total CPU time: " << tfcpu - t0cpu << " secs" << endl;
}

void read_xyz(char ifname[]){ // Reading .ltcp ////////////////////
	ifstream in_t0(ifname, ios::in);
	if(! in_t0.is_open()) error(1, "t0.xyz was not open"); // check

	int ntotal; in_t0 >> ntotal;
	if(ntotal != nx*ny*nz) error(1, "ntotal != nx*ny*nz", 2, ntotal, nx*ny*nz); // check

	char dump[50];
	in_t0 >> dump >> timestep >> realtime;
	cout << "the file is at " << timestep << " " << realtime << " (timestep time)" << endl;

	for(int a=0; a<ntotal; a ++){
		if(in_t0.eof()) error(1, "reach end of file before finish reading all data");

		int state_in, x, y, z, ix, iy, iz, dir, head;
		in_t0 >> state_in >> x >> y >> z;

		if(state_in != 1 && state_in != -1) in_t0 >> ix >> iy >> iz;
		if(state_in == 2 || state_in == -2 || state_in == 3) in_t0 >> dir >> head;

		states[x][y][z]= state_in;
		if(Ttype==state_in) N_check ++;
	}
	in_t0.close();
	cout << "t0.xyz file reading completed" << endl;
}

void cal_den(){
	int N_Ttype= 0;

	for(int i=0; i<nx; i ++){
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				if(Ttype==states[i][j][k]){
					N_Ttype ++;
					den[i] ++;
				}
	}}}

	if(N_check != N_Ttype) error(2, "N Ttype inconsistent", N_check, N_Ttype);
	cout << "number of targeted type: " << N_Ttype << endl;

	double dis= 0.5*sqrt(3);
	fprintf(out_den, "%f %f %lld %e\n", -(nx/2.0)*dis, den[0]*1.0/ny/nz, timestep, realtime);
	for(int i=1; i<nx; i ++) fprintf(out_den, "%f %f\n", (i-nx/2.0)*dis, den[i]*1.0/ny/nz);
}
