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
#define SUB 0.0000000001 // a small number to avoid error
const double vbra[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
const double sigma= 0.05; // standard deviation of the gaussian function: vibration distance 
const int  N_totbox= 8; // total boxes
const double r0= 0.1 - SUB; // smallest sphere radius of rdf

double dr= 0.01; 
double r_trunc;
int nr; 

int type_A;
int type_B;

const int nx=  64;
const int ny=  64;
const int nz=  64;
const int nall= nx*ny*nz;

char name_in_t0[20];
// parameters //

// global variables
int states[nx][ny][nz];
vector < vector<int> > xlist;
vector < vector<int> > ylist;
vector < vector<int> > zlist;
vector < vector<double> > glist;
vector < vector<int> > Alist(3);

FILE * out;
FILE * out_nor;
// global variables

//######################### functions #########################//
double cal_rho();
void read_t0();
void make_cal_lists();
void cal_rdf(double rho);
// Sub-functions for making cal lists
void cal_dis(double dis[], int a2, int b2, int c2);
double gaussian(double x_);
// Sub-function for cal rdf
int pbc(int x_, int nx_);
//######################### functions #########################//

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	cout << "Calculation is starting..." << endl;
	cout << "The system size is " << nx << " x " << ny << " x " << nz << endl;

	// READ IMPORTED PARAMETERS
	if(nArg != 5){
		cout << "Calculate RDF by using linked-cell method" << endl;
		cout << "WARNING!! have pre-assumed ltc are fixed BCC" << endl;
		cout << "Error: Parameter required: <name_in_ltcp> <truncation radius> <A_TYPE> <B_TYPE>" << endl;
		exit(1);
	}
	strcpy(name_in_t0, Arg[1]);
	r_trunc= atof(Arg[2]);
	type_A = atoi(Arg[3]);
	type_B = atoi(Arg[4]);
	cout << "\nImported parameters:" << endl << "dr= " << dr << ", truncation radius= " << r_trunc << endl;
	cout << "Atom type A is: " << type_A << ", and type B is: " << type_B << endl;

	// CAL PARAMETERS
	nr= r_trunc/dr + 10; // much larger than actually need 
	double rho= cal_rho();
	cout << "\nCalculated parameters: " << endl << "nr= " << nr << ", density is: " << rho << endl;

	// INITIALIZE VECTOR ARRAYS
	for(int i=0; i<nr; i ++){
		xlist.push_back(vector<int>());
		ylist.push_back(vector<int>());
		zlist.push_back(vector<int>());
		glist.push_back(vector<double>());
	}

	// OPEN OUTPUT FILES  
	char name[13]  = "out_nonn.rdf";
	out= fopen(name, "w");
	if(NULL==out) error(1, "out file was not open");

	char name_nor[8]  = "out.rdf";
	out_nor= fopen(name_nor, "w");
	if(NULL==out_nor) error(1, "nor out file was not open");
	
	// CALCULATIONS
	make_cal_lists();
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

	if(r_trunc > nx*0.5*sqrt(3)) error(1, "truncation radius is larger than the system size");

	return nall / (v1[0]*cross[0] + v1[1]*cross[1] + v2[2]*cross[2]);
}

void make_cal_lists(){
	cout << "\nNow making calculation list..." << endl;

	int index= -1;
	for(int i=0; i<nx; i ++){
		cout << i << "/" << nx << endl;
		for(int j=0; j<ny; j ++){
			for(int k=0; k<nz; k ++){
				double dis[N_totbox];
				cal_dis(dis, i, j, k);

				for(int a=0; a<N_totbox; a ++){
					if(0==dis[a]) continue;
					if(dis[a]>r_trunc+1) continue; // seems redundant; just simply use the system size as the truncation radius?

					for(int b=0; b < nr; b++){
						double r_shell= r0 + b*dr + 0.5*dr;
						double gvalue= gaussian(dis[a]-r_shell);

						if(gvalue != 0){
							xlist[b].push_back(i);
							ylist[b].push_back(j);
							zlist[b].push_back(k);
							glist[b].push_back(gvalue);
						}
					}
				}
	}}}

	cout << "\nList making completed!!" << endl;
}

void cal_dis(double dis[], int a2, int b2, int c2){ // index pbc distance from the point (0, 0, 0)
	int index= -1;
	for(int a=-1; a<=0; a++){
		for(int b=-1; b<=0; b++){
			for(int c=-1; c<=0; c++){
				index ++;
				
				int a2i= a2 + a*nx;
				int b2i= b2 + b*ny;
				int c2i= c2 + c*nz;
				double dx= (0-a2i)*vbra[0][0] + (0-b2i)*vbra[1][0] + (0-c2i)*vbra[2][0];
				double dy= (0-a2i)*vbra[0][1] + (0-b2i)*vbra[1][1] + (0-c2i)*vbra[2][1];
				double dz= (0-a2i)*vbra[0][2] + (0-b2i)*vbra[1][2] + (0-c2i)*vbra[2][2];

				dis[index]= sqrt(dx*dx + dy*dy + dz*dz);
	}}}
}

double gaussian(double x_){
	if(abs(x_) > 3*sigma)	return 0;
	else 			return dr * ( exp(-x_*x_/2/sigma/sigma) / sqrt(2*PI*sigma*sigma) );
}

void read_t0(){ // Reading t0.ltcp
	cout << "\nNow reading t0.ltcp ...";
	
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

		states[i][j][k]= state_in;
		if(type_A==state_in){
			Alist[0].push_back(i);
			Alist[1].push_back(j);
			Alist[2].push_back(k);
		}
	}
	cout << "t0.ltcp file reading completed" << endl;
	in_t0.close();
}

void cal_rdf(double rho){
	cout << "\nNow calculating rdf ..." << endl;

	double sum[nr]; // sum of all ltcp included within a certain shell
	for(int a=0; a < nr; a++) sum[a]= 0;
	
	for(int i=0; i<Alist[0].size(); i ++){
		for(int a=0; a < nr; a++){
			for(int b=0; b<xlist[a].size(); b ++){
				int x2= pbc(Alist[0].at(i)+xlist[a].at(b), nx);
				int y2= pbc(Alist[1].at(i)+ylist[a].at(b), ny);
				int z2= pbc(Alist[2].at(i)+zlist[a].at(b), nz);

				if(type_B==states[x2][y2][z2]) sum[a] += glist[a].at(b);
			}
		}
		
		if(0==i%1000) cout << i << "/" << Alist[0].size() << endl;
	
	}

	cout << "\nRDF calculations completed!! Now output to the file..." << endl; 

	for(int a=0; a < nr; a++){
		double r_print= r0 + 0.5*dr + a*dr;
		double r1= r0 + dr* a;
		double r2= r0 + dr*(a+1);
		double volume = 4*PI/3*(pow(r2,3) - pow(r1,3));
		fprintf(out, "%f %f\n", r_print, sum[a] / Alist[0].size());
		fprintf(out_nor, "%f %f\n", r_print, sum[a] / Alist[0].size() / volume / rho);
	}
}

int pbc(int x_, int nx_){ // Periodic Boundary Condition
	if	(x_<-nx_ || x_>=2*nx_)
		error(1, "(pbc) input x is out of bound", 2, x_, nx_);

	if	(x_<0) 		return (x_ + nx_);
	else if	(x_<nx_)	return  x_;
	else			return (x_ - nx_);
}
