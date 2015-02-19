#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>

using namespace std;

double cal_dis(double x1, double y1, double z1, double x2, double y2, double z2){
	double dx= x1 - x2;
	double dy= y1 - y2;
	double dz= z1 - z2;

	return sqrt(dx*dx + dy*dy + dz*dz);
}

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
#define MAX_TYPE 10
#define DEF_CLTR 9 // The definition of minimum number of ltcps for a cluster
const double vbra[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
const int n1nbr= 8;
const int v1nbr[8][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1},
			{-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1},
			{ 1,  1,  1}, {-1, -1, -1}	       };

const int nx= 64;
const int ny= 64;
const int nz= 64;
// parameters //

// global variables
double realtime= 0; // real time in the simulation (unit: s) 
int timestep= 0;
int states[nx*ny*nz];
double x[nx*ny*nz];
double y[nx*ny*nz];
double z[nx*ny*nz];

double csd[nx*ny*nz]= {0};

FILE * out_csd; 
// global variables

//######################### functions #########################//
// Read files
void read_xyz(char ifname[]);
// Calculations
void cal_csd();
//######################### functions #########################//

int main(int nArg, char *Arg[]){
	int t0cpu= time(0);

	cout << "Calculation is starting..." << endl;
	cout << "The system size is " << nx << " x " << ny << " x " << nz << endl;

	// READING FILES
	if(nArg != 2) error(0, "nArg must be 2\nUse csd.exe <.xyz>\n");
	read_xyz(Arg[1]);
	// READING FILES
	
	// OPEN OUTPUT FILES  
	out_csd= fopen("out.csd", "w");
	if(NULL==out_csd) error(1, "out_csd was not open");
	// OPEN OUTPUT FILES  
	
	// CALCULATIONS
	cal_csd();
	// CALCULATIONS

	int tfcpu= time(0);
	cout << "\nCalculation is completed!! Total CPU time: " << tfcpu - t0cpu << " secs" << endl;
}

void read_xyz(char ifname[]){ // Reading t0.xyz ////////////////////
	ifstream in_t0(ifname, ios::in);
	if(! in_t0.is_open()) error(1, "t0.xyz was not open"); // check

	int ntotal; in_t0 >> ntotal; in_t0.ignore(); 
	if(ntotal != nx*ny*nz) error(1, "ntotal != nx*ny*nz", 2, ntotal, nx*ny*nz); // check
	
	string line2; getline(in_t0, line2);

	for(int a=0; a<ntotal; a ++){
		if(in_t0.eof()) error(1, "reach end of file before finish reading all data");
	
		in_t0 >> states[a] >> x[a] >> y[a] >> z[a];
	}
	in_t0.close();
	cout << "t0.xyz file reading completed" << endl;
}

void cal_csd(){
	vector < vector <double> > ximg(nx*ny*nz, vector<double>(27)); // coordinates of all atoms in image boxes
	vector < vector <double> > yimg(nx*ny*nz, vector<double>(27)); // coordinates of all atoms in image boxes
	vector < vector <double> > zimg(nx*ny*nz, vector<double>(27)); // coordinates of all atoms in image boxes

	cout << "Now calculating CSD..." << endl;

	// construct image boxes
	for(int i=0; i<nx*ny*nz; i ++){ // put a=0, b=0, c=0 in [0] so that it scans first and faster
		ximg[i][0]= x[i];
		yimg[i][0]= y[i];
		zimg[i][0]= z[i];
	}
	int i_image= 0;
	for(int a=-1; a<=1; a++){
		for(int b=-1; b<=1; b++){
			for(int c=-1; c<=1; c++){
				if((a==0) && (b==0) && (c==0)) continue;
				i_image ++;
				
				for(int i=0; i<nx*ny*nz; i ++){
					ximg[i][i_image]= x[i] + a*nx*vbra[0][0] + b*ny*vbra[1][0] + c*nz*vbra[2][0];
					yimg[i][i_image]= y[i] + a*nx*vbra[0][1] + b*ny*vbra[1][1] + c*nz*vbra[2][1];
					zimg[i][i_image]= z[i] + a*nx*vbra[0][2] + b*ny*vbra[1][2] + c*nz*vbra[2][2];
				}
	}}}
	cout << "image boxes construction completed" << endl;

	// calculate distance and get 8 vectors
	cout << "calculating csd:" << endl;
	for(int i= 0; i<nx*ny*nz; i++){
#define TOL 0.01
		double xi= x[i];
		double yi= y[i];
		double zi= z[i];
		vector < vector<double> > vec(3, vector<double>(0));
		vector <int> statej; // states of j

		for(int a= 0; a<27; a++){ // search for 8 nearest neighbors
			for(int j= 0; j<nx*ny*nz; j++){
				double xj= ximg[j][a];
				double yj= yimg[j][a];
				double zj= zimg[j][a];

				double dis= cal_dis(xi, yi, zi, xj, yj, zj);

				if((dis<(0.866+TOL)) && (dis>(0+TOL))){
					vec[0].push_back(xj-xi);
					vec[1].push_back(yj-yi);
					vec[2].push_back(zj-zi);
					statej.push_back(states[j]);
				}

				if(8==vec[0].size()) goto search_end;
			}
		}
		cout << "(xi, yi, zi): " << xi << " " << yi << " " << zi << endl;
		error(1, "can't find 8 neighbors: ", 1, vec[0].size());
search_end: ;
	   	while(vec[0].size() !=0){
	   		double min_vdiff2= 99999; // smallest diff(v1-v2)**2
			int    min_id= 0;
			for(int a= 1; a<vec[0].size(); a ++){
				double vdiff2= pow((vec[0][0]+vec[0][a]), 2) + pow((vec[1][0]+vec[1][a]), 2) + pow((vec[2][0]+vec[2][a]), 2);
				if(vdiff2<min_vdiff2){
					min_vdiff2= vdiff2;
					min_id    = a;
				}
			}

			if((0==statej[0]) && (0==statej[min_id])) ;
			else if(0==statej[0])		csd[i] += (vec[0][min_id]*vec[0][min_id] + vec[1][min_id]*vec[1][min_id] + vec[2][min_id]*vec[2][min_id]);
			else if(0==statej[min_id])	csd[i] += (vec[0][0]*vec[0][0]		 + vec[1][0]*vec[1][0]		 + vec[2][0]*vec[2][0]);
			else				csd[i] += min_vdiff2;
			
			vec[0].erase(vec[0].begin()+min_id); vec[0].erase(vec[0].begin()); 
			vec[1].erase(vec[1].begin()+min_id); vec[1].erase(vec[1].begin());
			vec[2].erase(vec[2].begin()+min_id); vec[2].erase(vec[2].begin()); 
			statej.erase(statej.begin()+min_id); statej.erase(statej.begin());
		}

		if((i+1)%(ny*nz)==0) cout << "Progress: " << i+1 << "/" << nx*ny*nz << endl;
	}

#define MLTPLR 10
	int sum_csd[MLTPLR*5]={0};
	for(int i=0; i<nx*ny*nz; i ++){
		int int_csd;
		if(0==csd[i])	int_csd=0;
		else		int_csd= (int) (csd[i]*MLTPLR+1);
		if(int_csd>=MLTPLR*5) error(1, "int_csd too large", 1, int_csd);

		sum_csd[int_csd] ++;
	}

	fprintf(out_csd, "%f %d\n", 0.0, sum_csd[0]);
	for(int a=1; a<MLTPLR*5; a ++) fprintf(out_csd, "%f %d\n", (a-0.5)/MLTPLR, sum_csd[a]);
}
