/* 	Calculate cluster related properties:
	1. csize: cluster counts with sizes in a specific region (1, >1, >5, >10, >30)
	2. csave: cluster size average (Define a cluster when # of atoms is more than 9)
	3. csind: cluster individual: all cluster sizes at several steps
	4.   lce: local chemical environment, percentage of solute atoms surrounding the vacancy during a period of time 
	5.   sro: short range order, 1-(zij/z/cj)
	6.   lro: long range order, abs[(-1)^i * sigma_i] / N

	Code structure:
	1. read_t0
	2. read_history
		a. update_cid(loop)
		b. sum csize
			(1) output csize
			(2) output csave
			(3) output csind
		c. cal sro
		d. cal lro
		e. cal lce

*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>

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

int pbc(int x_, int nx_){ // Periodic Boundary Condition
	if	(x_<-nx_ || x_>=2*nx_)
		error(1, "(pbc) input x is out of bound", 2, x_, nx_);

	if	(x_<0) 		return (x_ + nx_);
	else if	(x_<nx_)	return  x_;
	else			return (x_ - nx_);
}

// parameters //
#define MAX_TYPE 10
#define DEF_CLTR 9 // The definition of minimum number of ltcps for a cluster
const double vbra[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
const int n1nbr= 8;
const int v1nbr[8][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1},
			{-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1},
			{ 1,  1,  1}, {-1, -1, -1}	       };

const int nx=  63;
const int ny=  63;
const int nz=  63;

const int Ttype= -1;

const int sample_cltr=		1e5; 
const int cycle_out_csind=	1e7;

const int sample_msd=		1e3;
const int sample_sro=		1e3;
const int sample_lro=		1e3;

const int sample_lce=		1e3;
const int cycle_lce=		1e6; // MC steps that the period of the lce calculations (output an point of lce)

const int itype_sro= -1;
const int jtype_sro= 1;

const char name_in_t0[50]=  "./t0.ltcp";
const char name_in_his[50]= "./history.info";
// parameters //

// global variables
double realtime= 0; // real time in the simulation (unit: s) 
int timestep= 0;
int states[nx][ny][nz];
int ntt= 0; 
int cid= 0;
int id_cltr[nx*ny*nz]= {0}; // indicate the cluster id that a Ttype ltc(index) is labeled with (0 means not Ttype)
int v0[3];
int vpos[3];
int vi[3]= {0}; // image box of the vacancy

FILE * out_csize; // in sum_csize
FILE * out_csave; // in sum_csize
FILE * out_csind; // in write_csind
FILE * out_msd;   // in cal_msd
FILE * out_lce;   // in cal_lce
FILE * out_sro;   // in cal_sro
FILE * out_lro;   // in cal_lro
// global variables

//######################### functions #########################//
// Read files
void read_t0();
void read_his_cal();
// Calculations
void update_vpos(int iltcp);
void update_cid(int x, int y, int z); // input one ltc point and renew cluster id around it
void sum_csize();
void write_csind(int cid, int *N_in_cltr);
void cal_msd();
void cal_lce();
void cal_sro();
void cal_lro();
//######################### functions #########################//

int main(){
	int t0cpu= time(0);

	cout << "Calculation is starting..." << endl;
	cout << "The TARGETED type is " << Ttype << endl;
	cout << "The system size is " << nx << " x " << ny << " x " << nz << endl;

	// OPEN OUTPUT FILES  
	char name_csize[20]= "out.csize";
	char name_csave[20]= "out.csave";
	char name_csind[20]= "out.csind";
	char name_msd[20]  = "out.vmsd";
	char name_lce[20]  = "out.lce";
	char name_sro[20]  = "out.sro"; 
	char name_lro[20]  = "out.lro";
	
	out_csize= fopen(name_csize, "w");
	out_csave= fopen(name_csave, "w");
	out_csind= fopen(name_csind, "w");
	out_msd=   fopen(name_msd,   "w");
	out_lce=   fopen(name_lce,   "w");
	out_sro=   fopen(name_sro,   "w");
	out_lro=   fopen(name_lro,   "w");

	if(NULL==out_csize) error(1, "out_csize was not open");
	if(NULL==out_csave) error(1, "out_csave was not open");
	if(NULL==out_csind) error(1, "out_csind was not open");
	if(NULL==out_msd)   error(1, "out_msd   was not open");
	if(NULL==out_lce)   error(1, "out_lce   was not open");
	if(NULL==out_sro)   error(1, "out_sro   was not open");
	if(NULL==out_lro)   error(1, "out_lro   was not open");
	// OPEN OUTPUT FILES  

	// CALCULATIONS
	read_t0();
	read_his_cal();
	// CALCULATIONS

	int tfcpu= time(0);
	cout << "\nCalculation is completed!! Total CPU time: " << tfcpu - t0cpu << " secs" << endl;
}

void read_t0(){ // Reading t0.ltcp ////////////////////
	ifstream in_t0(name_in_t0, ios::in);
	if(! in_t0.is_open()) error(1, "in_t0 was not open"); // check

	int ntotal; in_t0 >> ntotal; 
	if(ntotal != nx*ny*nz) error(1, "ntotal != nx*ny*nz", 2, ntotal, nx*ny*nz); // check
	
	char line2[8]; in_t0 >> line2;
	if(strcmp(line2, "t0.ltcp") !=0) error(1, "the second line of t0.ltcp != t0.ltcp"); // check

	for(int a=0; a<ntotal; a ++){
		int state_in, i, j, k;
		if(in_t0.eof()) error(1, "reach end of file before finish reading all data");
		in_t0 >> state_in >> i >> j >> k;
		if(state_in > MAX_TYPE) error(1, "(in_t0) state input is larger than MAX_TYPE", 1, state_in); // check
		states[i][j][k]= state_in;

		if(state_in==Ttype) ntt ++;
	}
	cout << "ntt= " << ntt << endl;
	in_t0.close();
	cout << "t0.ltcp file reading completed; ntt= " << ntt << endl;

	// identify clusters and give them ids at t0
	int vcheck= 0;
	for(int a=0; a<ntotal; a ++){
		int x= (int) ((a)/nz)/ny;
		int y= (int) ((a)/nz)%ny;
		int z= (int)  (a)%nz;

		if(0==states[x][y][z]){
			vcheck ++;
			v0[0]= x; vpos[0]= x;
			v0[1]= y; vpos[1]= y;
			v0[2]= z; vpos[2]= z;
		}

		update_cid(x, y, z);
	}
	if(vcheck !=1) error(1, "(read_t0) vacancy larger than 1, need recoding");
	
	sum_csize();
}

void read_his_cal(){ // Reading history.info and calculating /////////
	ifstream in_hi(name_in_his, ios::in);
	if(! in_hi.is_open()) error(1, "in_t0 was not open");
	cout << "history.info is now reading to calculate csize..." << endl;

	vector <int> i_will; //the indexs list which will update cid later
	while(! in_hi.eof()){
		for(int a=1; a<=2; a++){
			int state_in, i_ltcp;
			in_hi >> state_in >> i_ltcp;
			if(0==state_in) update_vpos(i_ltcp);

			int i= (int) (i_ltcp/nz)/ny;
			int j= (int) (i_ltcp/nz)%ny;
			int k= (int)  i_ltcp%nz;
			states[i][j][k]= state_in;
			
			i_will.push_back(i_ltcp);
		}
		
		char c_dt[4];
		double dt= -1;
		in_hi >> c_dt >> dt;
		if(dt==-1) break;
		if(0 != strcmp(c_dt, "dt:"))  error(1, "(reading history) the format is incorrect");
	
		timestep ++;
		realtime += dt;

		for(int a=0; a<i_will.size(); a++){
			int x1= (int) (i_will.at(a)/nz)/ny;
			int y1= (int) (i_will.at(a)/nz)%ny;
			int z1= (int)  i_will.at(a)%nz;

			for(int b=0; b<n1nbr; b ++){
				int x2= pbc(x1+v1nbr[b][0], nx);
				int y2= pbc(y1+v1nbr[b][1], ny);
				int z2= pbc(z1+v1nbr[b][2], nz);
				int index= x2*ny*nz + y2*nz + z2;
				
				update_cid(x2, y2, z2);
			}
		}

		// properties calculations
		if(0==timestep%sample_cltr) sum_csize();
		if(0==timestep%sample_msd)  cal_msd();
		if(0==timestep%sample_sro)  cal_sro();
		if(0==timestep%sample_lro)  cal_lro();
		if(0==timestep%sample_lce)  cal_lce();
		// properties calculations
		
		i_will.clear();
		
		if(0==timestep%100000) cout << "T: " << timestep << " " << realtime << endl;
	}
}

void update_vpos(int iltcp){
	if(0==states[iltcp]) error(1, "(update_vpos) the ltcp where the vacancy jumps into is with a state 0");

	int x= (int) (iltcp/nz)/ny;
	int y= (int) (iltcp/nz)%ny;
	int z= (int)  iltcp%nz;
	
	if((x-vpos[0])==(nx-1)) vi[0] --; else if((x-vpos[0])==(-nx+1)) vi[0] ++; else if(abs(x-vpos[0])>1) error(1, "(update_vpos) the vacancy moves larger than 1");
	if((y-vpos[1])==(ny-1)) vi[1] --; else if((y-vpos[1])==(-ny+1)) vi[1] ++; else if(abs(y-vpos[1])>1) error(1, "(update_vpos) the vacancy moves larger than 1");
	if((z-vpos[2])==(nz-1)) vi[2] --; else if((z-vpos[2])==(-nz+1)) vi[2] ++; else if(abs(z-vpos[2])>1) error(1, "(update_vpos) the vacancy moves larger than 1");
	
	vpos[0]= x;
	vpos[1]= y;
	vpos[2]= z;
}

void update_cid(int x, int y, int z){ 
	int i_ltcp= x*ny*nz + y*nz + z;
	
	vector <int> i_ing;  //	the indexs which are being counted

	if(states[x][y][z]==Ttype){
		cid ++;	
		
		i_ing.push_back(i_ltcp);
		id_cltr[i_ltcp]= cid;

		for(int a=0; a<i_ing.size(); a++){ // !! be careful here: i_ing.size() increases during the iteration and automatically changes the condition
			int x1= (int) (i_ing.at(a)/nz)/ny;
			int y1= (int) (i_ing.at(a)/nz)%ny;
			int z1= (int)  i_ing.at(a)%nz;

			for(int b=0; b<n1nbr; b ++){
				int x2= pbc(x1+v1nbr[b][0], nx);
				int y2= pbc(y1+v1nbr[b][1], ny);
				int z2= pbc(z1+v1nbr[b][2], nz);
				int index= x2*ny*nz + y2*nz + z2;
	
				if(states[x2][y2][z2]==Ttype){
					bool IsNew= true;
					for(int c=0; c<i_ing.size(); c++)
						if((index) == i_ing.at(c)) IsNew=false;

					if(IsNew){
						i_ing.push_back(index);
						id_cltr[index]= cid;
				
					}
				}
			}
		}
	}
	else id_cltr[i_ltcp]= 0;
}

void sum_csize(){
	int Ncsize[5]= {0}; // # of clusters for single, l>1, l>5, l>10, l>30
	int Ncltr= 0;       // # of clusters, not atoms in clusters
	int sumNincltr= 0;  // total number of Ttype in clusters 

	int N_in_cltr[cid+1]; // for a certain cluster id, the number of counts
	for(int i=0; i<cid+1; i ++) N_in_cltr[i]= 0;

	int N_check1= 0;
	for(int i=0; i<nx*ny*nz; i ++){
		if(id_cltr[i] != 0){
			N_in_cltr[id_cltr[i]] ++;
			N_check1 ++;
		}
	}

	// calculated cluster related properties
	write_csind(cid, N_in_cltr);
	// calculated cluster related properties

	int N_check2= 0;
	for(int j=1; j<=cid; j ++){
		if(0==N_in_cltr[j]) continue;

		if(N_in_cltr[j] >= DEF_CLTR){ 
			Ncltr ++;
			sumNincltr += N_in_cltr[j];
		}
		
		N_check2 += N_in_cltr[j];

		if(N_in_cltr[j]>30) Ncsize[4] ++;
		if(N_in_cltr[j]>10) Ncsize[3] ++;
		if(N_in_cltr[j]> 5) Ncsize[2] ++;
		if(N_in_cltr[j]> 1) Ncsize[1] ++;
		if(N_in_cltr[j]==1) Ncsize[0] ++;
	}

	if(N_check1 != ntt) error(2, "number inconsistent for check1", 2, N_check1, ntt);
	if(N_check2 != ntt) error(2, "number inconsistent for check2", 2, N_check2, ntt);

	double Nave;
	if(0==Ncltr) Nave= 0;
	else 	     Nave= (double) sumNincltr / Ncltr; 
	
	fprintf(out_csize, "%d %e %d %d %d %d %d\n", timestep, realtime, Ncsize[0], Ncsize[1], Ncsize[2], Ncsize[3], Ncsize[4]);
	fprintf(out_csave, "%d %e %d %f\n", timestep, realtime, Ncltr, Nave);
}
	
void write_csind(int cid, int *N_in_cltr){
	if(cycle_out_csind%sample_cltr !=0) error(1, "cycle_out_csind must be a multiplier of sample_cltr", 2, cycle_out_csind, sample_cltr);
	
	if(0==timestep%cycle_out_csind){
		vector <int> size_cltr; //the indexs list which will update cid later
		vector <int> amount_cltr; //the amounts for the cluster sizes 

		for(int j=1; j<=cid; j ++){
			if(0==N_in_cltr[j]) continue;
		
			for(int a=0; a<size_cltr.size(); a ++){
				if(N_in_cltr[j]==size_cltr.at(a)){ 
					amount_cltr.at(a);
					goto endloop;
				}
			}
			
			size_cltr.push_back(N_in_cltr[j]);
			amount_cltr.push_back(1);
endloop:;
		}

		for(int a=0; a<size_cltr.size(); a ++){
			fprintf(out_csind, "%d %e %d %d\n", timestep, realtime, size_cltr.at(a), amount_cltr.at(a));
		}
	}
}

void cal_msd(){
	int v1= vpos[0] + vi[0]*nx;
	int v2= vpos[1] + vi[1]*ny;
	int v3= vpos[2] + vi[2]*nz;

	double dx= (v1-v0[0])*vbra[0][0] + (v2-v0[1])*vbra[1][0] + (v3-v0[2])*vbra[2][0];
	double dy= (v1-v0[0])*vbra[0][1] + (v2-v0[1])*vbra[1][1] + (v3-v0[2])*vbra[2][1];
	double dz= (v1-v0[0])*vbra[0][2] + (v2-v0[1])*vbra[1][2] + (v3-v0[2])*vbra[2][2];

	double msd= dx*dx + dy*dy + dz*dz;
	fprintf(out_msd, "%d %e %f\n", timestep, realtime, msd);
}

void cal_lce(){
	static int count_lce= 0; // count how many samples
	static int count_vcc= 0; // count how many vacancies
	static int ttsum= 0;     // sum of how many targeted types ltcp
	static double rtsum= 0;  // sum of realtime

	count_lce ++;
	rtsum += realtime;
	
	for(int i= 0; i<nx; i++){
		for(int j= 0; j<ny; j++){
			for(int k= 0; k<nz; k++){
				if(0==states[i][j][k]){
					count_vcc ++;

					for(int b=0; b<n1nbr; b ++){
						int x2= pbc(i+v1nbr[b][0], nx);
						int y2= pbc(j+v1nbr[b][1], ny);
						int z2= pbc(k+v1nbr[b][2], nz);
					
						if(Ttype==states[x2][y2][z2]) ttsum ++;
					}
				}
	}}}

	if(0==timestep%cycle_lce){
		if(0==count_lce) error(1, "LCE has no sample, cant divided by 0\n");

		double tsout= timestep - 0.5 * cycle_lce;
		double rtout= rtsum/count_lce;
		double lce= (double) ttsum / count_vcc / (double) n1nbr;
		fprintf(out_lce, "%.0f %e %f\n", tsout, rtout, lce);
	
		count_lce= 0;
		count_vcc= 0;
		ttsum= 0;
		rtsum= 0;
	}
}

void cal_sro(){
	int sum_i=  0; // sum of  i atoms 
	int sum_j=  0; // sum of  j atoms 
	int sum_ij= 0; // sum of ij bonds

	for(int i=0; i<nx; i++){ 
		for(int j=0; j<ny; j++){ 
			for(int k=0; k<nz; k++){ 
				if(jtype_sro == states[i][j][k]) 
					sum_j ++;
			
				if(itype_sro == states[i][j][k]){
					sum_i ++;

					for(int b=0; b<n1nbr; b ++){
						int x2= pbc(i+v1nbr[b][0], nx);
						int y2= pbc(j+v1nbr[b][1], ny);
						int z2= pbc(k+v1nbr[b][2], nz);
		
						if(states[x2][y2][z2]== jtype_sro) sum_ij ++;
					}	
				}
	}}}
	double cj=  (double)sum_j/nx/(double)ny/nz; 	// stoichiometric composition
	double zij= (double)sum_ij/sum_i;		// ij bond number for one atom

	double sro= 1 - zij/n1nbr/cj;
	fprintf(out_sro, "%d %e %f\n", timestep, realtime, sro);
}

void cal_lro(){
	int sum_lro= 0;

	for(int i=0; i<nx; i++){ 
		for(int j=0; j<ny; j++){ 
			for(int k=0; k<nz; k++){ 
				int sign_ltc;
				if(0==((i+j+k)%2)) sign_ltc=  1;
				else		   sign_ltc= -1;

				sum_lro += sign_ltc * states[i][j][k];
	}}}
	double lro= (double)abs(sum_lro) /nx /ny /nz;

	fprintf(out_lro, "%d %e %f\n", timestep, realtime, lro);
}
