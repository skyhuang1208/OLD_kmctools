#include <iostream>
#include <fstream>

using namespace std;

// parameters //
const double vbra[3][3]= {{-0.5,  0.5,  0.5}, { 0.5, -0.5,  0.5}, { 0.5,  0.5, -0.5}};
const int n1nbr= 8;
const int v1nbr[8][3]= {{ 1,  0,  0}, { 0,  1,  0}, { 0,  0,  1},
			{-1,  0,  0}, { 0, -1,  0}, { 0,  0, -1},
			{ 1,  1,  1}, {-1, -1, -1}	       };

ofstream of_nn("./ucell.xyz"); 
ofstream of_nnnn("./nnnn.xyz"); 


void out_nn(int i, int j, int k){
	double x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0];
	double y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1];
	double z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2];
	
	of_nn << "1 " << x << " " << y << " " << z << endl;
}

void out_nnnn(int i, int j, int k){
	double x= i*vbra[0][0] + j*vbra[1][0] + k*vbra[2][0];
	double y= i*vbra[0][1] + j*vbra[1][1] + k*vbra[2][1];
	double z= i*vbra[0][2] + j*vbra[1][2] + k*vbra[2][2];
	
	of_nnnn << "1 " << x << " " << y << " " << z << endl;
}

int main(){
	cout << "Calculation is starting..." << endl;

	of_nn << "9\n\n";
	of_nnnn << "64\n\n";
	
	out_nn(5, 5, 5);
	for(int a=0; a<n1nbr; a ++){
		int i= 5+v1nbr[a][0];
		int j= 5+v1nbr[a][1];
		int k= 5+v1nbr[a][2];
		out_nn(i, j, k);
	
		for(int b=0; b<n1nbr; b ++){
			int i2= i+v1nbr[b][0];
			int j2= j+v1nbr[b][1];
			int k2= k+v1nbr[b][2];
			out_nnnn(i2, j2, k2);
		}
	}
}
