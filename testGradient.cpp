// 2D Crank-Nicolson
// UNITS: x-> cm, y->um

#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
using namespace std;

//const double hbar = 1.05*pow(10,-34);
//const double m = 9.11*pow(10,-31);
//const double hbar_m = 0.1157675; // um2/ns;
const double hbar_m = 0.115772423; // um2/ns;
const complex<double> i(0.0,1.0);
const double pi = 3.1415926536;
const double vx = 18; // cm/ns

const complex <double> D = i*hbar_m/(2.); 
//const complex <double> D = 0.5;

const double tolerance = 1e-13;
const double yMin = -10, yMax = 10; // um
const double tMin = 0, tMax = 2; // ns

const int ny = 700;
const int nt = 700;

const double hy = (double)(yMax-yMin)/ny;
const double ht = (double)(tMax-tMin)/nt;

const complex<double> dy = D*ht/(2*hy*hy);

const int ne = 50;

bool done = false;

complex <double> matrix[nt+1][ny+1];
complex <double> matrix2[nt+1][ny+1];

//double matrixTraj[nt][ne]; // x=vt, we calculate y using the guiding eq.
//double matrixY[nt][ne];

//////////////////////////////////
//////////////////////////////////
////// Function declarations /////
//////////////////////////////////
//////////////////////////////////

double psi0(double y);
double gaussian(double x, double mu, double sigma);
void equalMatrix(complex <double> matrix[nt+1][ny+1], complex <double> matrix2[nt+1][ny+1]);
void progress(float percentage);
double getRho(complex <double> psi);
complex <double> getGradient(complex <double> matrix[nt+1][ny+1],int time, int pos,double step);
complex <double> getP(complex <double> matrix[nt+1][ny+1],int time, int pos);
int findIndex(double num, int nNum, double min, double max);

//////////////////////////////////
//////////////////////////////////
//////////// Program /////////////
//////////////////////////////////
//////////////////////////////////

int main(void){

	// initial conditions
	for(int yi=0;yi<=ny;yi++){
		double y = (double) yMin + yi*hy;
		matrix[0][yi]=sin(y);
	}

	for(int yi=1;yi<ny;yi++){
		matrix2[0][yi]=getGradient(matrix,0,yi,hy);
	}
	
	// Export Data
	cout<<"\nExporting data to file testGrad.dat...\n";
	ofstream testFile ("testGrad.dat");
	for(int yi=0;yi<ny;yi++){
		float percentage = (float) (yi+1)/ny;
		progress(percentage);
		double y = (double) yMin + yi*hy;
		testFile << y << "\t" << real(matrix[0][yi]) << "\t" << real(matrix2[0][yi]) << "\n";
	}
	testFile.close();
	std::cout<<"\n";

	return 0;
}


//////////////////////////////////
//////////////////////////////////
////// Function definitions //////
//////////////////////////////////
//////////////////////////////////

/*complex <double> getP(complex <double> matrix[nt+1][ny+1],int time, int pos){
	complex <double> current;

	current = conj(matrix[time][pos])*getGradient(matrix,time,pos,hy) - matrix[time][pos]*conj(getGradient(matrix,time,pos,hy));
	current = -current*hbar_m*i/2.;

	return current/getRho(matrix[time][pos]);	
}*/

complex <double> getGradient(complex <double> matrix[nt+1][ny+1],int time, int pos,double step){
	return (matrix[time][pos+1]-matrix[time][pos-1])/(2.*step);
}

void progress(float percentage){
	int barWidth = 70;
    std::cout << "[";
    int pos = barWidth * percentage;
    //std::cout << pos;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(100*percentage) << "%\r";
    std::cout.flush();
}

int findIndex(double num, int nNum, double min, double max){
	double hnum = (double) (max-min)/nNum;

	for(int i=0;i<nNum;i++){
		double numTemp = (double) min+i*hnum;

		//cout<<abs(num-numTemp)<<"\n";
		if(abs(num-numTemp)<=hnum/2){
			//cout<<num<<"\t"<<numTemp<<"\t"<<hnum/2<<"\n";
			return i;
		}
	}
	cout<<"Index not found! num="<<num<<", range ["<<min<<","<<max<<"], nNum="<<nNum<<"\n";
	return 0;
}