// 2D Crank-Nicolson
// UNITS: x-> cm, y->um

#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <time.h> 
#include "Random64.h"
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

const double tolerance = 1e-9;
const double yMin = -10, yMax = 10; // um
const double tMin = 0, tMax = 2; // ns

const int ny = 2000;
const int nt = 500;

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
//////////// Classes /////////////
//////////////////////////////////
//////////////////////////////////

class Electron {
	double posY[nt];

public:
	void setPos(double y, int ti);
	double getPos(int ti);
};

void Electron::setPos(double y, int ti){
	posY[ti] = y;
};

double Electron::getPos(int ti){
	return posY[ti];
};

//////////////////////////////////
//////////////////////////////////
//////////// Program /////////////
//////////////////////////////////
//////////////////////////////////

int main(void){

	// Creating electrons
	Electron elecMatrix[ne];

	Crandom RanP(time(NULL));
  
  	/*for(int i=0; i<10;i++){
    	double num = (RanP.r())+0;
    	cout<<num<<"\n";
  	}*/

    int ei=0;
    int yi=ny/2;
    while(ei<ne/2){
    	bool done=0;
    	while(yi<=ny & ~done){
    		double y = yMin+yi*hy;
    		double num = (RanP.r());
    		//cout<<y<<"\t"<<psi0(y)<<"\n";
    		if(num<psi0(y)){
    			//cout<<y<<"\n";
    			elecMatrix[ei].setPos(y,0);
    			elecMatrix[ei+ne/2].setPos(-y,0);
    			ei++;
    			done=true;
    			//cout<<"inside\t"<<ei<<"\n";
    		}
    		yi++;
    		if(y>1){
    			yi=ny/2;
    		}
    	}
    	//cout<<ei<<"\n";
    }
    	/*for(int yi=ny/2;yi<=ny;yi++){
    		
    		//cout<<num<<"\n";
    		
    	}*/
    //}


	// initial conditions
	for(int yi=0; yi<=ny;yi++){
		double y = (double) yMin + yi*hy;
		matrix[0][yi]=psi0(y);
	}

	// Crank-Nicolson 
	cout<<"Executing Crank-Nicolson for |psiA+psiB|^2. "<<ny<<"x"<<nt<<"...\n";

	// matrix2 = matrix
	equalMatrix(matrix,matrix2);
	for(int ti=1;ti<=nt;ti++){
		done = false;

		//double t = tMin + ti*ht;
		// show progress bar
		float percentage = (float) ti/nt;
		progress(percentage);
		
		// relaxation of the system	
		while(!done){
			done = true;

			// matrix = matrix2
			equalMatrix(matrix2,matrix);

			for(int yi=0;yi<=ny;yi++){
				double y = yMin + yi*hy;
													
				matrix2[ti][yi] = (dy/(1.+(2.*dy)))*(matrix[ti][yi+1]+matrix[ti][yi-1]+matrix[ti-1][yi+1]+matrix[ti-1][yi-1])
									+((1.-2.*dy)/(1.+2.*dy))*matrix[ti-1][yi];
				
				/*if( (t-0<tolerance) && ((y>0.59) || (y<-0.59) || ( (y>-0.41)&&(y<0.41) ) ) ){
					matrix2[ti][yi] = 0;
					matrix[ti][yi] = 0;
				}*/


				if( abs(matrix[ti][yi]-matrix2[ti][yi]) > tolerance ){
					done = done & false;
				}
				else
					done = done & true;
			}			
		}		
	}
	std::cout<<"\n";
	// END OF Crank-Nicolson for |psiA+psiB|

	/*for(int ei=0;ei<ne/2;ei++){
		double initPos = (0.5-0.15)+(2*0.3/ne)*ei; //0.09, 0.18 // 0.2,0.4
		elecMatrix[ei].setPos(initPos,0);
		elecMatrix[ei+ne/2].setPos(-initPos,0);
	}*/
	

	// Electron trajectories
	cout<<"Calculating electron trajectories...\n";

	//cout<<findIndex(elecMatrix[0].getPos(0),ny,yMin,yMax)<<"\n";

	//cout<<elecMatrix[1].getPos(0)<<"\n";

	for(int ei=0;ei<ne;ei++){
		float percentage = (float) (ei+1)/ne;
		progress(percentage);
		//double initPos = 0;

		for(int ti=1;ti<nt;ti++){
			double oldY = elecMatrix[ei].getPos(ti-1);
			int yIndex = findIndex(oldY,ny,yMin,yMax);

			//double newY = oldY;
			//cout<<yIndex<<"\n";
			double newY = real(getP(matrix2,ti,yIndex))*ht+elecMatrix[ei].getPos(ti-1);

			//cout<<oldY<<"\t"<<newY<<"\n";
			elecMatrix[ei].setPos(newY,ti);
			//cout<<ei<<"\t"<<ti<<"\t"<<elecMatrix[1].getPos(0)<<"\t"<<newY<<"\n";
		}
	}

	// END of electron trajectories

	// Export Data
	cout<<"\nExporting data to file traj.dat...\n";
	ofstream trajFile ("trajsquare.dat");
	for(int ei=0;ei<ne;ei++){
		float percentage = (float) (ei+1)/ne;
		progress(percentage);
		for(int ti=0;ti<=nt;ti++){
			//matrixFile << ti <<"\n";

			double t = (double) tMin + ti*ht;
			trajFile << t << "\t" << vx*t << "\t" << elecMatrix[ei].getPos(ti) << "\n";
		}
		trajFile << "\n";
	}
	trajFile.close();
	std::cout<<"\n";

	// Export Data
	cout<<"Exporting data to file psi2.dat...\n";
	ofstream psiFile ("psi2square.dat");
	for(int ti=0;ti<=nt;ti++){
		float percentage = (float) ti/nt;
		progress(percentage);
		//matrixFile << ti <<"\n";
		for (int yi=0; yi<=ny; yi++){
			double t = (double) tMin + ti*ht;
			double y = (double) yMin + yi*hy;
			psiFile << t << "\t" << vx*t << "\t" << y << "\t" << (double) pow(abs(matrix2[ti][yi]),2) << "\n";
		}
		psiFile << "\n";
	}
	psiFile.close();
	std::cout<<"\n";
	// END OF Export Data

	return 0;
}

//////////////////////////////////
//////////////////////////////////
////// Function definitions //////
//////////////////////////////////
//////////////////////////////////

double psi0(double y){
	double sigmay = 0.09; //um
	double mu = 0.5; // um
	//return 0.5*gaussian(y,-mu,sigmay)+0.5*gaussian(y,mu,sigmay);
	return (fabs(y-mu)<sigmay)? 1 : (fabs(y+mu)<sigmay)? 1 : 0; // barrier
	//return 10*cos(y);	
}

double gaussian(double x, double mu, double sigma){
	return exp(-0.5*pow((x-mu)/sigma,2))/sqrt(sigma*sqrt(2*pi));
}

void equalMatrix(complex <double> matrix[nt+1][ny+1], complex <double> matrix2[nt+1][ny+1]){
	for(int ti=0; ti<=nt; ti++){
		for(int yi=0; yi<=ny;yi++){
			matrix2[ti][yi]=matrix[ti][yi];
		}
	}
}

double getRho(complex <double> psi){
	return pow(abs(psi),2);
}

complex <double> getP(complex <double> matrix[nt+1][ny+1],int time, int pos){
	complex <double> current;

	current = (conj(matrix[time][pos])*getGradient(matrix,time,pos,hy)) - (matrix[time][pos]*conj(getGradient(matrix,time,pos,hy)));
	current = -current*hbar_m*i/2.;

	return current/getRho(matrix[time][pos]);	
}

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