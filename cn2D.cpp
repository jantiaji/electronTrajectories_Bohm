// 2D Crank-Nicolson

#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
using namespace std;

const double hbar = 1;
const double m = 1;
const complex<double> i(0.0,1.0);
const double pi = 3.1415926535897;

const complex <double> D = i*hbar/(2.*m);
//const complex <double> D = 0.5;


const double tolerance = 1e-2;
const double xMin = -50, xMax = 20;
const double yMin = -10, yMax = 10;
const double tMin = 0, tMax = 100;

const int nx = 65;
const int ny = 50;
const int nt = 75;

const double hx = (double)(xMax-xMin)/nx;
const double hy = (double)(yMax-yMin)/ny;
const double ht = (double)(tMax-tMin)/nt;

const complex<double> dx = D*ht/(2*hx*hx);
const complex<double> dy = D*ht/(2*hy*hy);

bool done = false;

// Function Definitions
double psi0(double x, double y);
void equalMatrix(complex <double> matrix[nt+1][nx+1][ny+1], complex <double> matrix2[nt+1][nx+1][ny+1]);

int main(void){

	//std::cout << (i*i) <<"\n";
	//std::cout<<"test \n";
	complex <double> matrix[nt+1][nx+1][ny+1];
	complex <double> matrix2[nt+1][nx+1][ny+1];

	// initial conditions
	for(int xi=0; xi<=nx;xi++){
		for(int yi=0; yi<=ny;yi++){
			double x = (double) xMin + xi*hx;
			double y = (double) yMin + yi*hy;
			matrix[0][xi][yi]=psi0(x,y);
		}
	}

	// matrix2 = matrix
	equalMatrix(matrix,matrix2);

	// Crank-Nicolson
	for(int ti=1;ti<=nt;ti++){
		done = false;
		std::cout<<ti<<"\n";	
		// relax		
		while(!done){
			done = true;

			// matrix = matrix2
			equalMatrix(matrix2,matrix);
			
			for(int xi=0;xi<=nx;xi++){
				for(int yi=0;yi<=ny;yi++){

					double x = xMin + xi*hx;
					double y = yMin + yi*hy;
					
					complex <double> C = dx*(matrix[ti][xi+1][yi]+matrix[ti][xi-1][yi]+matrix[ti-1][xi+1][yi]+matrix[ti-1][xi-1][yi])+
					dy*(matrix[ti][xi][yi+1]+matrix[ti][xi][yi-1]+matrix[ti-1][xi][yi+1]+matrix[ti-1][xi][yi-1]);
					
					matrix2[ti][xi][yi] = ((C)-matrix[ti-1][xi][yi]*(2.*dx+2.*dy-1.))/(2.*dx+2.*dy+1.);
					
					// border conditions
					matrix2[ti][0][yi] = 0;//psi0(xMin,y);
					matrix2[ti][nx][yi] = 0; //psi0(xMax,y);
					matrix2[ti][xi][0] = 0; //psi0(x,yMin);
					matrix2[ti][xi][ny] = 0; //psi0(x,yMax);

					matrix[ti][0][yi] = 0;//psi0(xMin,y);
					matrix[ti][nx][yi] = 0; //psi0(xMax,y);
					matrix[ti][xi][0] = 0; //psi0(x,yMin);
					matrix[ti][xi][ny] = 0; //psi0(x,yMax);

					matrix2[ti][0][0] = psi0(xMin,yMin);
					matrix2[ti][nx][ny] = psi0(xMax,yMax);
					matrix2[ti][0][ny] = psi0(xMin,yMax);
					matrix2[ti][nx][0] = psi0(xMax,yMin);

					matrix[ti][0][0] = psi0(xMin,yMin);
					matrix[ti][nx][ny] = psi0(xMax,yMax);
					matrix[ti][0][ny] = psi0(xMin,yMax);
					matrix[ti][nx][0] = psi0(xMax,yMin);

					if( abs(matrix[ti][xi][yi]-matrix2[ti][xi][yi]) > tolerance ){
						done = done & false;
					}
					else
						done = done & true;
				}			
			}
		}
		
			
	}
	// END OF Crank-Nicolson
	
	// Export Data
	ofstream matrixFile ("psi2_2D.dat");
	for(int ti=0;ti<=nt;ti++){
		//matrixFile << ti <<"\n";
		for (int xi=0; xi<=nx; xi++){
			for (int yi=0; yi<=ny; yi++){
				double t = (double) tMin + ti*ht;
				double x = (double) xMin + xi*hx;
				double y = (double) yMin + yi*hy;
				matrixFile << t << "\t" << x << "\t" << y << "\t" << (double) pow(abs(matrix2[ti][xi][yi]),2) << "\n";
			}
			//matrixFile << "\n";
		}
		matrixFile << "\n";
	}
	matrixFile.close();

	return 0;
}

double psi0(double x, double y){
	double sigmax = 0.5;
	double sigmay = 2.0;
	return 2*exp(-pow(x/sigmax,2)-pow((y-4)/sigmay,2)) + 2*exp(-pow(x/sigmax,2)-pow((y+4)/sigmay,2));
	//return 10*cos(x);
	//return 4*exp(-pow(x+2,2)-pow(y+2,2))+2*exp(-pow(x,2)-pow(y,2))+exp(-pow(x-2,2)-pow(y-2,2));	
}

void equalMatrix(complex <double> matrix[nt+1][nx+1][ny+1], complex <double> matrix2[nt+1][nx+1][ny+1]){
	for(int ti=0; ti<=nt; ti++){
		for(int xi=0; xi<=nx;xi++){
			for(int yi=0; yi<=ny;yi++){
				matrix2[ti][xi][yi]=matrix[ti][xi][yi];
			}
		}
	}
}
