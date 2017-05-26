#include<iostream>
#include<cmath>
using namespace std;

//Constantes del generador aleatorio
class Crandom{
  unsigned long long u,v,w;
public:
  Crandom(unsigned long long j);
  unsigned long long int64();
  double r() {return 5.42101086242752217E-20 * int64();}
  unsigned int int32(){return (unsigned int) int64();};
  double exponencial(float tau);
  double gauss(float mu,float sigma);
};
Crandom::Crandom(unsigned long long j){
    v=4101842887655102017LL; w=1;
    u = j ^ v; int64();
    v = u; int64();
    w = v; int64();
  }
unsigned long long Crandom::int64() {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    return (x + v) ^ w;
  }

double Crandom::exponencial(float tau){
  return -tau*log(r());
}
double Crandom:: gauss(float mu,float sigma){
  return sigma*sqrt(-2*log(r()))*cos(2*M_PI*r())+mu;
}


