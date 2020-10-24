#include<cstdlib>
#include <iostream>
using namespace std;
template <class ntype>

class Particle{
  double v[3]={};
public:
  void set(double *val){
    for(int i =0;i<3;i++){
      v[i]=val[i];
    }
  }

  void show(){
    cout <<"Vector position"<<endl;
    for(int i =0;i<3;i++){
      cout <<v[i]<<'\n';
    }
  }
};
/*
class Molecule: public Particle{
  double u[3][3];

public:
  void
};
*/

int main(){
  Particle <double> p;
  double x[3]={1.,1.,1.};
  p.set(x);
  p.show();
}
