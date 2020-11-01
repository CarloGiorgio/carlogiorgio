# include <iostream>
# include <math.h>
# include <cstdlib>
#include <time.h>
# define pi atan(1)*4
using namespace std;

class MyVec{

  double *v;

public:
  MyVec(){
    v=new double[3];
    for(int i=0;i<3;i++){
      v[i]=0.;
    }
  }
  MyVec(double x,double y,double z){
    v=new double[3];
    v[0]=x;
    v[1]=y;
    v[3]=z;
  }

  void init(double *x){
    for(int i=0;i<3;i++){
      v[i]=x[i];
    }
  }
  void rand_1(){
    for(int i=0;i<3;i++){
      this->v[i]=1-2*drand48();
    }
  }

  void copy_vec(double* x){
    for(int i=0;i<3;i++){
      x[i]=this->v[i];
    }
  }
//moltiplication with a scalar
  friend MyVec operator *(double x,MyVec&para){
    MyVec temp;
    for(int i=0;i<3;i++){
      temp.v[i]=para.v[i]*x;
    }
    return (temp);
  }
//commutative of moltiplication
  friend MyVec operator *(MyVec& para,double x){
    return (x*para);
  }

  friend MyVec operator /(MyVec&para,double x){
    MyVec temp;
    for(int i=0;i<3;i++){
      temp.v[i]=para.v[i]/x;
    }
    return (temp);
  }

//Addition and subtraction
  friend MyVec operator +(MyVec & para,MyVec& parat){
    MyVec temp;
    for(int i=0;i<3;i++){
      temp.v[i]=para.v[i]+parat.v[i];
    }
    return temp;
  }
  friend MyVec operator -(MyVec & para,MyVec& parat){
    MyVec temp;
    for(int i=0;i<3;i++){
      temp.v[i]=para.v[i]-parat.v[i];
    }
    return temp;
  }
//inner product
  friend double operator ^ ( MyVec& para,MyVec &temp){
    double x=0.;
    for(int i=0;i<3;i++){
      x+=temp.v[i]*para.v[i];
    }
    return x;
  }
  //cross product
  friend MyVec operator &&(MyVec &para,MyVec &parat){
    MyVec temp;
    for(int i=0;i<3;i++){
      temp.v[i]=para.v[(i+1)%3]*parat.v[(i+2)%3]-para.v[(i+2)%3]*parat.v[(i+1)%3];
    }
    return temp;
  }
  void printV(){
    for(int i=0;i<3;i++){
      cout<<this->v[i]<<" ";
    }
    cout<<endl;
  }
};


class Particle: public MyVec{
  MyVec p;
public:
  Particle(){}

  Particle(double x,double y,double z){
    double v[3];
    v[0]=x;
    v[1]=y;
    v[2]=z;
    p.init(v);
  }
};


class Molecule: public Particle{
public:
  Particle p;
  MyVec o;
  Molecule(){

    double r1,r2;
    double o_r[3]={0.};
    r1=2*pi*drand48();
    r2=pi*drand48();
    o_r[0]=cos(r1)*sin(r2);
    o_r[1]=sin(r1)*sin(r2);
    o_r[2]=cos(r2);
    o.init(o_r);
  }

  Molecule(double x,double y,double z){
    double v[3];
    v[0]=x;
    v[1]=y;
    v[2]=z;
    o.init(v);
  }
};

void G_S_method(){
  MyVec u1,u2,u3,temp,temp1;
  u1.rand_1();
  u2.rand_1();
  u1=u1/(sqrt(u1^u1));
  temp=(u1^u2)*u1;
  u2=u2-temp;
  u2=u2/sqrt(u2^u2);
  u3=u1&&u2;
  u3=u3/sqrt(u3^u3);
  u1.printV();
  u2.printV();
  u3.printV();
  cout<<endl<<(u1^u1)<<endl<<(u2^u2)<<endl<<(u2^u1)<<endl<<(u3^u1)<<endl;
}

int main(){
  srand48(time(NULL));
  double x;
  MyVec u,u1,u3,u4,u5;
  Molecule m;
  m.o.printV();
  x=m.o^m.o;
  cout<<x<<endl;
  m.p.printV();
  cout<<endl;
  u.rand_1();
  u1.rand_1();
  u3=u1-u;
  u4=u1+u;
  u5=u*4;
  u.printV();
  u1.printV();
  u3.printV();
  u4.printV();
  u5.printV();
  cout<<endl;
  G_S_method();
  return 0;
}
