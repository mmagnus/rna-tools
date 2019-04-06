/***********************************************************
  Compute the error function, adapted from Numeric Receipt.
  --fding
***/
#include <cmath>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
using namespace std;
static const int ITMAX=100;
static const double EPS=3.0e-7;

double gammln(double xx){
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;
  
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}

void gser(double* gamser, double a, double x, double* gln){
  int n;
  double sum,del,ap;
  
  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) cerr << "x less than 0 in routine GSER"<< endl;
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ap += 1.0;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    cerr << "a too large, ITMAX too small in routine GSER" << endl;
    return;
  }
}

void gcf(double* gammcf, double a, double x, double *gln){
  int n;
  double gold=0.0,g,fac=1.0,b1=1.0;
  double b0=0.0,anf,ana,an,a1,a0=1.0;
  
  *gln=gammln(a);
  a1=x;
  for (n=1;n<=ITMAX;n++) {
    an=static_cast<double>(n);
    ana=an-a;
    a0=(a1+a0*ana)*fac;
    b0=(b1+b0*ana)*fac;
    anf=an*fac;
    a1=x*a0+anf*a1;
    b1=x*b0+anf*b1;
    if (a1) {
      fac=1.0/a1;
      g=b1*fac;
      if (fabs((g-gold)/g) < EPS) {
	*gammcf=exp(-x+a*log(x)-(*gln))*g;
	return;
      }
      gold=g;
    }
  }
  cerr << "a too large, ITMAX too small in routine GCF" << endl;
}

double gammp(double a, double x){
  double gamser,gammcf,gln;
  
  if (x < 0.0 || a <= 0.0) cerr << "Invalid arguments in routine GAMMP" << endl;
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return 1.0-gammcf;
  }  
}

double erf(double x){
  return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);  
}

int main(int argc, char* argv[]){
  if(argc<3){
    cout << "usage: q-value.linux rnaLength predictedRMSD " << endl;
    exit(1);
  }
  int n = atol(argv[1]);
  double rmsd = atof(argv[2]);
  //cout << n << " " << rmsd << endl;
  /*<rmsd> = a*N^0.41-b*/
  double a_wo = 4.6*sqrt(2.0), b_wo = 9.1*sqrt(2.0);
  double a_w  = 3.6*sqrt(2.0), b_w  = 11.2*sqrt(2.0);
  //double std = 1.6*sqrt(2.0);
  double std = 1.8;
  
  double ave_wo = a_wo*pow(n,0.41)-b_wo;
  double ave_w  = a_w*pow(n,0.41)-b_w;
  
  double z_wo = (rmsd - ave_wo)/std;
  double z_w  = (rmsd - ave_w)/std;
  double sqrt2 = sqrt(2.0);
  
  double q_wo = (1.0+erf(z_wo/sqrt2))/2.0;
  double q_w = (1.0+erf(z_w/sqrt2))/2.0;
  
  cout << "Without base pair constraints" << endl;
  char buf[1000];
  sprintf(buf, " The average RMSD by chance: %.1f ", ave_wo); cout << buf << endl;
  sprintf(buf, " Z-score of the prediction: %.2f", z_wo); cout << buf << endl;
  if(q_wo<1.0e-6){
  q_wo = 1.0e-06;
  sprintf(buf, " p-value of the prediction: < %.2le", q_wo); cout << buf << endl;
  }
  else{
  sprintf(buf, " p-value of the prediction: %.2le", q_wo); cout << buf << endl;
  }

  cout << endl; 
  cout << "With base pair constraints" << endl;
  sprintf(buf, " The average RMSD by chance: %.1f ", ave_w); cout << buf << endl;
  sprintf(buf, " Z-score of the prediction: %.2f", z_w); cout << buf << endl;
  if(q_w<1.0e-6){
  q_w = 1.0e-6;
  sprintf(buf, " p-value of the prediction: < %.2le", q_w); cout << buf << endl;
  }
  else{
  sprintf(buf, " p-value of the prediction: %.2le", q_w); cout << buf << endl;
  }
  
}
