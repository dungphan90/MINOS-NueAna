#include "TYHadronWeight.h"
#include "TMath.h"

double TYDefPar[5] = {1.20448,0.0066884,0.589213,0.0191205,1};

double TYHadronWeight(double w2, double totpt, int npi0, int res, double *par)
{
  //Params:
  //par[0] =  offset for A
  //par[1] =  slope for A
  //par[2] =  offset for B
  //par[3] =  slope for B
  //par[4] =  scale factor for npi0
 
  double wei = 1.0;
  if ( res != 1003) return wei;

  double A = 0.2289*TMath::Power(w2+18.14,-1.17)*(1+32.51*w2); //mean
  double B = 0.1886*TMath::Power(w2,0.3501);   //sigma

  double Ap = (par[0]+w2*par[1])*A;
  double Bp = (par[2]+w2*par[3])*B;

  double ExpTerm = -TMath::Power((totpt-Ap)/Bp,2)+TMath::Power((totpt-A)/B,2); 
  wei =wei*B/Bp*TMath::Exp(ExpTerm/2.0);

  double C = 0.01951+0.5909*TMath::Log(w2);
  double Cp = (par[4])*C;
  wei = wei*TMath::Exp(-(Cp-C))*TMath::Power(Cp/C,npi0);

  if (wei<0) wei = 1;
  if (wei<1e-5) wei = 0;
  if (wei>10) wei = 10;

  return wei;
}

double TYHadronWeight(NueRecord *nr, double *par)
{
   double w2 = double(nr->mctrue.w2);

   return TYHadronWeight(w2, nr->shi.totpt, nr->shi.npi0, 
                    nr->mctrue.resonanceCode, par);
}

double TYHadronWeight(NueRecord *nr)
{
  return TYHadronWeight(nr, TYDefPar);
}

double TYHadronWeight(double w2, double totpt, int npi0, int res)
{
  return TYHadronWeight(w2,totpt,npi0,res, TYDefPar);
}

