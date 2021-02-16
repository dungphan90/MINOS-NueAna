#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TH1.h"
#include <iostream>

using namespace std;

void Test(double bl=1290.0);
void DrawMultiple(double bl=1290.0);
double ElecAppear_Cervera(double *,double *);
double MuDisappear(double *,double *);
double MuSurvive(double *,double *);
double ElecAppear_Winter(double *,double *);

void Test(double bl){

  TCanvas *can = new TCanvas("TestCan","TestCan");
  can->SetLogx();

  TF1 *form = new TF1("form",ElecAppear_Cervera,0.05,100,10);
  form->SetParameter(0,bl); //baseline (km)
  form->SetParameter(1,1.0); //sinsq_2th23
  form->SetParameter(2,1.0); //sinsq_2th12
  form->SetParameter(3,0.04); //sinsq_2th13
  form->SetParameter(4,0.0026); //dmsq23
  form->SetParameter(5,5.0e-5); //dmsq12
  form->SetParameter(6,2.65); //density
  form->SetParameter(7,0); //cp phase
  form->SetParameter(8,1); //anti-nu
  form->SetParameter(9,0); //print messages?
  form->SetNpx(1000);
  form->Draw();

}

void DrawMultiple(double bl){

  TCanvas *can = new TCanvas("MultiCan","MultiCan");
  can->cd();

  //bl,th23,th12,th13,dm23,dm12,density,delta,nu/nubar,messages
  
  double sinsq_vals[3] = {0.86,1.0,0.04};
  double dmsq_vals[2] = {7.3e-5,2.e-3};
  double density = 2.65;

  can->SetLogx();

  TF1 *form1 = new TF1("form1",ElecAppear_Cervera,0.05,6,10);
  form1->SetParameters(bl, sinsq_vals[1], sinsq_vals[0], sinsq_vals[2], 
		       dmsq_vals[1], dmsq_vals[0], density,
		       0., 1, 0);
  form1->SetNpx(1000);
  form1->SetLineColor(1);
  form1->SetLineStyle(1);

  TF1 *form2 = new TF1("form2",ElecAppear_Cervera,0.05,6,10);
  form2->SetParameters(bl, sinsq_vals[1], sinsq_vals[0], sinsq_vals[2], 
		       dmsq_vals[1], dmsq_vals[0], density,
		       TMath::Pi()/4., 1, 0);
  form2->SetNpx(1000);
  form2->SetLineColor(3);
  form2->Draw("");
  form1->Draw("same");

  TF1 *form3 = new TF1("form3",ElecAppear_Cervera,0.05,6,10);
  form3->SetParameters(bl, sinsq_vals[1], sinsq_vals[0], sinsq_vals[2], 
		       dmsq_vals[1], dmsq_vals[0], density,
		       -TMath::Pi()/4., 1, 0);
  form3->SetNpx(1000);
  form3->SetLineColor(2);
  form3->SetLineStyle(2);
  //form3->Draw("same");

  TF1 *form4 = new TF1("form4",ElecAppear_Cervera,0.05,6,10);
  form4->SetParameters(bl, sinsq_vals[1], sinsq_vals[0], sinsq_vals[2], 
		       dmsq_vals[1], dmsq_vals[0], density,
		       TMath::Pi()/4., -1, 0);
  form4->SetNpx(1000);
  form4->SetLineColor(4);
  form4->SetLineStyle(1);
  form4->Draw("same");

  TF1 *form5 = new TF1("form5",ElecAppear_Cervera,0.05,6,10);
  form5->SetParameters(bl, sinsq_vals[1], sinsq_vals[0], sinsq_vals[2], 
		       dmsq_vals[1], dmsq_vals[0], density,
		       3*TMath::Pi()/4., 1, 0);
  form5->SetNpx(1000);
  form5->SetLineColor(2);
  form5->SetLineStyle(1);
  form5->Draw("same");

  TLegend *leg = new TLegend(0.7,0.45,0.87,0.63,"");
  leg->AddEntry(form1,"NO, #delta_{cp}=0^{o}","l");
  leg->AddEntry(form2,"NO, #delta_{cp}=45^{o}","l");
  //leg->AddEntry(form3,"NO, #delta_{cp}=-45^{o}","l");
  leg->AddEntry(form5,"NO, #delta_{cp}=135^{o}","l");
  leg->AddEntry(form4,"RO, #delta_{cp}=45^{o}","l");
  leg->SetBorderSize(0);
  leg->Draw();

  form2->SetMaximum(0.2);
  form2->GetHistogram()->SetXTitle("#nu Energy (GeV)");
  form2->GetHistogram()->SetYTitle("P(#nu_{#mu} #rightarrow #nu_{e})");
  form2->GetHistogram()->SetTitle("#nu_{e} Appearance Probability");

  char bltext[256];
  sprintf(bltext,"Baseline=%.0f km",bl);
  TLatex *tex0 = new TLatex(0.6,0.8,bltext);
  tex0->SetNDC();
  tex0->SetTextSize(0.033);
  tex0->Draw();

  char sinsqtext[256];
  sprintf(sinsqtext,"sin^{2}2#theta_{ij}(12,23,13)=%.2f, %.2f, %.2f",
	  sinsq_vals[0],sinsq_vals[1],sinsq_vals[2]);
  TLatex *tex1 = new TLatex(0.6,0.75,sinsqtext);
  tex1->SetNDC();
  tex1->SetTextSize(0.033);
  tex1->Draw();

  char dmsqtext[256];
  sprintf(dmsqtext,"dm^{2}_{ij}(21,32)=%.1e, %.1e eV^{2}",
	  dmsq_vals[0],dmsq_vals[1]);
  TLatex *tex2 = new TLatex(0.6,0.7,dmsqtext);
  tex2->SetNDC();
  tex2->SetTextSize(0.033);
  tex2->Draw();

  char mdtext[256];
  sprintf(mdtext,"Matter density=%.2f g/cm^{3}",density);
  TLatex *tex3 = new TLatex(0.6,0.65,mdtext);
  tex3->SetNDC();
  tex3->SetTextSize(0.033);
  tex3->Draw();

}

double ElecAppear_Cervera(double *x,double *par)
{

  //x[0] = E
  //Params:
  //par[0] = L
  
  //par[1] = sin^2_2th23
  //par[2] = sin^2_2th12
  //par[3] = sin^2_2th13
  
  //par[4] = dm23^2
  //par[5] = dm12^2

  //par[6] = density

  //par[7] = d_cp

  //par[8] = +/-;

  double E = x[0]; //energy
  double L = par[0]; //baseline
  double plusminus = int(par[8]);

  //standard rock:
  double density = par[6]; //g/cm^{3}
  double z_a = 0.5; //average Z/A

  double A_av = 6.02214199e23; //avogadro's number
  double invCmToeV = 1.97e-5; //convert 1/cm to eV
  double invCmToGeV = 1.97e-14; //convert 1/cm to GeV
  double invKmToeV = 1.97e-10; //convert 1/km to eV

  double Gf = 1.166391e-5; //fermi constant (GeV^{-2})
  double ne = z_a*A_av*density; //electron density #/cm^{3}
  double ne_natunits = ne*invCmToeV*invCmToGeV*invCmToGeV;
  //electron density with units Gev^{2} eV
  //Gev^{2} to cancel with GeV^{-2} in Gf

  double sinsq_2th23 = par[1]; //atm
  double sinsq_2th12 = par[2]; //sol
  double sinsq_2th13 = par[3]; //chooz
  
  double cossq_2th23 = 1.-sinsq_2th23;
  double cossq_2th13 = 1.-sinsq_2th13;
  double cossq_2th12 = 1.-sinsq_2th12;

  double cos_th23 = TMath::Cos(TMath::ACos(sqrt(cossq_2th23))/2.);
  double sin_th23 = sqrt(1. - cos_th23*cos_th23);

  double cos_th13 = TMath::Cos(TMath::ACos(sqrt(cossq_2th13))/2.);
  double sin_th13 = sqrt(1. - cos_th13*cos_th13);

  double cos_th12 = TMath::Cos(TMath::ACos(sqrt(cossq_2th12))/2.);
  double sin_th12 = sqrt(1. - cos_th12*cos_th12);

  double d_cp = par[7];
  
  double dmsq_23 = par[4];
  double dmsq_12 = par[5];
  double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}

  double Delta23 = dmsq_23/(2.*E*1e9); //eV
  double Delta12 = dmsq_12/(2.*E*1e9);
  double Delta13 = dmsq_13/(2.*E*1e9);
  
  double A = sqrt(2.)*Gf*ne_natunits; //eV
  double B = fabs(A - plusminus*Delta13); //eV
  double J = cos_th13*sqrt(sinsq_2th12)*sqrt(sinsq_2th13)*sqrt(sinsq_2th23);

  double p1 = sin_th23*sin_th23*sinsq_2th13*TMath::Power(Delta13/B,2)
    *TMath::Power(TMath::Sin(B*L/(invKmToeV*2.)),2);
  
  double p2 = 0;
  double p3 = 0;
  if(density!=0){
    p2 = cos_th23*cos_th23*sinsq_2th12*TMath::Power(Delta12/A,2)
      *TMath::Power(TMath::Sin(A*L/(invKmToeV*2.)),2);
    
    p3 = J*Delta12*Delta13*TMath::Sin(A*L/(invKmToeV*2.))
      *TMath::Sin(B*L/(invKmToeV*2.))
      *TMath::Cos(plusminus*d_cp - Delta13*L/(invKmToeV*2.))/(A*B);
  }

  if(par[9]==1){

    cout << "Energy = " << E << endl;
    cout << "Baseline = " << L << endl;
    cout << "hierarchy = " << plusminus << endl;
    cout << "SinSq(2theta23) = " << sinsq_2th23 << endl;
    cout << "Sin(theta23) = " << sin_th23 << endl;
    cout << "SinSq(2theta12) = " << sinsq_2th12 << endl;
    cout << "Sin(theta12) = " << sin_th12 << endl;
    cout << "SinSq(2theta13) = " << sinsq_2th13 << endl;
    cout << "Sin(theta13) = " << sin_th13 << endl;
    
    cout << "CP phase = " << d_cp << endl;

    cout << "dmsq_23 = " << dmsq_23 << endl;
    cout << "dmsq_12 = " << dmsq_12 << endl;
    cout << "dmsq_13 = " << dmsq_13 << endl;

    cout << "Delta23 = " << Delta23 << endl;
    cout << "Delta13 = " << Delta13 << endl;
    cout << "Delta12 = " << Delta12 << endl;

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "J = " << J << endl;

    cout << "p1 = " << p1 << endl;
    cout << "p2 = " << p2 << endl;
    cout << "p3 = " << p3 << endl;
    
    cout << "Probability = " << p1+p2+p3 << endl;

  }
  if(p1+p2+p3>1) return 1;
  return p1+p2+p3;

}

double ElecAppear_Winter(double *x,double *par)
{

  //x[0] = E
  //Params:
  //par[0] = L
  
  //par[1] = sin^2_2th23
  //par[2] = sin^2_2th12
  //par[3] = sin^2_2th13
  
  //par[4] = dm23^2
  //par[5] = dm12^2

  //par[6] = density

  //par[7] = d_cp

  //par[8] = +/-;

  double E = x[0]; //energy
  double L = par[0]; //baseline
  double plusminus = int(par[8]);

  //standard rock:
  double density = par[6]; //g/cm^{3}
  double z_a = 0.5; //average Z/A

  double A_av = 6.02214199e23; //avogadro's number
  double invCmToeV = 2.0e-5; //convert 1/cm to eV
  double invCmToGeV = 2.0e-14; //convert 1/cm to GeV
  double invKmToeV = 1.97e-10; //convert 1/km to eV

  double Gf = 1.166391e-5; //fermi constant (GeV^{-2})
  double ne = z_a*A_av*density; //electron density #/cm^{3}
  double ne_natunits = ne*invCmToeV*invCmToGeV*invCmToGeV;
  //electron density with units Gev^{2} eV
  //Gev^{2} to cancel with GeV^{-2} in Gf 

  double sinsq_2th23 = par[1]; //atm
  double sinsq_2th12 = par[2]; //sol
  double sinsq_2th13 = par[3]; //chooz
  
  double cossq_2th23 = 1.-sinsq_2th23;
  double cossq_2th13 = 1.-sinsq_2th13;
  double cossq_2th12 = 1.-sinsq_2th12;

  double cos_th23 = sqrt((sqrt(cossq_2th23)+1.)/2.);
  double sin_th23 = sqrt(1. - cos_th23*cos_th23);

  double cos_th13 = sqrt((sqrt(cossq_2th13)+1.)/2.);
  double sin_th13 = sqrt(1. - cos_th13*cos_th13);

  double cos_th12 = sqrt((sqrt(cossq_2th12)+1.)/2.);
  double sin_th12 = sqrt(1. - cos_th12*cos_th12);

  double d_cp = par[7];
  double sin_dcp = TMath::Sin(d_cp);
  double cos_dcp = TMath::Cos(d_cp);
  
  double dmsq_23 = par[4];
  double dmsq_12 = par[5];
  double dmsq_13 = dmsq_23+dmsq_12; //ev^{2}
  
  double alpha = dmsq_12/dmsq_13; //dimensionless
  double Delta = dmsq_13*L/(invKmToeV*4.*E*1e9); //dimensionless
  double eta = cos_th13*sqrt(sinsq_2th12)*sqrt(sinsq_2th23);
  double A = plusminus*2.*sqrt(2.)*Gf*ne_natunits*E*1e9/dmsq_13;//dimensionless
  
  double p1 = sinsq_2th13*sin_th23*sin_th23
    *TMath::Power(TMath::Sin((A-1)*Delta)/(A-1),2);
  
  double p2 = 0;
  double p3 = 0;  
  double p4 = 0;
  if(density!=0) {
    
    p2 = alpha*sqrt(sinsq_2th13)*eta
      *sin_dcp*TMath::Sin(Delta)*TMath::Sin(A*Delta)
      *TMath::Sin((A-1)*Delta)/(A*(A-1));
    
    p3 = alpha*sqrt(sinsq_2th13)*eta
      *cos_dcp*TMath::Cos(Delta)*TMath::Sin(A*Delta)
      *TMath::Sin((A-1)*Delta)/(A*(A-1));
    
    p4 = alpha*alpha*cos_th23*cos_th23*sinsq_2th12
      *TMath::Sin(A*Delta)*TMath::Sin(A*Delta)/(A*A);
  }

  if(par[9]==1){

    cout << "Energy = " << E << endl;
    cout << "Baseline = " << L << endl;
    cout << "hierarchy = " << plusminus << endl;
    cout << "SinSq(2theta23) = " << sinsq_2th23 << endl;
    cout << "Sin(theta23) = " << sin_th23 << endl;
    cout << "SinSq(2theta12) = " << sinsq_2th12 << endl;
    cout << "Sin(theta12) = " << sin_th12 << endl;
    cout << "SinSq(2theta13) = " << sinsq_2th13 << endl;
    cout << "Sin(theta13) = " << sin_th13 << endl;
    
    cout << "CP phase = " << d_cp << endl;
    cout << "Sin(d_cp) = " << sin_dcp << endl;
    cout << "Cos(d_cp) = " << cos_dcp << endl;

    cout << "dmsq_23 = " << dmsq_23 << endl;
    cout << "dmsq_12 = " << dmsq_12 << endl;
    cout << "dmsq_13 = " << dmsq_13 << endl;

    cout << "alpha = " << alpha << endl;
    cout << "Delta = " << Delta << endl;
    cout << "eta = " << eta << endl;
    cout << "A = " << A << endl;

    cout << "p1 = " << p1 << endl;
    cout << "p2 = " << p2 << endl;
    cout << "p3 = " << p3 << endl;
    cout << "p4 = " << p4 << endl;
    
    cout << "Probability = " << p1+plusminus*p2+p3+p4 << endl;

  }

  return p1+plusminus*p2+p3+p4;

}

double MuDisappear(double *x,double *par)
{

  //x[0] = E
  //Params:
  //par[0] = L
  
  //par[1] = sin^2_2th23
  //par[2] = sin^2_2th12
  //par[3] = sin^2_2th13
  
  //par[4] = dm23^2
  //par[5] = dm12^2

  //par[6] = density

  //par[7] = d_cp

  //par[8] = +/-;

  double E = x[0]; //energy
  double L = par[0]; //baseline
  double plusminus = int(par[8]);

  //standard rock:
  double density = par[6]; //g/cm^{3}
  double z_a = 0.5; //average Z/A

  double A_av = 6.02214199e23; //avogadro's number
  double invCmToeV = 1.97e-5; //convert 1/cm to eV
  double invCmToGeV = 1.97e-14; //convert 1/cm to GeV
  double invKmToeV = 1.97e-10; //convert 1/km to eV

  double Gf = 1.166391e-5; //fermi constant (GeV^{-2})
  double ne = z_a*A_av*density; //electron density #/cm^{3}
  double ne_natunits = ne*invCmToeV*invCmToGeV*invCmToGeV;
  //electron density with units Gev^{2} eV
  //Gev^{2} to cancel with GeV^{-2} in Gf

  double sinsq_2th23 = par[1]; //atm
  double sinsq_2th12 = par[2]; //sol
  double sinsq_2th13 = par[3]; //chooz
  
  double cossq_2th23 = 1.-sinsq_2th23;
  double cossq_2th13 = 1.-sinsq_2th13;
  double cossq_2th12 = 1.-sinsq_2th12;

  double cos_th23 = TMath::Cos(TMath::ACos(sqrt(cossq_2th23))/2.);
  double sin_th23 = sqrt(1. - cos_th23*cos_th23);

  double cos_th13 = TMath::Cos(TMath::ACos(sqrt(cossq_2th13))/2.);
  double sin_th13 = sqrt(1. - cos_th13*cos_th13);

  double cos_th12 = TMath::Cos(TMath::ACos(sqrt(cossq_2th12))/2.);
  double sin_th12 = sqrt(1. - cos_th12*cos_th12);

  double d_cp = par[7];
  
  double dmsq_23 = par[4];
  double dmsq_12 = par[5];
  double dmsq_13 = dmsq_23+dmsq_12; //eV^{2}

  double Delta23 = dmsq_23/(2.*E*1e9); //eV
  double Delta12 = dmsq_12/(2.*E*1e9);
  double Delta13 = dmsq_13/(2.*E*1e9);
  
  double A = sqrt(2.)*Gf*ne_natunits; //eV
  double B = fabs(A - plusminus*Delta13); //eV
  double J = cos_th13*sqrt(sinsq_2th12)*sqrt(sinsq_2th13)*sqrt(sinsq_2th23);

  //  double p1 = sinsq_2th23*TMath::Power(cos_th13,4)
  //*TMath::Power(TMath::Sin(Delta13*L/(invKmToeV*2.)),2);  //numu->nutau

  double p1 = sinsq_2th23
    *TMath::Power(TMath::Sin(Delta23*L/(invKmToeV*2.)),2); //1 - numu survival

  if(par[9]==1){

    cout << "Energy = " << E << endl;
    cout << "Baseline = " << L << endl;
    cout << "hierarchy = " << plusminus << endl;
    cout << "SinSq(2theta23) = " << sinsq_2th23 << endl;
    cout << "Sin(theta23) = " << sin_th23 << endl;
    cout << "SinSq(2theta12) = " << sinsq_2th12 << endl;
    cout << "Sin(theta12) = " << sin_th12 << endl;
    cout << "SinSq(2theta13) = " << sinsq_2th13 << endl;
    cout << "Sin(theta13) = " << sin_th13 << endl;
    
    cout << "CP phase = " << d_cp << endl;

    cout << "dmsq_23 = " << dmsq_23 << endl;
    cout << "dmsq_12 = " << dmsq_12 << endl;
    cout << "dmsq_13 = " << dmsq_13 << endl;

    cout << "Delta23 = " << Delta23 << endl;
    cout << "Delta13 = " << Delta13 << endl;
    cout << "Delta12 = " << Delta12 << endl;

    cout << "A = " << A << endl;
    cout << "B = " << B << endl;
    cout << "J = " << J << endl;

    cout << "p1 = " << p1 << endl;
    
    cout << "Probability = " << p1 << endl;

  }

  return p1;

}


double MuSurvive(double *x,double *par)
{
  return 1. - MuDisappear(x,par); 
}
