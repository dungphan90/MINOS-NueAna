#include <iostream>
#include "TMatrix.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "math.h"
#include "TLeaf.h"

using namespace std;

void Covar(string file, double pidcut=0)
{
  TFile *f;
 	f =  TFile::Open(file.c_str());
  
  TTree *tree = (TTree *)f->Get("TrainTree");

  TTree *descriptor = (TTree *)f->Get("descriptor");

	descriptor->GetEntry(0);
  const int NVAR=descriptor->GetLeaf("npars")->GetValue();
  printf("using %d vars\n",NVAR);

  double pars[100];
  double pid;
  tree->SetMakeClass(1);
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("pars",1);
  tree->SetBranchAddress("pars",&pars);
  tree->SetBranchStatus("pid",1);
  tree->SetBranchAddress("pid",&pid);

  
  int z=0; 
  double means[NVAR];
  for(int i=0;i<NVAR;i++){
    means[i]=0.;
    pars[i]=0.;
  }
  while(tree->GetEntry(z)){
    if(z%10000==0){
      cout<<"On entry "<<z<<endl;
    }
    z++;
    if(pid<pidcut) continue;

	for(int i=0;i<NVAR;i++){
    	means[i]+=1.* pars[i];
    }
  }

  double covar[NVAR][NVAR];
  double vars[NVAR];
  for(int i=0;i<NVAR;i++){
    means[i]/=(1.*z);
    cout<<"Means "<<i<<" "<<means[i]<<endl;
    for(int j=0;j<NVAR;j++){
      covar[i][j]=0.;
    }
  }

  z=0; 
  while(tree->GetEntry(z)){
    if(z%10000==0){
      cout<<"On entry "<<z<<endl;
    }
    z++;
    if(pid<pidcut) continue;

	for(int i=0;i<NVAR;i++){
    	vars[i]=1.* pars[i];
    }
 
  
    for(int i=0;i<NVAR;i++){
      for(int j=0;j<NVAR;j++){
		covar[i][j]+=(vars[i]-means[i])*(vars[j]-means[j]);
      }
    }
  }

  //  TFile *g = new TFile(Form("covar-%d-highpid.root",det),"RECREATE");
  TFile *g = new TFile("covar.root","RECREATE");
  TH2D *covarhist = new TH2D("covarhist","",NVAR+2,-0.5,NVAR+.5,NVAR+2,-0.5,NVAR+.5);
  TMatrixT<double> mat(NVAR,NVAR);
  for(int i=0;i<NVAR;i++){
    for(int j=0;j<NVAR;j++){
      covar[i][j]/=(1.*z-1);
      mat[i][j]=covar[i][j];
      covarhist->SetBinContent(i+1,j+1,covar[i][j]);
      if(i==j){
		cout<<"Covars "<<i<<" "<<j<<" "<<covar[i][j]<<endl;
      }
    }
  }

  TH2D *corrhist = new TH2D("corrhist","",NVAR+2,-0.5,NVAR+.5,NVAR+2,-0.5,NVAR+.5);
  for(int i=0;i<NVAR;i++){
    for(int j=0;j<NVAR;j++){
      corrhist->SetBinContent(i+1,j+1,fabs(covar[i][j]/(sqrt(covar[i][i])*sqrt(covar[j][j]))));
    }
  }
  mat.Write();
  g->Write();
//  g->Close();
}

