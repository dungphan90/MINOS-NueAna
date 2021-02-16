#include <iostream>
#include <fstream>
#include <vector>
#include "TRandom.h"

using namespace std;
vector<vector<float> > RecursiveVFill(vector<vector<float> > v, vector<float> r);

void MakeRandFile(int nrows)
{

//eventually, I want random numbers from a ncol dim. Gauss.
//for now just generate ncol numbers from ncol independent gauss
   const int ncol=24;
   float mean[ncol]={1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,2.5e-3,1.0};
   float rms[ncol]={0.3,0.3,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.3,0.3,
		    0.3,0.3,0.3,0.3,
		    0.3,0.3,0.2,0.2};

   ofstream out("randnumbers.txt");

   for(int i=0;i<ncol;i++){
      if(i==0){
	 out<<mean[i];
      }
      else{
	 out<<" "<<mean[i];
      }
   }
   out<<endl;

   TRandom *r=new TRandom();
   for(int i=0;i<nrows;i++){
      for(int j=0;j<ncol;j++){
	 float rnum = abs(r->Gaus(mean[j],rms[j]*mean[j]));
	 if(j==0){
	    out<<rnum;
	 }
	 else{	    
	    out<<" "<<rnum;
	 }
      }
      out<<endl;
   }

   out.close();

}


void MakeUniformFile()
{

   const int ncol=24;
   float mean[ncol]={1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,2.5e-3,1.0};
   float rms[ncol]={0.3,0.3,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.3,0.3,
		    0.3,0.3,0.3,0.3,
		    0.3,0.3,0.2,0.2};

   ofstream out("uniformnumbers.txt");


   for(int nvar=0;nvar<ncol;nvar++){       
     if((nvar>=2&&nvar<14)||(nvar>14)){
       continue;
     }
     for(int nrows=0;nrows<11;nrows++){
       float r=mean[nvar]+(5-nrows)*0.1;
       for(int col=0;col<ncol;col++){
	 float o=mean[col];
	 if(col==nvar){
	   o=r;
	 }
	 if(col==0){
	   out<<o;
	 }
	 else{
	   out<<" "<<o;
	 }
       }
       out<<endl;
     }
   }

   out.close();

}


void MakeBigUniformFile()
{

   const int ncol=25;

/*

//FOR PRODUCTION REWEIGHTING
   float mean[ncol]={1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,2.175e-3,0.925,0.0};
   float rms[ncol]={0.03,0.03,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.04,0.04,
		    0.04,0.04,0.04,0.04,
		    0.04,0.04,0.07,0.07,.005};

   int NPTS[ncol]={3,3,0,0,
		   0,0,0,0,
		   0,0,0,0,
		   0,0,3,0,
		   0,0,0,0,
		   0,0,3,3,21};
*/

/*
//special for chris

   float mean[ncol]={1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,0.0,0.0,0.0};
   float rms[ncol]={0.011,0.011,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.015,0.015,
		    0.015,0.015,0.015,0.015,
		    0.015,0.015,0.07,0.07,.005};

   int NPTS[ncol]={17,17,0,0,
		   0,0,0,0,
		   0,0,0,0,
		   0,0,17,0,
		   0,0,0,0,
		   0,0,0,0,0};
*/

//Another special for chris

   float mean[ncol]={1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,1.,1.,
		     1.,1.,2.175e-3,0.925,0.0};
   float rms[ncol]={0.0,0.0,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.0,0.0,
		    0.0,0.0,0.07,0.01,.005};

   int NPTS[ncol]={0,0,0,0,
		   0,0,0,0,
		   0,0,0,0,
		   0,0,0,0,
		   0,0,0,0,
		   0,0,17,17,17};

//   ofstream out("biguniformnumbers-specialforchris.txt");
   ofstream out("biguniformnumbers-specialforchrisFD.txt");
   vector< vector<float> >vals;
   for(int i=ncol-1;i>=0;i--){
     //     cout<<"on col "<<i<<endl;
     //     cout<<"doing "<<NPTS[i]<<" loops "<<endl;
     vector<float> rvec;
     for(int j=0;j<NPTS[i];j++){       
       float r=mean[i]+((int)(NPTS[i]/2.)-j)*rms[i]*mean[i];
       if(i==ncol-1){
	 r=mean[i]+j*rms[i];
       }
       rvec.push_back(r);
       //       cout<<" added "<<r<<" to j="<<j<<endl;
     }
     if(rvec.size()==0){
       rvec.push_back(mean[i]);
     }
     vals=RecursiveVFill(vals,rvec);
   }

   cout<<"Number of rows "<<vals.size()<<endl;

   for(unsigned int i=0;i<vals.size();i++){
     //     cout<<"number of cols "<<vals[i].size()<<" in row"<<i<<endl;
     for(unsigned int j=0;j<vals[i].size();j++){
       float o=vals[i][j];
       if(j==0){
	 out<<o;
       }
       else{
	 out<<" "<<o;
       }
     }
     out<<endl;
   }

   out.close();

}

vector<vector<float> > RecursiveVFill(vector<vector<float> > v, vector<float> r)
{
  vector< vector<float> > nv;
  for(unsigned int i=0;i<r.size();i++){
    for(unsigned int j=0;j<v.size();j++){
      vector<float> nr;
      nr.push_back(r[i]);
      for(unsigned int k=0;k<v[j].size();k++){
	nr.push_back(v[j][k]);
      }
      nv.push_back(nr);
    }
    if(v.size()==0){
      vector<float> nrz;
      nrz.push_back(r[i]);
      nv.push_back(nrz);
    }
  }
  return nv;
}
