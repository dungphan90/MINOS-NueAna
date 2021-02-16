//automatically generate code for MC/Data comparison
#include "Riostream.h"
#include "string"
#include <vector>

void gen_md()
{
  vector<string> varName;
  vector<Float_t> beg;
  vector<Float_t> end;
  vector<Int_t> nbins;
  vector<string> gtype;
  string dumName, dumtype;
  Float_t dumstart, dumend;
  Int_t dumbins, nvar;
  ifstream ins;
  ins.open("AllParam.txt");
  
  nvar = 0;
  
  //read in the file
  while(!ins.eof()) {
    ins>>dumName>>dumstart>>dumend>>dumbins>>dumtype;
    if(!ins.eof()){
      varName.push_back(dumName);
      beg.push_back(dumstart);  end.push_back(dumend);
      nbins.push_back(dumbins);   gtype.push_back(dumtype);
      nvar++;
    }
  }
  cout<<nvar<<" variables read in"<<endl;
  ins.close();
  ofstream out;
  out.open("plot_md.C");
  out<<"#include \"TH1F.h\""<<endl;
  out<<"#include \"TCanvas.h\""<<endl;
  out<<endl;
  out<<"void plot_md()"<<endl;
  out<<"{"<<endl;
  out<<"  TH1F *htemp;"<<endl;
  out<<"  TFile file(\"testbig.root\");"<<endl;
  out<<endl;
  out<<"  TCanvas *c1 = new TCanvas(\"c1\",\"c1\",600,400);"<<endl;
  out<<endl;
  for(int i = 0; i<int(varName.size()); i++){
    out<<"  htemp = (TH1F*) file.Get(\"allcomp/"<<varName[i]<<"_d\");"<<endl;;
    out<<"  htemp->SetTitle(\""<<varName[i]<<"\");"<<endl;
    out<<"  htemp->Draw(\"e\");"<<endl;
    out<<"  htemp = (TH1F*) file.Get(\"allcomp/"<<varName[i]<<"_m_all\");"<<endl;
    out<<"  htemp->SetLineColor(2);"<<endl;
    out<<"  htemp->SetTitle(0);"<<endl;
    out<<"  htemp->Draw(\"same\");"<<endl;
    out<<"  c1->Print(\""<<varName[i]<<".gif\");"<<endl;
    out<<endl;
  }
  out<<endl;
  out<<"}"<<endl;
  out.close();

}
