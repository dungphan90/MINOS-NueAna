double CountPOT(string files);

void MakeGrid(Int_t whichPID, int whichSep = 0) 
{
  gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs();
  gSystem->Load("libNueAnaExtrapolation.so");

  MCReweight *mcr = &MCReweight::Instance();
  NeugenWeightCalculator *n = new NeugenWeightCalculator();
  mcr->AddWeightCalculator(n);

  string nearMCFiles = "./ManyFiles/AnaNue-n130*SysTrim*.root";
  string farMCFiles = "./ManyFiles/AnaNue-f21*SysTrim*.root";

  Double_t nearPOT = CountPOT(nearMCFiles.c_str());
  Double_t farPOT = CountPOT(farMCFiles.c_str())/3.0;

  TChain *nearChain = new TChain("NueMini");
  TChain *farChain = new TChain("NueMini");

  nearChain->Add(nearMCFiles.c_str());
  farChain->Add(farMCFiles.c_str());

  string nearCCMCFiles = "./ManyFiles/AnaNue-n130*-Trim-Mini.root";
  string farCCMCFiles = "./ManyFiles/AnaNue-f21*-Trim-Mini.root";
                                                                                
  Double_t nearCCPOT = CountPOT(nearCCMCFiles.c_str());
  Double_t farCCPOT = CountPOT(farCCMCFiles.c_str());
                                                                                
  TChain *nearCCChain = new TChain("NueMini");
  TChain *farCCChain = new TChain("NueMini");
                                                                                
  nearCCChain->Add(nearCCMCFiles.c_str());
  farCCChain->Add(farCCMCFiles.c_str());

  NueExtrapolationJB extrap;

  extrap.AddChain(farChain, farPOT);
  extrap.AddChain(nearChain, nearPOT);
  extrap.AddChain(nearCCChain, nearCCPOT, false);
  extrap.AddChain(farCCChain, farCCPOT, false);
  extrap.UseMCAsData(false);
 
  if(whichSep == 0)  extrap.SetSeparation("HOO");
  if(whichSep == 1)  extrap.SetSeparation("MRCC");

   string file;
     if( whichSep == 1){
         file = "DataHists/MRCC/SPECTRA_GAMMA.root";
     }
 
     if( whichSep == 0){
       file = "DataHists/HornOnOff/HornOnOff_NDspecta_TotError_15PercentBnue.root";
     }
 
     extrap.SetDataFileName(file);

  extrap.SetNueRecoBins(100,0.,100.);
  extrap.SetTrueBins(1200,0.,120.0);
  extrap.SetCCRecoBins(100, 0., 100.);
  extrap.SetExtrapMethod(3);

  Selection::Selection_t sel = whichPID;

  //add a standard entry:  (numu CC osc pars, no nue oscillations; SKZP)
  NueSystematic *nueSys = new NueSystematic("Standard");
  nueSys->AddSystematic(Systematic::kOscProb,2);
  nueSys->SetOscParams(0.554,0.7854,0.3,8.2e-5,2.4e-3,0,1);
  nueSys->AddSystematic(Systematic::kSKZP,0);
  extrap.AddNueSystematic(nueSys);

  int nDelta23 = 1; //25;
  int nDelta = 101; //200;
  int nUe3 = 101; //500;

  double d23start = 2.4e-3;
  double d23end = 2.5e-3;
  double dStart = 0.00;
  double dEnd = 2.0;
  double ue3start = 0.0;
  double ue3end = 0.4;

  double d23Step = 1;
  if (nDelta23 > 1) d23Step = (d23end-d23start)/(nDelta23-1);
  double dStep   = (dEnd-dStart)/(nDelta-1);
  double ue3Step = (ue3end-ue3start)/(nUe3-1);
  double pi = 3.1415926; 
 
  for(int i= 0; i< nDelta23; i++){
   for(int j = 0; j < nDelta; j++){
      for(int k = 0; k < nUe3; k++){

        double d23 = d23start + i*d23Step;
        double d   = (dStart + j*dStep)*pi;
        double ue3 = ue3start + k*ue3Step;

        double Th13 = TMath::ASin(TMath::Sqrt(ue3))/2.;
        double m2 = TMath::Sin(Th13)*TMath::Sin(Th13);

        char sysname[256]; sprintf(sysname,"AlterOsc_%i_%i_%i_N",i,j,k);
        nueSys = new NueSystematic(sysname);
        nueSys->AddSystematic(Systematic::kOscProb,2);
        nueSys->AddSystematic(Systematic::kSKZP,0);
        nueSys->SetOscParams(0.554,0.7854,Th13,8.2e-5,d23,d,1);
        extrap.AddNueSystematic(nueSys);

        char sysname[256]; sprintf(sysname,"AlterOsc_%i_%i_%i_I",i,j,k);
        nueSys = new NueSystematic(sysname);
        nueSys->AddSystematic(Systematic::kOscProb,2);
        nueSys->AddSystematic(Systematic::kSKZP,0);
        nueSys->SetOscParams(0.554,0.7854,Th13,8.2e-5,d23,d,-1);
        extrap.AddNueSystematic(nueSys);
      }
    }
  }

  string outfiletag = "ExtrapFile_" + string(Selection::AsString(sel));

  if(whichSep == 0)  outfiletag += "_HOO";
  if(whichSep == 1)  outfiletag += "_MRCC";
 
  outfiletag += ".root";

  extrap.SetOutputFile(outfiletag);
  extrap.RunExtrapStudy(sel);
  //write out the file:
}

double CountPOT(string file)
{
  TChain *nChain = new TChain("pottree");
                                                                                
  nChain->Add(file.c_str());
                                                                                
  double pot;
  nChain->SetMakeClass(1);
  nChain->SetBranchAddress("pot", &pot);
                                                                                
  Double_t POT = 0.0;
                                                                                
  for(int i=0; i < nChain->GetEntries(); i++)
  {
      nChain->GetEntry(i);
      POT += pot;
  }

  return POT*1.0e12;
}
