#include "Trimmer.h"
#include "NueAna/NueStandard.h"
#include "DataUtil/MCInfo.h"
#include "DataUtil/DataQualDB.h"
 
Trimmer::Trimmer(){ 
  outSet = false;
  fReweight = false;
  fBaseline = 735;
  fDeltaMSquare = 0.0027;
  fUe3Square = 0.025;
  fTheta23 = TMath::Pi()/4;

  fOverWritePOT = false;
  separatebyRunPeriod = false;
  isRHC = false;
}

void Trimmer::AddFiles(string infiles)
{
   files.push_back(infiles);
}

void Trimmer::SetOutputFile(string file)
{
   outSet = true;
   outFile = file;
}

void Trimmer::Config(Registry &r)
{
  fCuts.Config(r);
}

void Trimmer::SetDeltaMSquare(double dm2) { fDeltaMSquare = dm2; fReweight = true;}
void Trimmer::SetUe3Square(double dUe32) { fUe3Square = dUe32; fReweight = true;}
void Trimmer::SetTheta23(double t23) { fTheta23 = t23; fReweight = true;}
void Trimmer::SetBaseline(double bl) { fBaseline = bl; fReweight = true;}

void Trimmer::RunTrimmer()
{
  
  NueRecord *nr = new NueRecord();
  NuePOT *np = new NuePOT();

  TChain *chain = new TChain("ana_nue");
  chain->SetBranchAddress("NueRecord",&nr);
                                                                        
  TChain *pchain = new TChain("pottree");
  pchain->SetBranchAddress("NuePOT", &np);
 
  for(unsigned int i = 0; i < files.size(); i++)
  {
      chain->Add(files[i].c_str());
      pchain->Add(files[i].c_str());
  }
  
  if(!outSet){
     string file;
                                                                        
     if(pchain->GetEntries() > 0){
        pchain->GetEntry(0);
        file = pchain->GetFile()->GetName();
     }
     string minifile = file.substr(file.find_last_of("/")+1, file.find_last_of(".root")-file.find_last_of("/") - 5);
     minifile += "-Trim.root";
  
     outFile = minifile;
  }
  if(outFile == "-Trim.root"){
     cout<<"No input file found"<<endl; 
     return;
  }
  
  string out1,out2,out3,out11,out12,out13;
  if(separatebyRunPeriod)
  {
    string base = outFile.substr(0,outFile.size()-5);//strip off '.root' and add label for run period
    out1 = base + "_RunPeriod1.root";
    out2 = base + "_RunPeriod2.root";
    out3 = base + "_RunPeriod3.root";
    out11 = base + "_RunPeriod11.root";
    out12 = base + "_RunPeriod12.root";
    out13 = base + "_RunPeriod13.root";
    cout<<"Setting output to:"<<endl;
    cout<<out1<<endl;
    cout<<out2<<endl;
    cout<<out3<<endl;
    cout<<out11<<endl;
    cout<<out12<<endl;
    cout<<out13<<endl;


  }
  else
  {
    cout<<"Setting output to "<<outFile<<endl;
    out1 = outFile;
  }

  // And now the actual looping
  
  TFile *f1 = new TFile(out1.c_str(),"RECREATE");
  TFile *f2=0,*f3=0, *f11=0, *f12=0, *f13=0;
  if(separatebyRunPeriod)
  {
    f2 = new TFile(out2.c_str(),"RECREATE");
    f3 = new TFile(out3.c_str(),"RECREATE");

    f11 = new TFile(out11.c_str(),"RECREATE");
    f12 = new TFile(out12.c_str(),"RECREATE");
    f13 = new TFile(out13.c_str(),"RECREATE");

  }
  
  f1->cd();
  TTree *tree1 = new TTree("ana_nue","ana_nue");
  TTree::SetBranchStyle(1);
  TBranch* br1 = tree1->Branch("NueRecord","NueRecord", &nr );
  br1->SetAutoDelete(kFALSE);
  
  TTree *tree2=0;
  TTree *tree3=0;
  TBranch* br2=0;
  TBranch* br3=0;

  TTree *tree11=0;
  TTree *tree12=0;
  TTree *tree13=0;
  TBranch* br11=0;
  TBranch* br12=0;
  TBranch* br13=0;



  if(separatebyRunPeriod)
  {
    f2->cd();
    tree2 = new TTree("ana_nue","ana_nue");
    br2 = tree2->Branch("NueRecord","NueRecord", &nr );
    br2->SetAutoDelete(kFALSE);
    
    f3->cd();
    tree3 = new TTree("ana_nue","ana_nue");
    br3 = tree3->Branch("NueRecord","NueRecord", &nr );
    br3->SetAutoDelete(kFALSE);

    
    f11->cd();
    tree11 = new TTree("ana_nue","ana_nue");
    br11 = tree11->Branch("NueRecord","NueRecord", &nr );
    br11->SetAutoDelete(kFALSE);
    
    f12->cd();
    tree12 = new TTree("ana_nue","ana_nue");
    br12 = tree12->Branch("NueRecord","NueRecord", &nr );
    br12->SetAutoDelete(kFALSE);
    
    f13->cd();
    tree13 = new TTree("ana_nue","ana_nue");
    br13 = tree13->Branch("NueRecord","NueRecord", &nr );
    br13->SetAutoDelete(kFALSE);

  }
  
  bool isRunPeriod1=false;
  bool isRunPeriod2=false;
  bool isRunPeriod3=false;
  bool isRunPeriod11=false;
  bool isRunPeriod12=false;
  bool isRunPeriod13=false;
  
  Int_t n = chain->GetEntries();
  int count = 0;

  bool isMC = false;
  
  //for recalculating the pot tree for each run period
  double TotPOT1 = 0.0,TotPOT2 = 0.0,TotPOT3 = 0.0,TotPOT11 = 0.0, TotPOT12 = 0.0, TotPOT13 = 0.0;
  double TotPOT1_nocut = 0.0,TotPOT2_nocut = 0.0,TotPOT3_nocut = 0.0,TotPOT11_nocut = 0.0,TotPOT12_nocut = 0.0,TotPOT13_nocut = 0.0;
  Int_t nruns=0;
  Int_t nsnarls1=0,nsnarls2=0,nsnarls3=0;
  Int_t nsnarls11=0,nsnarls12=0,nsnarls13=0;
 
  int lastrun=-1;

  double potweight;
  double nddataweight;
  
  chain->GetEntry(0);
  ReleaseType::Release_t rel = nr->GetHeader().GetRelease();
  BeamType::BeamType_t beamtype = nr->GetHeader().GetBeamType();
  Detector::Detector_t det = nr->GetHeader().GetVldContext().GetDetector();
  
  for(int i=0;i<n;i++){
    if(i%10000==0) cout << 100*float(i)/float(n)
                       << "% done" << endl;
    chain->GetEvent(i);
    
    NueConvention::NueEnergyCorrection(nr);
    if(isRHC) NueConvention::RHCNueEnergyCorrection(nr);
    
    nr->anainfo.isRHC = isRHC;
    
    //NueStandard::ModifyANNPID(nr);
    NueStandard::SetE50PID(nr);    


    nr->eventq.passFarDetQuality = DataUtil::IsGoodData(nr->GetHeader().GetVldContext());
    nr->eventq.passNearDetQuality = nr->eventq.passFarDetQuality;
    
    if(i == 0){
      SimFlag::SimFlag_t sim = nr->GetHeader().GetVldContext().GetSimFlag();
      isMC = (sim == SimFlag::kMC);
    }
    
    if(separatebyRunPeriod && !isMC)
    {
      cout<<"'separatebyRunPeriod' option is not meant to be used for data - pottree would not be calculated correctly for data.  Qutting..."<<endl;
      return;
    }
    
    //sanity check
    if(NueStandard::PassesPOTStandards(nr))
    {
      if(nr->GetHeader().GetEventNo() == 0 || nr->GetHeader().GetEvents() == 0)
      {
        if(nr->bmon.bI < 0 && !isMC) cout<<"Unexpected behavior"<<endl;
      }
    }
   
    isRunPeriod1=false;
    isRunPeriod2=false;
    isRunPeriod3=false;
    isRunPeriod11=false;
    isRunPeriod12=false;
    isRunPeriod13=false;

    if(separatebyRunPeriod)
    {
      if(NueStandard::IsRun1(nr))
      {
        isRunPeriod1=true;
      }
      if(NueStandard::IsRun2(nr))
      {
        isRunPeriod2=true;
      }
      if(NueStandard::IsRun3(nr))
      {
        isRunPeriod3=true;
      }

      if(NueStandard::IsRun11(nr))
      {
        isRunPeriod11=true;
      }
      if(NueStandard::IsRun12(nr))
      {
        isRunPeriod12=true;
      }
      if(NueStandard::IsRun13(nr))
      {
        isRunPeriod13=true;
      }
    }
    
    if(separatebyRunPeriod)//recalculating pottree; this part works only for MC
    {
      if(nr->GetHeader().GetEventNo() == 0 || nr->GetHeader().GetEvents() == 0)//for each snarl
      {
        if(isRunPeriod1)
        {
          if(det==Detector::kNear)//this only works for ND MC
          {
            TotPOT1_nocut += MCInfo::GetMCPoT(det, beamtype, rel);
            TotPOT1 += MCInfo::GetMCPoT(det, beamtype, rel);
          }
          
	  nsnarls1++;//count the snarls per run period
        }
        if(isRunPeriod2)
        {
          if(det==Detector::kNear)//this only works for ND MC
          {
            TotPOT2_nocut += MCInfo::GetMCPoT(det, beamtype, rel);
            TotPOT2 += MCInfo::GetMCPoT(det, beamtype, rel);
          }
          
          nsnarls2++;
        }
        if(isRunPeriod3)
        {
          if(det==Detector::kNear)//this only works for ND MC
          {
            potweight = NueStandard::GetIntensityBeamWeight(nr);
            TotPOT3_nocut += (potweight*MCInfo::GetMCPoT(det, beamtype, rel));
            TotPOT3 += (potweight*MCInfo::GetMCPoT(det, beamtype, rel));
          }
          
          nsnarls3++;
        }


      if(isRunPeriod11)
        {
          if(det==Detector::kNear)//this only works for ND MC
          {
            TotPOT11_nocut += MCInfo::GetMCPoT(det, beamtype, rel);
            TotPOT11 += MCInfo::GetMCPoT(det, beamtype, rel);
          }
          
	  nsnarls11++;//count the snarls per run period
        }
        if(isRunPeriod12)
        {
          if(det==Detector::kNear)//this only works for ND MC
          {
            TotPOT12_nocut += MCInfo::GetMCPoT(det, beamtype, rel);
            TotPOT12 += MCInfo::GetMCPoT(det, beamtype, rel);
          }
          
          nsnarls12++;
        }
        if(isRunPeriod13)
        {
          if(det==Detector::kNear)//this only works for ND MC
          {
            TotPOT13_nocut += (MCInfo::GetMCPoT(det, beamtype, rel));
            TotPOT13 += (MCInfo::GetMCPoT(det, beamtype, rel));
          }
          
          nsnarls13++;
        }


      }
      
      if(nr->GetHeader().GetRun()!=lastrun)//counting runs - BUT not counting separately for each run period!!
      {
        lastrun = nr->GetHeader().GetRun();
        nruns++;
        
        if(det==Detector::kFar)//in FD MC, each run is a single file, and GetMCPoT returns pot/file
        {
          if(isRunPeriod1)
          {
            TotPOT1_nocut += MCInfo::GetMCPoT(det, beamtype, rel);
            TotPOT1 += MCInfo::GetMCPoT(det, beamtype, rel);
          }
          if(isRunPeriod2)
          {
            TotPOT2_nocut += MCInfo::GetMCPoT(det, beamtype, rel);
            TotPOT2 += MCInfo::GetMCPoT(det, beamtype, rel);
          }
          if(isRunPeriod3)
          {
            TotPOT3_nocut += MCInfo::GetMCPoT(det, beamtype, rel);
            TotPOT3 += MCInfo::GetMCPoT(det, beamtype, rel);
          }
          if(isRunPeriod11)
          {
            TotPOT11_nocut += MCInfo::GetMCPoT(det, beamtype, rel);
            TotPOT11 += MCInfo::GetMCPoT(det, beamtype, rel);
          }
          if(isRunPeriod12)
          {
            TotPOT12_nocut += MCInfo::GetMCPoT(det, beamtype, rel);
            TotPOT12 += MCInfo::GetMCPoT(det, beamtype, rel);
          }
          if(isRunPeriod13)
          {
            TotPOT13_nocut += MCInfo::GetMCPoT(det, beamtype, rel);
            TotPOT13 += MCInfo::GetMCPoT(det, beamtype, rel);
          }

        }
      }
    }
    
    if(fOverWritePOT && !isMC)//recalculating POT for data
    {
      if(NueStandard::PassesPOTStandards(nr))
      {
        if(nr->GetHeader().GetEventNo() == 0 || nr->GetHeader().GetEvents() == 0)
	{
	  if(det==Detector::kNear)
	  {
	    nddataweight = NueStandard::GetNDDataWeights(nr);
	    TotPOT1 += (nr->bmon.bI*nddataweight);
	  }
	  else//FD
	  {
	    TotPOT1 += nr->bmon.bI;
	  }
        }
      }
    }
    
    if(EvaluateCuts(nr)){
       if(fReweight && nr->mctrue.nuFlavor > -10)
       {
          int nuFlavor = nr->mctrue.nuFlavor;
          int  nonOsc = nr->mctrue.nonOscNuFlavor;
          float energy = nr->mctrue.nuEnergy;
                                                                                
          Float_t newWeight = NueConvention::Oscillate(nuFlavor, nonOsc, energy,735, fDeltaMSquare, fTheta23, fUe3Square);
                                                                                
          nr->mctrue.Ue3Squared = fUe3Square;
          nr->mctrue.DeltamSquared23 = fDeltaMSquare;
          nr->mctrue.Theta23  = fTheta23;
          nr->mctrue.fOscProb = newWeight;
       }
       
       if(separatebyRunPeriod)
       {
         if(isRunPeriod1)
         {
           tree1->Fill();
         }
         if(isRunPeriod2)
         {
           tree2->Fill();
         }
         if(isRunPeriod3)
         {
           tree3->Fill();
         }
         if(isRunPeriod11)
         {
           tree11->Fill();
         }
         if(isRunPeriod12)
         {
           tree12->Fill();
         }
         if(isRunPeriod13)
         {
           tree13->Fill();
         }
       }
       else
       {
         tree1->Fill();
       }
       count++;
    }
  }
  cout<<count<<" of "<<n<<" entries were passed"<<endl;
  
  f1->cd();
  NuePOT *total1 = new NuePOT();
  TTree *ptree1 = new TTree("pottree","pottree");
  TBranch* brp1 = ptree1->Branch("NuePOT","NuePOT", &total1 );
  brp1->SetAutoDelete(kFALSE);
  
  NuePOT *total2=0,*total3=0,*total11=0,*total12=0,*total13=0;
  TTree *ptree2=0,*ptree3=0,*ptree11=0,*ptree12=0,*ptree13=0;
  TBranch *brp2=0,*brp3=0,*brp11=0,*brp12=0,*brp13=0;
  if(separatebyRunPeriod)
  {
    f2->cd();
    total2 = new NuePOT();
    ptree2 = new TTree("pottree","pottree");
    brp2 = ptree2->Branch("NuePOT","NuePOT", &total2 );
    brp2->SetAutoDelete(kFALSE);
    
    f3->cd();
    total3 = new NuePOT();
    ptree3 = new TTree("pottree","pottree");
    brp3 = ptree3->Branch("NuePOT","NuePOT", &total3 );
    brp3->SetAutoDelete(kFALSE);


    f11->cd();
    total11 = new NuePOT();
    ptree11 = new TTree("pottree","pottree");
    brp11 = ptree11->Branch("NuePOT","NuePOT", &total11 );
    brp11->SetAutoDelete(kFALSE);

    f12->cd();
    total12 = new NuePOT();
    ptree12 = new TTree("pottree","pottree");
    brp12 = ptree12->Branch("NuePOT","NuePOT", &total12 );
    brp12->SetAutoDelete(kFALSE);
    
    f13->cd();
    total13 = new NuePOT();
    ptree13 = new TTree("pottree","pottree");
    brp13 = ptree13->Branch("NuePOT","NuePOT", &total13 );
    brp13->SetAutoDelete(kFALSE);


  }
  
  if(pchain->GetEntries() > 0){
    pchain->GetEntry(0);
    total1->beamtype = np->beamtype;
    if(separatebyRunPeriod)
    {
      total2->beamtype = np->beamtype;
      total3->beamtype = np->beamtype;
      total11->beamtype = np->beamtype;
      total12->beamtype = np->beamtype;
      total13->beamtype = np->beamtype;


    }
  }
  
  cout << pchain->GetEntries() << endl;
  cout << np->nsnarls << endl;

  for(int i = 0; i < pchain->GetEntries(); i++)
  {
     pchain->GetEntry(i);
     
     if(total1->beamtype != np->beamtype)
     {
       cerr<<"You are merging files of different beamtype - BAD"<<endl;
     }
     
     total1->nruns +=  np->nruns;
     total1->nsnarls +=  np->nsnarls;
     total1->pot += np->pot;
     total1->pot_nocut += np->pot_nocut;
  }

  if(fOverWritePOT && !isMC) total1->pot = TotPOT1;
  
  if(separatebyRunPeriod)//recalculating pottree
  {
    total1->pot = TotPOT1;
    total2->pot = TotPOT2;
    total3->pot = TotPOT3;
    
    total1->pot_nocut = TotPOT1_nocut;
    total2->pot_nocut = TotPOT2_nocut;
    total3->pot_nocut = TotPOT3_nocut;
    
    total1->nruns = nruns;
    total2->nruns = nruns;
    total3->nruns = nruns;
    
    total1->nsnarls = nsnarls1;
    total2->nsnarls = nsnarls2;
    total3->nsnarls = nsnarls3;


    total11->pot = TotPOT11;
    total12->pot = TotPOT12;
    total13->pot = TotPOT13;
    
    total11->pot_nocut = TotPOT11_nocut;
    total12->pot_nocut = TotPOT12_nocut;
    total13->pot_nocut = TotPOT13_nocut;
    
    total11->nruns = nruns;
    total12->nruns = nruns;
    total13->nruns = nruns;
    
    total11->nsnarls = nsnarls11;
    total12->nsnarls = nsnarls12;
    total13->nsnarls = nsnarls13;


  }
  
  //fill the pottree
  ptree1->Fill();
  if(separatebyRunPeriod)
  {
    ptree2->Fill();
    ptree3->Fill();

    ptree11->Fill();
    ptree12->Fill();
    ptree13->Fill();


  }
  
  //save the file(s)
  f1->cd();
  tree1->Write("ana_nue",TObject::kWriteDelete);
  ptree1->Write("pottree",TObject::kWriteDelete);
  if(separatebyRunPeriod)
  {
    f2->cd();
    tree2->Write("ana_nue",TObject::kWriteDelete);
    ptree2->Write("pottree",TObject::kWriteDelete);
    
    f3->cd();
    tree3->Write("ana_nue",TObject::kWriteDelete);
    ptree3->Write("pottree",TObject::kWriteDelete);


    f11->cd();
    tree11->Write("ana_nue",TObject::kWriteDelete);
    ptree11->Write("pottree",TObject::kWriteDelete);

    f12->cd();
    tree12->Write("ana_nue",TObject::kWriteDelete);
    ptree12->Write("pottree",TObject::kWriteDelete);
    
    f13->cd();
    tree13->Write("ana_nue",TObject::kWriteDelete);
    ptree13->Write("pottree",TObject::kWriteDelete);

  }
  
  f1->Close();
  if(separatebyRunPeriod)
  {
    f2->Close();
    f3->Close();

    f11->Close();
    f12->Close();
    f13->Close();

  }
  
  
  cout<<"Trimming completed with "<<np->pot<<"x10^12 POT"<<endl;
}

void Trimmer::SetCuts(string type, int level)
{
   cutSet = type;
   cutLevel = level;   
}

bool Trimmer::EvaluateCuts(NueRecord *nr)
{
	//just for merging files!
   if(cutSet == "Merge")
	return true;

   if(cutSet == "Fiducial") 
     return NueStandard::IsInFid(nr);

   if(cutSet == "Standard"){
      bool ret;
      switch(cutLevel){
        case 0: return NueStandard::PassesSelection(nr, Selection::kDataQual);
        case 1: return NueStandard::PassesSelection(nr, Selection::kFid);
        case 2: return NueStandard::PassesSelection(nr, Selection::kBasic);
        case 3: 
             ret = NueStandard::PassesSelection(nr, Selection::kBasic);
             ret = ret && NueStandard::PassesPreSelectionTrackCuts(nr);
             return ret;
        case 4: 
             ret = NueStandard::PassesSelection(nr, Selection::kBasic);
             ret = ret &&  NueStandard::PassesNonHEPreSelection(nr);
             return ret;
        case 5: return NueStandard::PassesSelection(nr, Selection::kPre);
      }
   }
   if(cutSet == "MRE" || cutSet == "MRCC"){
      bool dq = NueStandard::PassesSelection(nr, Selection::kDataQual);
      bool mrefid = NueStandard::PassesMREFiducial(nr);
      bool ret;
      switch(cutLevel){
        case 0: return dq && NueStandard::IsInFid(nr);
        case 1: return mrefid && dq && NueStandard::IsInFid(nr);
        case 2: 
          return NueStandard::PassesSelection(nr, Selection::kFid);//includes mrcc presel & mrcc fiducial
        case 3: return NueStandard::PassesSelection(nr,Selection::kBasic);
        case 4:
             ret = NueStandard::PassesSelection(nr,Selection::kBasic);
             ret = ret && NueStandard::PassesPreSelectionTrackCuts(nr);
             return ret;
        case 5:
             ret = NueStandard::PassesSelection(nr, Selection::kBasic);
             ret = ret &&  NueStandard::PassesNonHEPreSelection(nr);
             return ret;
        case 6: 
         return NueStandard::PassesSelection(nr, Selection::kPre);
      }
   }
   if(cutSet == "Presel")
      return NueStandard::PassesSelection(nr, Selection::kPre);

   if(cutSet == "PreselNoHE")
      return  NueStandard::PassesNonHEPreSelection(nr);

   if(cutSet == "Systematic"){
      // O - standard systematic envelope
      // 1 - standard systematic w/out HE cut
      if(!NueStandard::IsInFid(nr)) return false;

      switch(cutLevel){
        case 0: return NueStandard::PassesSysPreSelection(nr);
        case 1: return NueStandard::PassesSysPreSelectionNoHE(nr);
      }
   }


   if(cutSet == "CC")
      return NueStandard::PassesSelection(nr, Selection::kCC);
   
   if(cutSet == "GoodRun")
     return NueStandard::IsGoodRun(nr);
   
   if(cutSet == "AntiANN11")
   {
     if(NueStandard::PassesSelection(nr, Selection::kPre) && NueStandard::GetPIDValue(nr,Selection::kANN2PE)<0.55)
     {
       return true;
     }
     else
     {
       return false;
     }
   }
   
   if(cutSet == "AntiMCNN")
   {
     if(NueStandard::PassesSelection(nr, Selection::kPre) && NueStandard::GetPIDValue(nr,Selection::kMCNN)<0.55)
     {
       return true;
     }
     else
     {
       return false;
     }
   }
   if(cutSet=="MRCuts")
   {
     if(NueStandard::PassesSelection(nr, Selection::kDataQual) && NueStandard::PassesMRCCFiducial(nr) && NueStandard::PassesMRCCPreSelection(nr))
     {
       return true;
     }
     return false;
   }
   
  cout<<"Invalid Cut Level for "<<cutSet<<": "<<cutLevel<<endl;
  return false;
}
