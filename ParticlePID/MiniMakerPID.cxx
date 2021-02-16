#include "NueAna/ParticlePID/MiniMakerPID.h"
#include "NueAna/NueStandard.h"
#include "NueAna/ParticlePID/NueMiniAnaPID.h"

#include "DataUtil/MCInfo.h"
#include "DataUtil/DataQualDB.h"
#include "NueAna/AnnAna.h"

MiniMakerPID::MiniMakerPID(){ 
  outSet = false;
  fReweight = false;
  fBaseline = 735;
  fDeltaMSquare = 0.0027;
  fUe3Square = 0.025;
  fTheta23 = TMath::Pi()/4;

  fOverWritePOT = false;
  
  fIsMRCC=0;
}

void MiniMakerPID::AddFiles(string infiles)
{
   files.push_back(infiles);
}

void MiniMakerPID::SetOutputFile(string file)
{
   outSet = true;
   outFile = file;
}   

void MiniMakerPID::Config(Registry &r)
{
  fCuts.Config(r);
}

void MiniMakerPID::SetDeltaMSquare(double dm2) { fDeltaMSquare = dm2; fReweight = true;}
void MiniMakerPID::SetUe3Square(double dUe32) { fUe3Square = dUe32; fReweight = true;}
void MiniMakerPID::SetTheta23(double t23) { fTheta23 = t23; fReweight = true;}
void MiniMakerPID::SetBaseline(double bl) { fBaseline = bl; fReweight = true;}

void MiniMakerPID::RunMiniMakerPID()
{
  NueRecord *nr = new NueRecord();
  NuePOT *np = new NuePOT();
  NuePOT *total = new NuePOT();

  TChain *chain = new TChain("ana_nue");
  chain->SetBranchAddress("NueRecord",&nr);
                                                                        
  TChain *pchain = new TChain("pottree");
  pchain->SetBranchAddress("NuePOT", &np);
 
  int nchain=0;
  int npchain=0;

  for(unsigned int i = 0; i < files.size(); i++)
  {
      nchain += chain->Add(files[i].c_str());
      npchain += pchain->Add(files[i].c_str());
  }


  if(nchain != npchain || !nchain || !npchain)
  {   
      printf("missing tree in input file.... exiting with error!\n");
      exit(1);
  }

  if(pchain->GetEntries() > 0){
  	pchain->GetEntry(0);
    total->beamtype = np->beamtype;
  }
  
  if(!outSet){
     string file;
                                                                        
     if(pchain->GetEntries() > 0){
        pchain->GetEntry(0);
        total->beamtype = np->beamtype;
        file = pchain->GetFile()->GetName();
     }
     string minifile = file.substr(file.find_last_of("/")+1, file.find_last_of(".root")-file.find_last_of("/") - 5);
     minifile += "-Mini.root";
  
     outFile = minifile;
  }
  if(outFile == "-Mini.root"){
     cout<<"No input file found"<<endl; 
     return;
  }
                                                                         
  cout<<"Setting output to "<<outFile<<endl;

  // And now the actual looping
                                                                                
  TFile *save = new TFile(outFile.c_str(),"RECREATE");
  save->SetCompressionLevel(9);
  save->cd();
  TTree *tree = new TTree("NueMiniPID","NueMiniPID");
  TTree::SetBranchStyle(1);

  if(chain->GetEntries() > 0) chain->GetEvent(0);
  Detector::Detector_t det = nr->GetHeader().GetVldContext().GetDetector();
  ReleaseType::Release_t rel = nr->GetHeader().GetRelease();
  BeamType::BeamType_t beam = nr->GetHeader().GetBeamType();
                                                                                
  NueMiniAnaPID nma;
                                                                                
  NueMiniPID* nm = new NueMiniPID(beam, det, rel);

  TBranch* br = tree->Branch("NueMiniPID","NueMiniPID", &nm );
  br->SetAutoDelete(kFALSE);
                                                                        
  Int_t n = chain->GetEntries();
  int count = 0;

  bool isMC = false;
  double TotPOT = 0.0;

	NueStandard::SetDefaultOscParam();
	NueStandard::SetOscNoMatter(); 

  for(int i=0;i<n;i++){
    if(i%10000==0) cout << 100*float(i)/float(n)
                       << "% done" << endl;
    chain->GetEvent(i);

    if(i == 0){
      SimFlag::SimFlag_t sim = nr->GetHeader().GetVldContext().GetSimFlag();
      isMC = (sim == SimFlag::kMC);
    }       

	//reset it here
	nm->fPOT=0;



    nr->eventq.passFarDetQuality = DataUtil::IsGoodData(nr->GetHeader().GetVldContext());
    nr->eventq.passNearDetQuality = nr->eventq.passFarDetQuality;
    if( NueStandard::PassesPOTStandards(nr))
	{
    //if(NueStandard::PassesDataQuality(nr)) {
       if(nr->GetHeader().GetEventNo() == 0 ||
             nr->GetHeader().GetEvents() == 0){
          	TotPOT += nr->bmon.bI;
        	nm->fPOT=nr->bmon.bI;
	  if(nr->bmon.bI < 0 && !isMC)
	  {
	 	cout<<"Unexpected behavior"<<endl;
		exit(2);
	  }       
	}
    }


        //recalculate the ann to get ann14
        AnnAna z(*nr);
        z.Analyze();


	//fill it now, so the variables are available to EvaluateCuts
	nma.FillMini(nr, nm, fIsMRCC);

	
    if(EvaluateCuts(nr,nm)){

       if(fReweight && nr->mctrue.nuFlavor > -10)
       {
          //int nuFlavor = nr->mctrue.nuFlavor;
          //int  nonOsc = nr->mctrue.nonOscNuFlavor;
          //float energy = nr->mctrue.nuEnergy;
                                                                                
          Float_t newWeight = NueStandard::GetOscWeight(nr->mctrue.nuFlavor,nr->mctrue.nonOscNuFlavor,nr->mctrue.nuEnergy);//NueConvention::Oscillate(nuFlavor, nonOsc, energy,                                735, fDeltaMSquare, fTheta23, fUe3Square);
                                                                                
          nr->mctrue.Ue3Squared = fUe3Square;
          nr->mctrue.DeltamSquared23 = fDeltaMSquare;
          nr->mctrue.Theta23  = fTheta23;
          nr->mctrue.fOscProb = newWeight;
		
		  //fill again to pick up oscillation changes
	      nma.FillMini(nr, nm, fIsMRCC);
       }
       
		
       tree->Fill();
       count++;
    }else{
 /*    float x = nr->srevent.vtxX;
      float y = nr->srevent.vtxY;
      float r = TMath::Sqrt(x*x+y*y);

      cout<<"rejecting event"<<x<<"  "<<y<<"  "<<r<<"  "
          <<nr->srevent.vtxZ<<"  "
          <<NueStandard::PassesSelection(nr, Selection::kDataQual)
          <<"  "<<NueStandard::PassesSelection(nr, Selection::kFid)<<endl;
*/
    }
  }
  cout<<count<<" of "<<n<<" entries were passed"<<endl;
                                                                        
  TTree *ptree = new TTree("pottree","pottree");
  TBranch* br2 = ptree->Branch("NuePOT","NuePOT", &total );
  br2->SetAutoDelete(kFALSE);
                                                                        
  for(int i = 0; i < pchain->GetEntries(); i++)
  {
     pchain->GetEntry(i);
        
       cout <<"bt "<<total->beamtype <<" "<<np->beamtype<<"\n";                                                                 
     if(total->beamtype != np->beamtype)
        cerr<<"You are merging files of different beamtyp - BAD"<<endl;
                                                                        
     total->nruns +=  np->nruns;
     total->nsnarls +=  np->nsnarls;
     total->pot += np->pot;
  }

  if(fOverWritePOT && !isMC) total->pot = TotPOT;

  ptree->Fill();
                                                                        
  save->cd();
  tree->Write();
  ptree->Write();
 
  cout<<"Compression level of output file "<<save->GetCompressionLevel()<<"\n";
  
  cout<<"Trimming completed with "<<total->pot<<"x10^12 POT"<<endl;
  save->Close();
}

void MiniMakerPID::SetCuts(string type, int level)
{
   cutSet = type;
   cutLevel = level;   
}

bool MiniMakerPID::EvaluateCuts(NueRecord *nr, NueMiniPID *nm)
{
   NueConvention::NueEnergyCorrection(nr);  

   if(cutSet == "Fiducial") 
     return NueStandard::PassesSelection(nr, Selection::kFid);

   if(cutSet == "Standard"){
      bool ret;
      switch(cutLevel){
        case 0: return NueStandard::PassesSelection(nr, Selection::kDataQual);
        case 1: return NueStandard::PassesSelection(nr, Selection::kFid);
        case 2: 
             ret = NueStandard::PassesSelection(nr, Selection::kFid);
             ret = ret &&  NueStandard::PassesNonHEPreSelection(nr);
             return ret;
        case 3: return NueStandard::PassesSelection(nr, Selection::kPre);
	case 4: return nm->infid && nm->contained && ( fIsMRCC ? nm->mrcc_s : 1 );
	case 6: return NueStandard::PassesSelection(nr, Selection::kDataQual) && ((nm->infid && nm->contained ) || NueStandard::PassesSelection(nr, Selection::kFid) ) || nm->fPOT;


      }
   }
   if(cutSet == "MRE" || cutSet == "MRCC"){
      bool dq = NueStandard::PassesSelection(nr, Selection::kDataQual);
      bool mrefid = NueStandard::PassesMREFiducial(nr);
      bool mrepre = NueStandard::PassesMREPreSelection(nr) && mrefid && dq;
      bool ret;
      switch(cutLevel){
        case 0: return mrefid && dq;
        case 1: return mrepre;
        case 2: 
          return NueStandard::PassesSelection(nr, Selection::kFid) && mrepre;
        case 3:
             ret = NueStandard::PassesSelection(nr, Selection::kFid);
             ret = ret &&  NueStandard::PassesNonHEPreSelection(nr);
             return ret && mrepre;
        case 4: 
         return NueStandard::PassesSelection(nr, Selection::kPre) && mrepre;
      }
   }
   if(cutSet == "Presel")
      return NueStandard::PassesSelection(nr, Selection::kPre);

   if(cutSet == "PreselNoHE")
      return  NueStandard::PassesNonHEPreSelection(nr);

   if(cutSet == "Systematic")
      return NueStandard::PassesSysPreSelection(nr);

   if(cutSet == "CC")
      return NueStandard::PassesSelection(nr, Selection::kCC);


  cout<<"Invalid Cut Level for "<<cutSet<<": "<<cutLevel<<endl;
  return false;
}
