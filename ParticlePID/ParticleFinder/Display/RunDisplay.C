// Nue Display Macro
//-----------------------------------------------------------------------------
// Can be run in several ways:
//----------------------------
// SNTP files (Birch/Carrot or Cedar/Daikon):
//      loon NueAna/macros/RunNueDisplay.C blah.sntp.root
// or
//      loon 'NueAna/macros/RunNueDisplay.C("blah.sntp.root")'
//
// MRNT files (Birch/Carrot single file era):
//      loon NueAna/macros/RunNueDisplay.C blah.mrnt.root
// or
//      loon 'NueAna/macros/RunNueDisplay.C("blah.mrnt.root")' 
//
// SNTP + MRNT files (Cedar/Daikon split file era):
//      loon 'NueAna/macros/RunNueDisplay.C("blah.sntp.root","blah.mrnt.root")'
//
//------------------------------------------------------------------------------

class JobC;
class Mint;

JobC *jc = 0;
Mint *mint = 0;

void RunDisplay(string sntp_file="",string mrnt_file="")
{

 


  gStyle->SetPalette(1,0);
  //color code:
  //http://root.cern.ch/root/html/TColor.html#TColor:description
  int colorcode = 10; //10:white  21:grey
  gStyle->SetCanvasColor(colorcode);
  gStyle->SetPadColor(colorcode);
  gStyle->SetFrameFillColor(colorcode);
  gStyle->SetTitleFillColor(colorcode);
  gStyle->SetStatColor(colorcode);

  gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadDisplayLibs();


 
  	






  

  jc = new JobC;
  mint = new Mint(*jc);

  const int width = 1000, height = 680;

  jc->Msg.SetLevel("ParticleDisplay","Error");//Debug");
  jc->Msg.SetLevel("ParticleFinder","Error");
  jc->Msg.SetLevel("LongMuonFinder","Error"); 
  jc->Msg.SetLevel("ClusterManager","Error");
  jc->Msg.SetLevel("HitManager","Error");
  jc->Msg.SetLevel("Finder","Error");
  jc->Msg.SetLevel("Chain","Error");
  jc->Msg.SetLevel("PrimaryShowerFinder","Error");

  jc->Msg.SetLevel("ChainHelper","Error");

  jc->Msg.SetLevel("ShwFit","Error");
 
  // To configure the PageDisplay do something like:
  // PageDisplayConfig pdc;
  // pdc.use_run_snarl_entry = false;
  // PageDisplay* pd = mint->SpawnDisplay(width,height,&pdc);
  //
  // Or, just don't:
  PageDisplay* pd = mint->SpawnDisplay(width,height);

  jc->Path.Create("default","ParticleFinder::Reco  ParticleAna::Reco ParticleDisplay::Ana");
  jc->Input.Set("Format=input");
  jc->Input.Set("Streams=PO,PA,NtpSt");

      jc->Input.AddFile(sntp_file.c_str(),"NtpSt");
  

   
  
  jc->Path("default").Run(1);

}

     
