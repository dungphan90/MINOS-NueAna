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

void RunNueDisplay(string sntp_file="",string mrnt_file="")
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

  //gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadDisplayLibs();

  ifstream Test("$SRT_PRIVATE_CONTEXT/NueAna/macros/LoadLibs.C");
  if (Test){
    gROOT->LoadMacro("$SRT_PRIVATE_CONTEXT/NueAna/macros/LoadLibs.C"); 
  }
  else{
    gROOT->LoadMacro("$SRT_PUBLIC_CONTEXT/NueAna/macros/LoadLibs.C"); 
  }

  LoadDisplayLibs();

  jc = new JobC;
  mint = new Mint(*jc);

  const int width = 1000, height = 680;

  jc->Msg.SetLevel("NueDisplayModule","Info");

  // To configure the PageDisplay do something like:
  // PageDisplayConfig pdc;
  // pdc.use_run_snarl_entry = false;
  // PageDisplay* pd = mint->SpawnDisplay(width,height,&pdc);
  //
  // Or, just don't:
  PageDisplay* pd = mint->SpawnDisplay(width,height);

  jc->Path.Create("default","NueDisplayModule::Ana");
  jc->Input.Set("Format=input");
  jc->Input.Set("Streams=NtpSt,NtpMR,NtpOld,NtpSR,NtpMC,NtpTH,ana_nue");

  if(sntp_file != "") {
    if(mrnt_file == "") {
      jc->Input.AddFile(sntp_file.c_str(),"NtpSt");
    }
    else {
      jc->Input.DefineStream("NtpOld","NtpSt");
      jc->Input.AddFile(sntp_file.c_str(),"NtpOld");
      jc->Input.AddFile(mrnt_file.c_str(),"NtpSt");
      jc->Input.AddFile(mrnt_file.c_str(),"NtpMR");
    }
  }

  //set threshold for # of planes in track, above which no nue anal will be done.
  //cuts are done < or > not equal , energy is sigcor based
  // Cuts on max number of planes for the track 
  jc->Path("default").Mod("NueDisplayModule").Set("HiPlaneTrackCut=18");
  // Cuts on max number of tracklike planes
  //jc->Path("default").Mod("NueDisplayModule").Set("HiTrackLikeCut=18");
  
  // Cuts on min number of planes for the event
  //jc->Path("default").Mod("NueDisplayModule").Set("LoPlaneEventCut=-1");
  
  // Cuts on max and min total energy in sigcor
  //jc->Path("default").Mod("NueDisplayModule").Set("HiEnergyCut=-1");
  //jc->Path("default").Mod("NueDisplayModule").Set("LoEnergyCut=-1");
  
  // Cuts on max and min total energy for the shower in GeV
  // note that the min cut, is the cut for a single shower, is several
  // showers are present then the cut is half that value
  //jc->Path("default").Mod("NueDisplayModule").Set("HiEnergyShowerCut=-1");
  jc->Path("default").Mod("NueDisplayModule").Set("LoEnergyShowerCut=1.5");
  
  // cuts to produce hiPhStripCount and hiPhPlaneCount
  jc->Path("default").Mod("NueDisplayModule").Set("DPlaneCut=5"); // total number of planes to count strips/planes
  // A negative number in this variable will result zero counts
  jc->Path("default").Mod("NueDisplayModule").Set("PhStripCut=0.2"); // strip ph threshold

  jc->Path("default").Mod("NueDisplayModule").Set("PhPlaneCut=2"); //plane ph threshold
  jc->Path("default").Mod("NueDisplayModule").Set("ContPhPlaneCut=0.3"); //contiguous plane ph threshold

  jc->Path("default").Mod("NueDisplayModule").Set("LoPhNStripCut=-1"); // min number of strips in hiPhStripCount
  jc->Path("default").Mod("NueDisplayModule").Set("LoPhNPlaneCut=-1"); // min number of planes in hiPhPlaneCount
  
  // Use with PID file, to select events above a certain IsNue threshold.
  //jc->Path("default").Mod("NueDisplayModule").Set("PIDCut=0.5");
  
  //by changing ScanMode from 0 to 1, 
  //you don't have to click 'Log Details' to record your decisions.
  //jc->Path("default").Mod("NueDisplayModule").Set("ScanMode=1");
  
  //Switching on test mode removes the truth buttons and hides run/snarl info 
  //from the text box 
  //jc->Path("default").Mod("NueDisplayModule").Set("TestMode=1");
  
  //To just hide run/snarl info, but leave truth buttons
  //jc->Path("default").Mod("NueDisplayModule").Set("HideRunSnarl=1");
  
  //To switch on cluster info (not recommended for older ntuples)
  //jc->Path("default").Mod("NueDisplayModule").Set("DrawClu=1");
  
  //To switch on Interactive Reco
  //jc->Path("default").Mod("NueDisplayModule").Set("IntReco=1");
  
  jc->Path("default").Run(1);

}

     
