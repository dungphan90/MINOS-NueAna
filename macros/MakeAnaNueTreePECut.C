// Macro for making AnaNue trees
//-----------------------------------------------------------------------------
// Can be run in several ways:
//----------------------------
// SNTP files (Birch/Carrot or Cedar/Daikon):
//      loon NueAna/macros/MakeAnaNueTree.C blah.sntp.root
// or
//      loon 'NueAna/macros/MakeAnaNueTree.C("blah.sntp.root")'
//
// MRNT files (Birch/Carrot single file era):
//      loon NueAna/macros/MakeAnaNueTree.C blah.mrnt.root
// or
//      loon 'NueAna/macros/MakeAnaNueTree.C("blah.mrnt.root")'
//
// SNTP + MRNT files (Cedar/Daikon split file era):
//      loon 'NueAna/macros/MakeAnaNueTree.C("blah.sntp.root","blah.mrnt.root")'
//
// MCNN (LEM) file available for merging:
//      loon 'NueAna/macros/MakeAnaNueTree.C("blah.sntp.root","","mcnn_matches.root")'
// or                                                           
//      loon 'NueAna/macros/MakeAnaNueTree.C("blah.sntp.root","blah.mrnt.root","mcnn_matches.root")'
//
//------------------------------------------------------------------------------

void MakeAnaNueTreePECut(string sntp_file="",string mrnt_file="",string mcnn_file=""){

  // resolve any environment variables
  sntp_file=string(gSystem->ExpandPathName(sntp_file.c_str()));
  mrnt_file=string(gSystem->ExpandPathName(mrnt_file.c_str()));
  mcnn_file=string(gSystem->ExpandPathName(mcnn_file.c_str()));

  ////////////////////////////////////////////////////////
  //load necessary libraries
  gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs();
  
  JobC j;

  ////////////////////////////////////////////////////////
  //set messaging levels
  //change these to Debug, Info, Warning, Error to set different levels of output
  j.Msg.SetLevel("NueModule","Info");
  j.Msg.SetLevel("NueRecord","Info");
  j.Msg.SetLevel("ShwfitCalc","Info");
  j.Msg.SetLevel("HitCalcAna","Error");
  j.Msg.SetLevel("AngClusterAna","Error");
  j.Msg.SetLevel("AngClusterFitAna","Error");
  j.Msg.SetLevel("ANtpAnalysisInfoAna","Info");
  j.Msg.SetLevel("ANtpEventInfoAna","Debug");
  j.Msg.SetLevel("ANtpShowerInfoAna","Info");
  j.Msg.SetLevel("ANtpTrackInfoAna","Info");
  j.Msg.SetLevel("MSTCalcAna","Info");
  j.Msg.SetLevel("BeamSummary","Error");
  j.Msg.SetLevel("ANtpStNtpManipulator","Error");
  j.Msg.SetLevel("ANtpAnalysisInfoAna","Error");
  j.Msg.SetLevel("ANtpEventInfoAna","Error");
  j.Msg.SetLevel("FillAlg","Error");
//  j.Msg.SetLevel("NueModule","Debug");
  //j.Msg.SetLevel("NueRecord","Debug");
  //j.Msg.SetLevel("ShwfitAna","Debug");
  //j.Msg.SetLevel("BeamMonAna","Debug");
  //j.Msg.SetLevel("NueBeamMon","Debug");
  //j.Msg.SetLevel("MCFluxInfoAna","Debug");
  //j.Msg.SetLevel("NueFluxWeightsAna","Debug");
  //j.Msg.SetLevel("NueFluxWeights","Debug");
 j.Msg.SetLevel("PETrimmer","Warning");  //Debug, Fatal 
  ////////////////////////////////////////////////////////
  //set the analysis path
  JobCPath& main_path = j.Path.Create("Main",
                                      "SetKNNModule::Reco "
				      "PETrimmer::Reco "
				      "NueModule::Reco "
				      "MCNNMergeModule::Reco "
				      "NueBeamMon::Reco "
				      "Output::Put");
  
  ////////////////////////////////////////////////////////
  //Set Input Streams
  j.Input.Set("Streams=NtpSt,NtpSR,NtpMC,NtpTH,NtpMR,NtpOld,snarltree"); 

  if(sntp_file != "") {
    if(mrnt_file == "") {
      j.Input.AddFile(sntp_file.c_str(),"NtpSt");
    }
    else {
      j.Input.DefineStream("NtpOld","NtpSt");
      j.Input.AddFile(sntp_file.c_str(),"NtpOld");
      j.Input.AddFile(mrnt_file.c_str(),"NtpSt");
      j.Input.AddFile(mrnt_file.c_str(),"NtpMR");
    }
  }
  
  ////////////////////////////////////////////////////////
  //configure NueModule

  Registry r(false);
  r.Set("kNNPath", "/minos/app/nue/Releases/Griffin/Data/kNN_cedar/");
  main_path.Mod("SetKNNModule").Config(r);

  JobCModule& pe_mod = main_path.Mod("PETrimmer"); 
  float pecut=2.0;
  char tmp[100];
  sprintf(tmp,"PECut=%f",pecut);
  pe_mod.Set(tmp);// set the pe cut
  pe_mod.Set("updateEventEnergy=0");

  JobCModule& nue_mod = main_path.Mod("NueModule");

  //set filename of template hists for MST
  char *tmpltdir = getenv("SRT_PRIVATE_CONTEXT");
  if(tmpltdir=='\0'){
    tmpltdir="./";
  }
  char settmpname[10000];
  sprintf(settmpname,"MSTTmpltFile=%s/NueAna/data/templates.root",tmpltdir);
  nue_mod.Set(settmpname);
  
  //Set options for MDA posterior probability calculation
  //probability threshold:
  nue_mod.Set("MdaThreshCut=0.90");
  //location of SAS calibration coefficients.
  char *path = getenv("SRT_LOCAL");
  char MdaFile[1000];
  sprintf(MdaFile,"MdaSASFile=%s/NueAna/data/Mda_Coeff_Far_CedarDaikon_Chimaera.dat",path);
  nue_mod.Set(MdaFile);  

  // Set MCNN input file from argument; set other MCNN parameters.
  // Just copying the kNN approach here of making a Registry object.
  // The PDF file should probably become a CVS "proxy" file.  For now,
  // it's just sitting in /minos/app.
  Registry rr(false);
  rr.Set("InputMCNNFile",mcnn_file.c_str());
  rr.Set("IsLinfixMC", 0);
  rr.Set("NumBestMatches", 50);
  rr.Set("YCut", 0.9);
  rr.Set("PDFFile", "/minos/app/nue/Releases/Griffin/Data/LEM/pdf/MCNNpdf_2D_0.5GeVbins_MCNNv2_3.0pe_3.0pelib.root");
  j.Path("Main").Mod("MCNNMergeModule").Config(rr);

  ////////////////////////////////////////////////////////
  //Configure output
  
  //set the streams that will be written
  j.Path("Main").Mod("Output").Cmd("DefineStream ana_nue NueRecord");
  j.Path("Main").Mod("Output").Set("Streams=ana_nue");
  
  //Get the run number from the first Input file name
  //doesn't have to be this way, I just thought I'd be tricky
  path = getenv("NUE_OUT_DIR");
  if(path=='\0'){
    path="./";
  }
  string filename;
  if(mrnt_file != "") filename = mrnt_file;
  else if(sntp_file != "") filename = sntp_file;
  else { 
    JobCEnv& jce = JobCEnv::Instance();
    filename = jce.GetFileName(0);
  }
  string rnumstring = filename.substr(filename.find_last_of("/")+2,8);
  string subrunstring = filename.substr(filename.find_last_of("/")+11,4);
  
  string file = filename.substr(filename.find_last_of("/")+1, 
				filename.find_last_of(".root") - 
				filename.find_last_of("/")-5);
  
  char setfilename[5000];
  sprintf(setfilename,"FileName=%sAnaNue-%s%s-PECut-%f.root",
	  path,((mcnn_file=="")?"":"mcnn2-"),file.c_str(),pecut);
  cout<<"Saving AnaNue Tree to "<<setfilename<<endl;
  j.Path("Main").Mod("Output").Set(setfilename);
  
  char pottreefn[1000];
  sprintf(pottreefn,"POTTreeFileName=%sAnaNue-%s%s-PECut-%f.root",path,((mcnn_file=="")?"":"mcnn2-"),file.c_str(),pecut);
  j.Path("Main").Mod("NueBeamMon").Set(pottreefn);
  
  ////Configure Flux Reweighting objects
  //figure out beam type from filename
  int beam=2;
  
  //figure out which flux file directory to use from beam
  //set fixmuparents to non zero if you want to use the mupi trees to change the
  //tp info for muons.  set it to zero if you dont.
  j.Path("Main").Mod("NueModule").Set("FixMuParents=0");
  //string ffbasepath="MuPiDir=/stage/minos-data6/vahle/MuParTrees/fluka05_"+flxsfx;
  //string ffbasepath="MuPiDir=/stage/minos-data6/cbs/gnumi/v18/muonInfo/fluka05_"+flxsfx;
  string ffbasepath = "MuPiDir=/afs/fnal.gov/files/data/minos/d19/muonInfo/v18/fluka05_";
  
  cout<<"ffbasepath: "<<ffbasepath<<endl;
  j.Path("Main").Mod("NueModule").Set(ffbasepath.c_str());
 

  j.Path("Main").Mod("NueModule").Set("ContPhPlaneCut=0.3"); //contiguous plane ph threshold
    
  ////////////////////////////////////////////////////////
  //Run
  //j.Input.Next(74);  
  j.Path("Main").Run();
    
  ////////////////////////////////////////////////////////
  //report
  j.Path("Main").Report();
  
  ////////////////////////////////////////////////////////
   
}
