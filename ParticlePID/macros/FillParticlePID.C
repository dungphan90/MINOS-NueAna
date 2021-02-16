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

void FillParticlePID(string sntp_file=""){

  // resolve any environment variables
  sntp_file=string(gSystem->ExpandPathName(sntp_file.c_str()));
 
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
  j.Msg.SetLevel("XTalkFilter","Error");
//  j.Msg.SetLevel("NueModule","Debug");
  //j.Msg.SetLevel("NueRecord","Debug");
  //j.Msg.SetLevel("ShwfitAna","Debug");
  //j.Msg.SetLevel("BeamMonAna","Debug");
  //j.Msg.SetLevel("NueBeamMon","Debug");
  //j.Msg.SetLevel("MCFluxInfoAna","Debug");
  //j.Msg.SetLevel("NueFluxWeightsAna","Debug");
  //j.Msg.SetLevel("NueFluxWeights","Debug");
 j.Msg.SetLevel("PETrimmer","Warning");  //Debug, Fatal 




  j.Msg.SetLevel("ParticleFinder","Error");
  j.Msg.SetLevel("LongMuonFinder","Error");
  j.Msg.SetLevel("ClusterManager","Error");
  j.Msg.SetLevel("HitManager","Error");
  j.Msg.SetLevel("Finder","Error");
  j.Msg.SetLevel("Chain","Error");
  j.Msg.SetLevel("PrimaryShowerFinder","Error");
  j.Msg.SetLevel("ShwFit","Error");
  j.Msg.SetLevel("ChainHelper","Error");
  j.Msg.SetLevel("TruthCompareAna","Error");






  ////////////////////////////////////////////////////////
  //set the analysis path
  	JobCPath& main_path = j.Path.Create("Main");

	JobCPath& nueana = j.Path.Create("NueAna",	
				"PIDEval::Reco "//fill steves pids
			      	"Output::Put"
	);
	main_path.Attach(&nueana);
  


	////////////////////////////////////////////////////////
  //Set Input Streams
  j.Input.Set("Streams=ana_nue"); 

 	  
  ////////////////////////////////////////////////////////
  //configure NueModule

  ////////////////////////////////////////////////////////
  //Configure output
  
  //set the streams that will be written
  j.Path("NueAna").Mod("Output").Cmd("DefineStream ana_nue NueRecord");
 // j.Path("NueAna").Mod("Output").Cmd("DefineStream  PO ParticleObjectHolder");
 // j.Path("NueAna").Mod("Output").Cmd("DefineStream  PA PRecord");
  
  j.Path("NueAna").Mod("Output").Set("Streams=ana_nue");
  
  
  
  //Get the run number from the first Input file name
  //doesn't have to be this way, I just thought I'd be tricky
  path = getenv("NUE_OUT_DIR");
  if(path=='\0'){
    path="./";
  }
  string filename;
  if(sntp_file != "") filename = sntp_file;
  else { 
    JobCEnv& jce = JobCEnv::Instance();
    filename = jce.GetFileName(0);
  }
  string rnumstring = filename.substr(filename.find_last_of("/")+2,8);
  string subrunstring = filename.substr(filename.find_last_of("/")+11,4);
  
  string file = filename.substr(filename.find_last_of("/")+1, 
				filename.find_last_of(".root") - 
				filename.find_last_of("/")-5);
  

  char ext[100];
  //sprintf(ext,"-PIDFILL");
  sprintf(ext,"");


  char setfilename[5000];
  sprintf(setfilename,"FileName=%s%s%s.root",
	  path,file.c_str(),ext);
  cout<<"Saving AnaNue Tree to "<<setfilename<<endl;
  j.Path("NueAna").Mod("Output").Set(setfilename);
  
  char pottreefn[1000];
  sprintf(pottreefn,"POTTreeName=%s%s%s.root",
	  path,file.c_str(),ext);
  
  j.Path("NueAna").Mod("PIDEval").Set(pottreefn);

  sprintf(pottreefn,"POTTreeNameIn=%s",
          filename.c_str());

  j.Path("NueAna").Mod("PIDEval").Set(pottreefn);



 
  ////////////////////////////////////////////////////////
  //Run
  //j.Input.Next(74);  
  j.Path("Main").Run();
    
  ////////////////////////////////////////////////////////
  //report
  j.Path("Main").Report();
  
  ////////////////////////////////////////////////////////

  
   
}
