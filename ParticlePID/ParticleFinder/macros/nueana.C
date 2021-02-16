void nueana(int nSkip=0, int nRun=-1){

  gROOT->LoadMacro("NueAna/ParticlePID/ParticleFinder/macros/LoadLibs.C"); LoadLibs();


JobC j;

j.Input.Set("Format=input");
j.Input.Set("Streams=NtpSt");
j.Path.Create("test",
		"ParticleFinder::Reco "
		"ParticleAna::Reco "
	//	"ParticleTruthMaker::Reco "   //run after particle finder to try to match found particles to true... (and can add records for events not recoed)
	//	"ParticleReport::Reco "

				      "NueModule::Reco "
	//			      "NueBeamMon::Reco "
		"Output::Put "
		);

  ////////////////////////////////////////////////////////
  //Configure output
  
  //set the streams that will be written
//  j.Path("Main").Mod("Output").Cmd("DefineStream ana_nue NueRecord");
 // j.Path("Main").Mod("Output").Set("Streams=ana_nue");

//j.Path("test").Mod("Output").Set("FileName=test_out.root");
//j.Path("test").Mod("Output").Cmd("DefineStream  NtpSt NtpStRecord");
j.Path("test").Mod("Output").Cmd("DefineStream  PO ParticleObjectHolder");
j.Path("test").Mod("Output").Cmd("DefineStream  PA PRecord");
j.Path("test").Mod("Output").Cmd("DefineStream  BeamMon ParticleBeamMon");
j.Path("test").Mod("Output").Cmd("DefineStream ana_nue NueRecord");
j.Path("test").Mod("Output").Set("Streams=PA,BeamMon,PO,ana_nue");
//j.Path("test").Mod("Output").Set("Streams=PO");

  ////////////////////////////////////////////////////////
  //configure NueModule
JobCModule& nue_mod = j.Path("test").Mod("NueModule");

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


  string filename;
  
    JobCEnv& jce = JobCEnv::Instance();
    filename = jce.GetFileName(0);

  string rnumstring = filename.substr(filename.find_last_of("/")+2,8);
  string subrunstring = filename.substr(filename.find_last_of("/")+11,4);

  string file = filename.substr(filename.find_last_of("/")+1,
                                filename.find_last_of(".root") -
                                filename.find_last_of("/")-5);

  char setfilename[5000];
  sprintf(setfilename,"FileName=PO-%s.root",
          file.c_str());
  cout<<"Saving AnaNue Tree to "<<setfilename<<endl;
  j.Path("test").Mod("Output").Set(setfilename);

  sprintf(setfilename,"POTTreeFileName=PA-%s.root",
          file.c_str());
  j.Path("test").Mod("ParticleAna").Set(setfilename);
  
  
  
  ////Configure Flux Reweighting objects
  //figure out beam type from filename
  int beam=2;
  
  //figure out which flux file directory to use from beam
  //set fixmuparents to non zero if you want to use the mupi trees to change the
  //tp info for muons.  set it to zero if you dont.
  j.Path("test").Mod("NueModule").Set("FixMuParents=0");
  //string ffbasepath="MuPiDir=/stage/minos-data6/vahle/MuParTrees/fluka05_"+flxsfx;
  //string ffbasepath="MuPiDir=/stage/minos-data6/cbs/gnumi/v18/muonInfo/fluka05_"+flxsfx;
  string ffbasepath = "MuPiDir=/afs/fnal.gov/files/data/minos/d19/muonInfo/v18/fluka05_";
  
  cout<<"ffbasepath: "<<ffbasepath<<endl;
  j.Path("test").Mod("NueModule").Set(ffbasepath.c_str());
  

    j.Path("test").Mod("NueModule").Set("ContPhPlaneCut=0.3"); //contiguous plane ph threshold

//j.Input.Next(34);
//j.Path("test").Run(1);

j.Input.Next(nSkip);
j.Path("test").Run(nRun);

j.Path("test").Report();
}
