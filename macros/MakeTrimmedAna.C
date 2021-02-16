{

JobC j;

////////////////////////////////////////////////////////
//load necessary libraries
gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 

////////////////////////////////////////////////////////
//set messaging levels
//change these to Debug, Info, Warning, Error to set different levels of output
//j.Msg.SetLevel("Trimmer","Debug");
//j.Msg.SetLevel("Io","Debug");
//j.Msg.SetLevel("Per","Debug");
////////////////////////////////////////////////////////
//set the analysis path
JobCPath& main_path = j.Path.Create("Main",
	      "TrimModule::Reco "
              "Output::Put");

////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=ana_nue"); 

////////////////////////////////////////////////////////
//Configure output

//j.Input.AddFile("nue/Ana*.root","ana_nue");
//j.Input.AddFile("numu/Ana*.root","ana_nue");
//j.Input.AddFile("nutau/Ana*.root","ana_nue");

JobCModule& nue_mod = main_path.Mod("TrimModule"); 

// set threshold for # of planes in track, above which no nue anal will be done.
// cuts are done <= or >=, energy is sigcor based
 nue_mod.Set("HiPlaneTrackCut=-1"); // max number of planes in track
//nue_mod.Set("HiEnergyCut=150"); // max total energy in Meu
//nue_mod.Set("HiTrackLikeCut=16"); // max number of tracklike planes

//Trim event outside Fiducial and Containment volumes
// nue_mod.Set("FiducialCut=1");
// nue_mod.Set("ContainmentCut=1");
// min shower in Gev, if more than one shoshower cut is done at half the value
//nue_mod.Set("LoEnergyShowerCut=1.05");


//nue_mod.Set("LoPlaneEventCut=-1"); // min number of planes for event
//nue_mod.Set("LoEnergyCut=-1"); // min total energy in Meu
//nue_mod.Set("HiEnergyShowerCut=-1");// max shower energy in GeV

// Cut on min total ph (in sigcor or sigmap) for either prong shower or track
//nue_mod.Set("PhProngCut=5000");   

// cuts on min number of strips or planes in count above a ph treshold
//nue_mod.Set("LoPhNStripCut=14"); 
//nue_mod.Set("LoPhNPlaneCut=-1"); 

// Select different samples for training purposes:
// To eliminate one, just change to -1, defaults are as written
//By default all types will pass through to remove a specific type
nue_mod.Set("CutOnClasses=1");  //default option
nue_mod.Set("AddSignal=0");
nue_mod.Set("AddSignal=1"); 
nue_mod.Set("AddSignal=2");
nue_mod.Set("AddSignal=3");
nue_mod.Set("AddSignal=4");
//Remove the line above for the classes you don't want

// Select Resonance Codes
//nue_mod.Set("SetSignalResCode=1001");
//nue_mod.Set("AddSignalResCode=1002");

// Select EmFrac (truth) greater or equal than (default=-1)
//nue_mod.Mod("NuePrint").Set("HiSigEmFracCut=-1");


//Reweighting for oscillations
//nue_mod.Set("ReWeight=1");
//nue_mod.Set("DeltaMSquare=0.0025");
//nue_mod.Set("Ue3Square=0.01");
//nue_mod.Set("Theta23=TMath::Pi()/4");


//set the streams that will be written
j.Path("Main").Mod("Output").Cmd("DefineStream ana_nue NueRecord NueRecord_trim");
j.Path("Main").Mod("Output").Set("Streams=ana_nue");

//set output name
char *path = getenv("NUE_OUT_DIR");
if(path=='\0'){
   path="./";
}

JobCEnv& jce = JobCEnv::Instance();
string filename = jce.GetFileName(0);
string file = filename.substr(filename.find_last_of("/")+1, filename.find_last_of(".root")-filename.find_last_of("/")-5);

char setfilename[5000];
sprintf(setfilename,"FileName=%s-Trim.root", file.c_str());
cout<<"Saving AnaNue Tree to "<<setfilename<<endl;
j.Path("Main").Mod("Output").Set(setfilename);

//j.Path("Main").Mod("Output").Set("AutoSaveBytes=10000000");



////////////////////////////////////////////////////////
//Run
j.Path("Main").Run();
//j.Path("Main").RunNpass(50);
//j.Path("Main").RunNin(50);

////////////////////////////////////////////////////////
//report
j.Path("Main").Report();

////////////////////////////////////////////////////////

}
