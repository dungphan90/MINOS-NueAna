{

JobC j;

///////////////////////////////////////////////////////
//load necessary libraries
gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 

////////////////////////////////////////////////////////

//set messaging levels
//change these to Debug, Info, Warning, Error to set different levels of output
j.Msg.SetLevel("NuePrint","Info");
//j.Input.AddFile("/raid1/minos/boehm/Data/June3Files/near/AnaNue-n1430201*.root","ana_nue");
//j.Input.AddFile("/data/raid1/minos/msanchez/Minos/Data/VarSelSample/AnaNue-*.root","ana_nue");
j.Input.AddFile("~msanchez/Minos/testing-20050908/AnaNue-*.root","ana_nue");
////////////////////////////////////////////////////////
//set the analysis path
j.Path.Create("Main",
	      "NuePrint::Ana ");
////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=ana_nue"); 
////////////////////////////////////////////////////////
// Setup configuration
j.Path("Main").Mod("NuePrint").Set("printAll=0");
j.Path("Main").Mod("NuePrint").Set("training=0");
j.Path("Main").Mod("NuePrint").Set("outFilename=varselection.pat");
j.Path("Main").Mod("NuePrint").Set("outFormat=SPR");
j.Path("Main").Mod("NuePrint").Set("excludeVars=mctrue,bmon,anainfo,begPlane,end,vtx,fHeader,dcos,trace,srevent.event,srevent.triggerPass,stripsPerPlane");


// Cut configuration
// Cuts on max number of planes for the track 
j.Path("Main").Mod("NuePrint").Set("HiPlaneTrackCut=25");   
// Cut on min total ph (in sigcor or sigmap) for either prong shower or track
j.Path("Main").Mod("NuePrint").Set("PhProngCut=5000");   
// cut on max event energy in Meu 
j.Path("Main").Mod("NuePrint").Set("MeuEnergyCut=150");   
// Cuts on max number of tracklike planes
j.Path("Main").Mod("NuePrint").Set("HiTrackLikeCut=18");
// Cuts on min number of planes for the event
//j.Path("Main").Mod("NuePrint").Set("LoPlaneEventCut=-1");
// Cuts on max and min total energy in Meu
//j.Path("Main").Mod("NuePrint").Set("HiEnergyCut=-1");
//j.Path("Main").Mod("NuePrint").Set("LoEnergyCut=-1");
// Cuts on max and min total energy for the shower in GeV
//j.Path("Main").Mod("NuePrint").Set("HiEnergyShowerCut=-1");
//j.Path("Main").Mod("NuePrint").Set("LoEnergyShowerCut=1.5");
// cuts on min number of strips or planes in count above a ph treshold
//j.Path("Main").Mod("NuePrint").Set("LoPhNStripCut=-1"); 
//j.Path("Main").Mod("NuePrint").Set("LoPhNPlaneCut=-1"); 

//Request oscillation recalculation (if set to 1, default 0)
// This also prints a weight other than 1 for each event 
j.Path("Main").Mod("NuePrint").Set("recalcOsc=0");

// Select a different default string (as an integer for unknown values)

//j.Path("Main").Mod("NuePrint").Set("defaultTrkVal=-1");
//Set alternative oscillation parameters
//j.Path("Main").Mod("NuePrint").Set("DM2=0.002175");
//j.Path("Main").Mod("NuePrint").Set("Theta23=0.5905");
//j.Path("Main").Mod("NuePrint").Set("UE32=0.01");

// Select different samples for training purposes:
// To eliminate one, just change to -1, defaults are as written
//j.Path("Main").Mod("NuePrint").Set("SigClass=2");
//j.Path("Main").Mod("NuePrint").Set("BgNCClass=0");
//j.Path("Main").Mod("NuePrint").Set("BgNumuClass=1");
//j.Path("Main").Mod("NuePrint").Set("BgNutauClass=3");
//j.Path("Main").Mod("NuePrint").Set("BgBNueClass=4");

// Select Resonance Codes less or equal than (default=1004):
//j.Path("Main").Mod("NuePrint").Set("SigResCode=1004");

// Select EmFrac (truth) greater or equal than (default=-1)
//j.Path("Main").Mod("NuePrint").Set("EmFrac=-1");

// note only set this file if you want to print a specific set of vars
// note exclusion precedes inclusion, format as AllParam.txt
// j.Path("Main").Mod("NuePrint").Set("includeFile=myParam.txt");
////////////////////////////////////////////////////////
// Run
j.Path("Main").Run();
//j.Path("Main").RunNpass(50);
//j.Path("Main").RunNin(50);

////////////////////////////////////////////////////////
//report
j.Path("Main").Report();

////////////////////////////////////////////////////////



}
