{

JobC j;

////////////////////////////////////////////////////////
//load necessary libraries
gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 
gSystem->Load("libHistMan.so");

////////////////////////////////////////////////////////
//set messaging levels
//change these to Debug, Info, Warning, Error to set different levels of output
//j.Msg.SetLevel("NueReadwPID","Debug");
j.Msg.SetLevel("CompareAll","Info");
j.Msg.SetLevel("HistMan","Info");

//j.Input.AddFile("/data/raid1/minos/boehm/Data/AnaNue-24110001_0000.root", "ana_nue");
//j.Input.AddFile("/data/raid1/minos/boehm/Data/AnaNue-24100001_0000.root", "ana_nue");
//j.Input.AddFile("/data/raid1/minos/boehm/Data/AnaNue-24130001_0000.root", "ana_nue");


//j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/boehm/NueAna/20050921/numu/AnaNue-f24100001_0000.root", "ana_nue");
//j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/boehm/NueAna/20050921/nutau/AnaNue-f24130001_0000.root", "ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/numu/AnaNue-f24100001_0000.root", "ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/nutau/AnaNue-f24130001_0000.root", "ana_nue");
// j.Input.AddFile("/afs/fnal.gov/files/data/minos/d85/AnaNue/R1.16/MC/near/AnaNue-n14301001_0000.root","ana_nue");
 j.Input.AddFile("nf/*.root","ana_nue");


//Get the run number from the first Input file name
//doesn't have to be this way, I just thought I'd be tricky
char *path = getenv("NUE_OUT_DIR");
if(path=='\0'){
   path="./";
}

JobCEnv& jce = JobCEnv::Instance();
string filename = jce.GetFileName(0);
string rnumstring = filename.substr(filename.find_last_of("-")+1,8);
string subrunstring = filename.substr(filename.find_last_of("-")+10,4);


char setfilename[5000];
sprintf(setfilename,"%sHistallcomp-%s_%s.root",
        path,rnumstring.c_str(),subrunstring.c_str());
char outcommand[5000];
sprintf(outcommand, "OutputFile=%s", setfilename);


////////////////////////////////////////////////////////
//set the analysis path
j.Path.Create("Main",
	      "CompareAll::Ana "
//              "Output::Put"
	     );
////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=ana_nue"); 
//j.Path("Main").Mod("Output").Set(setfilename);
j.Path("Main").Mod("CompareAll").Set(outcommand);

////////////////////////////////////////////////////////
//Run
j.Path("Main").Run();
//j.Path("Main").RunNpass(50);
//j.Path("Main").RunNin(5);

////////////////////////////////////////////////////////
//report
j.Path("Main").Report();

////////////////////////////////////////////////////////

//HistMan *hm = new HistMan("allcomp");
//hm->WriteOut(setfilename);


}

