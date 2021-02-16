{

JobC j;

////////////////////////////////////////////////////////
//load necessary libraries
gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 

////////////////////////////////////////////////////////
//set messaging levels
//change these to Debug, Info, Warning, Error to set different levels of output
j.Msg.SetLevel("FillPIDFromText","Info");

////////////////////////////////////////////////////////
//set the analysis path
j.Path.Create("Main",
	      "FillPIDFromText::Reco "
	      "Output::Put");

////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=NtpSR"); 

////////////////////////////////////////////////////////
//configure NueModule
//tj's pid file from mc

//tj's pid file from challenge set

//alexs
j.Path("Main").Mod("FillPIDFromText").Set("TextFile=/afs/fnal.gov/files/home/room3/vahle/MinosSoft/test/NueAna/data/ASPID.txt");

//set decider
j.Path("Main").Mod("FillPIDFromText").Set("Decider=2");

////////////////////////////////////////////////////////
//Configure output

//set the streams that will be written
j.Path("Main").Mod("Output").Cmd("DefineStream pid_nue NuePID");
j.Path("Main").Mod("Output").Set("Streams=pid_nue");

//Get the run number from the first Input file name
//doesn't have to be this way, I just thought I'd be tricky
char *path = getenv("NUE_OUT_DIR");
if(path=='\0'){
   path="./";
}
JobCEnv& jce = JobCEnv::Instance();
string filename = jce.GetFileName(0);
string rnumstring = filename.substr(filename.find_last_of("_")-9,9);
string det=rnumstring.substr(0,1);
int itgt = atoi((rnumstring.substr(3,1)).c_str());
int flv = atoi((rnumstring.substr(4,1)).c_str());
string sflv="";
if(flv==0){ sflv="beam"; }
else if(flv==1){ sflv="nue";}
else if(flv==3){ sflv="nutau";}
else { sflv="ukn";}
int run = atoi((rnumstring.substr(5,4)).c_str());
char setfilename[5000];
//sprintf(setfilename,"FileName=%sPIDNue-%s-%s-%d-%d.root",
//	path,det.c_str(),sflv.c_str(),itgt,run);
sprintf(setfilename,"FileName=%sPIDNue-%s.root",path,rnumstring.c_str());
cout<<"Saving PIDNue Tree to "<<setfilename<<endl;
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
