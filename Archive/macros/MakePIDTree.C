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
j.Path("Main").Mod("FillPIDFromText").Set("TextFile=testpid.txt");
//set threshold for # of planes in track, above which no nue anal will be done.
j.Path("Main").Mod("FillPIDFromText").Set("Decider=99");

//configure selection for decider 6 (Tony/Dan)
// j.Path("Main").Mod("FillPIDFromText").Set("SelFlav=2");
// j.Path("Main").Mod("FillPIDFromText").Set("SelRes=2");


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
string rnumstring = filename.substr(filename.find_last_of("/")+2,8);
string subrunstring = filename.substr(filename.find_last_of("/")+11,4);
char setfilename[5000];
sprintf(setfilename,"FileName=%sAnaNue-%s_%s.root",
	path,rnumstring.c_str(),subrunstring.c_str());
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
