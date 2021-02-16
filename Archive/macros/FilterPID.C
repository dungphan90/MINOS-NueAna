{

JobC j;

////////////////////////////////////////////////////////
//load necessary libraries

gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 

////////////////////////////////////////////////////////
//set messaging levels
//change these to Debug, Info, Warning, Error to set different levels of output
//j.Msg.SetLevel("FilterPID","Debug");
//j.Msg.SetLevel("FilterPID","Debug");
//j.Msg.SetLevel("Io","Debug");
//j.Msg.SetLevel("Per","Debug");
////////////////////////////////////////////////////////
//set the analysis path
j.Path.Create("Main",
	      "FilterPID::Reco "
              "Output::Put");

////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=ana_nue,pid_nue"); 

////////////////////////////////////////////////////////
//Configure output

//set the streams that will be written
j.Path("Main").Mod("Output").Cmd("DefineStream ana_nue_filt NueRecord NueRecord_filt");
j.Path("Main").Mod("Output").Set("Streams=ana_nue_filt");

//set output name
char *path = getenv("NUE_OUT_DIR");
if(path=='\0'){
   path="./";
}
JobCEnv& jce = JobCEnv::Instance();
string filename = jce.GetFileName(0);
string rnumstring = filename.substr(filename.find_last_of("_")-9,9);
string subrunstring = filename.substr(filename.find_last_of("_")+1,4);
char setfilename[5000];
sprintf(setfilename,"FileName=%sAnaNue-%s_%s.filt.root",
	path,rnumstring.c_str(),subrunstring.c_str());
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
