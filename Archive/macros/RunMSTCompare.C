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
j.Msg.SetLevel("CompareMST","Info");
j.Msg.SetLevel("HistMan","Info");

////////////////////////////////////////////////////////
//set the analysis path
j.Path.Create("Main",
	      "CompareMST::Ana ");
////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=ana_nue"); 

////////////////////////////////////////////////////////
//Run
j.Path("Main").Run();
//j.Path("Main").RunNpass(50);
//j.Path("Main").RunNin(5);

////////////////////////////////////////////////////////
//report
j.Path("Main").Report();

////////////////////////////////////////////////////////

//Get the run number from the first Input file name
//doesn't have to be this way, I just thought I'd be tricky
char *path = getenv("NUE_OUT_DIR");
if(path=='\0'){
   path="./";
}
JobCEnv& jce = JobCEnv::Instance();
string filename = jce.GetFileName(0);
string rnumstring = filename.substr(filename.find_last_of("_")-9,9);
string subrunstring = filename.substr(filename.find_last_of("_")+1,4);
char setfilename[5000];
sprintf(setfilename,"%smstcomp-%s_%s.root",
	path,rnumstring.c_str(),subrunstring.c_str());
HistMan *hm = new HistMan("mstcomp");
hm->WriteOut(setfilename);

}
