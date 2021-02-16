{

JobC j;

////////////////////////////////////////////////////////
//load necessary libraries
gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs();
 

////////////////////////////////////////////////////////
//set messaging levels
//change these to Debug, Info, Warning, Error to set different levels of output
j.Msg.SetLevel("PETrimmer","Debug");  //Debug, Fatal


////////////////////////////////////////////////////////
//set the analysis path
JobCPath& main_path = j.Path.Create("Main",
	      "PETrimmer::Reco "
              "Output::Put");

////////////////////////////////////////////////////////
//Set Input Streams
j.Input.Set("Streams=NtpSt"); 

////////////////////////////////////////////////////////
//Configure output

JobCModule& nue_mod = main_path.Mod("PETrimmer"); 


float pecut=2.0;

char tmp[100];
sprintf(tmp,"PECut=%f",pecut);

 nue_mod.Set(tmp); // max number of planes in track


nue_mod.Set("updateEventEnergy=1");


//set the streams that will be written
j.Path("Main").Mod("Output").Cmd("DefineStream NtpSt NtpStRecord");
j.Path("Main").Mod("Output").Set("Streams=NtpSt");

//set output name
char *path = getenv("NUE_OUT_DIR");
if(path=='\0'){
   path="./";
}

JobCEnv& jce = JobCEnv::Instance();
string filename = jce.GetFileName(0);
string file = filename.substr(filename.find_last_of("/")+1, filename.find_last_of(".root")-filename.find_last_of("/")-5);

char setfilename[5000];
sprintf(setfilename,"FileName=%s-PECut-%f.root", file.c_str(),pecut);
cout<<"Saving AnaNue Tree to "<<setfilename<<endl;
j.Path("Main").Mod("Output").Set(setfilename);



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
