{
/////NOTE: Run on all nd files at once so that we can compute the normalization correctly
//if NueReweightFast finds no events in a filterd ana_nue tree, it will not write out 
//an file, then it makes it difficult to figure out the proper normalization!!!!



////////////////////////////////////////////////////////
//load necessary libraries
gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 
gSystem->Load("libNueAnaReweight.so");

//gSystem->Load("/home/vahle/MinosSoft/test/NueAna/macros/NueReweightFast_C.so");
gSystem->Load("/afs/fnal.gov/files/home/room3/vahle/MinosSoft/test/NueAna/macros/NueReweightFast_C.so");
//JobC j;
////////////////////////////////////////////////////////
//Get the Input file names
char *path = getenv("RW_OUT_DIR");
if(path=='\0'){
  path="./";
}

JobCEnv& jce = JobCEnv::Instance();
string filename = jce.GetFileName(0);
 
string rnumstring = filename.substr(filename.find_last_of("_")-9,9);
string subrunstring = filename.substr(filename.find_last_of("_")+1,4);
char setfilename[5000];
sprintf(setfilename,"%s/rwtree-%s_%s.root",path,rnumstring.c_str(),subrunstring.c_str());
cout<<"Saving rwtree to "<<setfilename<<endl;
 
TChain nue("ana_nue_filt","ana_nue_filt");
for(int i=0;i<jce.GetNfile();i++){
   cout<<"Adding "<<jce.GetFileName(i)<<endl;
   nue.Add(jce.GetFileName(i));
}

char *rpath = getenv("SRT_PRIVATE_CONTEXT");
if(rpath=='\0'){
  rpath=getenv("SRT_PUBLIC_CONTEXT");
  if(rpath=='\0'){
    rpath="./";
  }
}
char rndfile[5000];
//sprintf(rndfile,"%s/NueAna/data/biguniformnumbers.txt",rpath);
//sprintf(rndfile,"%s/NueAna/data/biguniformnumbers-specialforchris.txt",rpath);
sprintf(rndfile,"%s/NueAna/data/biguniformnumbers-specialforchrisFD.txt",rpath);


cout<<"Starting reweight at: "<<endl;
system("date");
NueReweightFast(rndfile,&nue,setfilename);
cout<<"Ending reweight at: "<<endl;
system("date");


}
