{

JobC j;

////////////////////////////////////////////////////////
//load necessary libraries
gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs(); 
gSystem->Load("libNueAnaReweight.so");

////////////////////////////////////////////////////////
//set messaging levels
//change these to Debug, Info, Warning, Error to set different levels of output
//j.Msg.SetLevel("NueReweight","Debug");
j.Msg.SetLevel("NueReweight","Info");
//j.Msg.SetLevel("NueReweight","Fatal");
j.Msg.SetLevel("NeugenWC","Warning");
//  j.Msg.SetLevel("MCReweight","Debug");
////////////////////////////////////////////////////////
j.Path.Create("Main",
	      "NueReweight::Reco ");   

j.Input.Set("Streams=ana_nue_filt"); 

////////////////////////////////////////////////////////
//set filename
char *path = getenv("RW_OUT_DIR");
if(path=='\0'){
  path="./";
}
JobCEnv& jce = JobCEnv::Instance();
string filename = jce.GetFileName(0);
string rnumstring = filename.substr(filename.find_last_of("_")-9,9);
//string newfilename = spath+"/rwtest-"+rnumstring+".root";
char newfilename[50000];
sprintf(newfilename,"%s/rwtest-%s.root",
	path,rnumstring.c_str());
cout<<"Saving file to "<<newfilename<<endl;

//figure out what kind of run we're looking at
//can't do it in NueReweight, because there may be no passed events
//in the file, so we'll never get a header
int run=atoi((filename.substr(filename.find_last_of("_")-8,8)).c_str());
int subrun=atoi((filename.substr(filename.find_last_of("_")+1,4)).c_str());
DetectorType::Detector_t det;
NueRW::FileType_t ftype;
string d=filename.substr(filename.find_last_of("_")-9,1);
if(d=="n"||d=="N"){
  det=DetectorType::kNear;
  ftype=NueRW::kBEAM;
}
else if(d=="f"||d=="F"){
  det=DetectorType::kFar;
  if(run>=24100000&&run<24110000){
    ftype=NueRW::kBEAM;
  }
  else if(run>=24110000&&run<24120000){
    ftype=NueRW::kNUE;
  }
  else if(run>=24130000&&run<24140000){
    ftype=NueRW::kTAU;
  }
  else{
    ftype=NueRW::kUnknown;
  }
}
else{
  det=DetectorType::kUnknown;
}

cout<<"Run: "<<run<<" subrun "<<subrun
    <<" detector "<<DetectorType::AsString(det)
    <<" ftype "<<NueRW::AsString(ftype)<<endl;

TFile *f = new TFile(newfilename,"RECREATE");
TTree *tree = new TTree("rwtree","rwtree");
NueRW *nrw = new NueRW();
tree->Branch("NueRW","NueRW",&nrw,64000,99);
delete nrw;
nrw=0;

////////////////////////////////////////////////////////
MomNavigator mom=j.Mom;
JobCPath jcp = j.Path("Main");

////////////////////////////////////////////////////////
//Run

system("date");
//  for(int i=0;i<500;i++){
//for(int i=0;i<34;i++){
for(int i=0;i<2675;i++){
  //   system("date");
  char rdn[100];
  sprintf(rdn,"RandNumber=%d",i);
  j.Path("Main").Mod("NueReweight").Set(rdn);
  if(i%100==0){
    cout<<"*******************************"<<rdn<<endl;
  }
  
  //   j.Path("Main").RunNin(10);
  j.Path("Main").Run();
  
  //   cout<<"ran main path"<<endl;
  NueReweight *R = dynamic_cast<NueReweight *>(jcp.GetModule("NueReweight"));
  nrw = R->rw;
  nrw->fRun=run;
  nrw->fSubRun=subrun;
  nrw->fDet=det;
  nrw->fFileType=ftype;
  nrw->randrow=i;
 //   nrw->Print();
  tree->Fill();
  
  R->Reset();
  j.Input.GoToFile(filename.c_str());
  j.Input.PrevFile();
  
  //    cout<<"Starting next loop"<<endl;
}
system("date");
////////////////////////////////////////////////////////
f->Write();
f->Close();

////////////////////////////////////////////////////////
//report
j.Path("Main").Report();

////////////////////////////////////////////////////////
}
