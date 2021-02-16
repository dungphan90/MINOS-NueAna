void chaintest(int nSkip=0, int nRun=-1){

  gROOT->LoadMacro("NueAna/ParticlePID/ParticleFinder/macros/LoadLibs.C"); LoadLibs();


JobC j;

j.Input.Set("Format=input");
j.Input.Set("Streams=NtpSt");
j.Path.Create("test",
		"ParticleFinder::Reco "
	//	"ParticleTruthMaker::Reco "   //run after particle finder to try to match found particles to true... (and can add records for events not recoed)
	//	"ParticleReport::Reco "
		"ParticleAna::Reco "
		"Output::Put "
		);

//j.Path("test").Mod("Output").Set("FileName=test_out.root");
j.Path("test").Mod("Output").Cmd("DefineStream  NtpSt NtpStRecord");
j.Path("test").Mod("Output").Cmd("DefineStream  PO ParticleObjectHolder");
j.Path("test").Mod("Output").Cmd("DefineStream  PA PRecord");
//j.Path("test").Mod("Output").Set("Streams=NtpSt,PO");
j.Path("test").Mod("Output").Set("Streams=PO,PA");

j.Msg.SetLevel("ParticleFinder","Warning");
j.Msg.SetLevel("ChainHelper","Warning");

  string filename;
  
    JobCEnv& jce = JobCEnv::Instance();
    filename = jce.GetFileName(0);

  string rnumstring = filename.substr(filename.find_last_of("/")+2,8);
  string subrunstring = filename.substr(filename.find_last_of("/")+11,4);

  string file = filename.substr(filename.find_last_of("/")+1,
                                filename.find_last_of(".root") -
                                filename.find_last_of("/")-5);

  char setfilename[5000];
  sprintf(setfilename,"FileName=PO-%s.root",
          file.c_str());
  cout<<"Saving AnaNue Tree to "<<setfilename<<endl;
  j.Path("test").Mod("Output").Set(setfilename);


//j.Input.Next(34);
//j.Path("test").Run(1);

j.Input.Next(nSkip);
j.Path("test").Run(nRun);

j.Path("test").Report();
}
