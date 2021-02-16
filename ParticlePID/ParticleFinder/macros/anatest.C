void anatest(int nSkip=0, int nRun=-1){

  gROOT->LoadMacro("NueAna/ParticlePID/ParticleFinder/macros/LoadLibs.C"); LoadLibs();



JobC j;

j.Input.Set("Format=input");
j.Input.Set("Streams=PO");
j.Path.Create("test",
	//	"ParticleTruthMaker::Reco "
	//	"ParticleFinder::Reco "
	//	"ParticleReport::Reco "
		"ParticleAna::Reco "
		"Output::Put "
		);

//j.Path("test").Mod("Output").Set("FileName=test_out.root");
//j.Path("test").Mod("Output").Cmd("DefineStream  NtpSt NtpStRecord");
j.Path("test").Mod("Output").Cmd("DefineStream  PA PRecord");
j.Path("test").Mod("Output").Cmd("DefineStream  PO ParticleObjectHolder");
j.Path("test").Mod("Output").Cmd("DefineStream  BeamMon ParticleBeamMon");
j.Path("test").Mod("Output").Set("Streams=PO,PA,BeamMon");


j.Msg.SetLevel("ParticleFinder","Error");
j.Msg.SetLevel("ParticleAna","Debug");

  string filename;
  
    JobCEnv& jce = JobCEnv::Instance();
    filename = jce.GetFileName(0);

  string rnumstring = filename.substr(filename.find_last_of("/")+2,8);
  string subrunstring = filename.substr(filename.find_last_of("/")+11,4);

  string file = filename.substr(filename.find_last_of("/")+1,
                                filename.find_last_of(".root") -
                                filename.find_last_of("/")-5);

  char setfilename[5000];
  sprintf(setfilename,"FileName=PA-%s.root",
          file.c_str());
  cout<<"Saving PRecord Tree to "<<setfilename<<endl;
  j.Path("test").Mod("Output").Set(setfilename);

  sprintf(setfilename,"POTTreeFileName=PA-%s.root",
          file.c_str());
  j.Path("test").Mod("ParticleAna").Set(setfilename);

//j.Input.Next(34);
//j.Path("test").Run(1);

j.Input.Next(nSkip);
j.Path("test").Run(nRun);

j.Path("test").Report();
}
