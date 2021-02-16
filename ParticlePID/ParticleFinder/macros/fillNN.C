void fillNN(int nSkip=0, int nRun=-1){

  gROOT->LoadMacro("NueAna/ParticlePID/ParticleFinder/macros/LoadLibs.C"); LoadLibs();



JobC j;

j.Input.Set("Format=input");
j.Input.Set("Streams=PA,PO,BeamMon");
j.Path.Create("test",
		"PIDEval::Reco "
		"Output::Put "
		);


j.Path("test").Mod("Output").Cmd("DefineStream  PA PRecord");
j.Path("test").Mod("Output").Cmd("DefineStream  PO ParticleObjectHolder");
j.Path("test").Mod("Output").Cmd("DefineStream  BeamMon ParticleBeamMon");
j.Path("test").Mod("Output").Set("Streams=PO,PA,BeamMon");


  string filename;
  
    JobCEnv& jce = JobCEnv::Instance();
    filename = jce.GetFileName(0);

  string rnumstring = filename.substr(filename.find_last_of("/")+2,8);
  string subrunstring = filename.substr(filename.find_last_of("/")+11,4);

  string file = filename.substr(filename.find_last_of("/")+1,
                                filename.find_last_of(".root") -
                                filename.find_last_of("/")-5);

  char setfilename[5000];
  sprintf(setfilename,"FileName=%s-NN.root", file.c_str());
  cout<<"Saving PRecord Tree to "<<setfilename<<endl;
  j.Path("test").Mod("Output").Set(setfilename);

	//need to copy the pot tree!
  sprintf(setfilename,"InFileName=%s", filename.c_str());
  j.Path("test").Mod("PIDEval").Set(setfilename);
  sprintf(setfilename,"OutFileName=%s-NN.root", file.c_str());
  j.Path("test").Mod("PIDEval").Set(setfilename);


j.Input.Next(nSkip);
j.Path("test").Run(nRun);



j.Path("test").Report();
}
