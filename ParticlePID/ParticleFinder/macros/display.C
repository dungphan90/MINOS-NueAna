void display(int nSkip=0, int nRun=1){


  gROOT->LoadMacro("NueAna/ParticlePID/ParticleFinder/macros/LoadLibs.C"); LoadLibs();



gSystem->Load("libParticleFinder.so");
gSystem->Load("libParticleDisplay.so");
gSystem->Load("libParticleTruth.so");
gSystem->Load("libParticleReporter.so");

JobC j;

j.Input.Set("Format=input");
j.Input.Set("Streams=NtpSt");
j.Path.Create("test",
//		"ParticleTruthMaker::Reco "
		"ParticleFinder::Reco "
		"ParticleDisplay::Reco "
//		"ParticleReport::Reco "
		"Output::Put "
		);
		
j.Msg.SetLevel("ParticleFinder","Debug");

j.Path("test").Mod("Output").Set("FileName=test_out.root");
j.Path("test").Mod("Output").Cmd("DefineStream  NtpSt NtpStRecord");
j.Path("test").Mod("Output").Cmd("DefineStream  PO ParticleObjectHolder");
j.Path("test").Mod("Output").Set("Streams=NtpSt,PO");

//j.Input.Next(34);
//j.Path("test").Run(1);

j.Input.Next(nSkip);
j.Path("test").Run(nRun);

j.Path("test").Report();
}
