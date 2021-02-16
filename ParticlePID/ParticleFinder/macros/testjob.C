void testjob(int nSkip=0, int nRun=1){
JobC j;

j.Input.Set("Format=input");
j.Input.Set("Streams=NtpSt");
j.Path.Create("test",
		"ParticleTruthMaker::Reco "
		"ParticleFinder::Reco "
	//	"ParticleReport::Reco "
		"Output::Put "
		);

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
