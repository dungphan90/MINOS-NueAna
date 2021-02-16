{
JobC j;

j.Input.Set("Format=input");
j.Input.Set("Streams=NtpSt, PO");
j.Path.Create("test",
		"ParticleFinder::Reco "
		"Output::Put "
		);

j.Path("test").Mod("Output").Set("FileName=test_out_n.root");
j.Path("test").Mod("Output").Cmd("DefineStream  NtpSt NtpStRecord");
j.Path("test").Mod("Output").Cmd("DefineStream  PO ParticleObjectHolder");
j.Path("test").Mod("Output").Set("Streams=NtpSt,PO");

//j.Path("test").Run(10);

j.Path("test").Run();

j.Path("test").Report();
}
