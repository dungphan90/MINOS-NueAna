
{

        gSystem->CompileMacro("$SRT_PRIVATE_CONTEXT/NueAna/ParticlePID/Analysis/macros/Stage2.C");
	Stage2 s;
	s.SetInputFile("staged.root");
	s.SetOutputFile("staged.root");
	s.SetOverwrite(1);
	s.SetPlotDirectory("plots");


	s.Run();

}
