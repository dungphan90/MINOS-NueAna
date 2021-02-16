
{

        gSystem->CompileMacro("/data/minos/NueAnalysis/Code/MiniPlotMaker/Stage1.C");

	Stage1 s;
	s.SetInputFile("staged.root");
	s.SetOutputFile("staged.root");
	s.SetOverwrite(1);

	s.Run();

}
