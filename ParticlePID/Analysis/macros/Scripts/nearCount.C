
{

        gSystem->CompileMacro("/data/minos/NueAnalysis/Code/MiniPlotMaker/Stage2.C");
	Stage2 s;
	s.SetInputFile("staged.root");
	s.SetOutputFile("staged.root");
	s.SetOverwrite(1);
	s.SetPlotDirectory("/data/minos/NueAnalysis/Scripts/MiniPlotMaker/plots");


	s.Run();
s.CountComponentsByCutLevelNear();
}
