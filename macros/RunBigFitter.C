{
  gROOT->LoadMacro("NueAna/macros/LoadLibs.C"); LoadLibs();
  gSystem->Load("libNueAnaExtrapolation.so");

  NueFitter  nf;
  nf.AddFile("Grid_1e4_MinNum_Fixed.root");
  nf.AddFile("Grid_1e4_Set2_MinNum.root");
  nf.AddFile("Grid_1e4_Set3_MinNum.root");
  nf.AddFile("Grid_1e4_Set4_MinNum.root");  

  nf.AddContour(68);
  nf.AddContour(90);
  nf.PerformFit(40, "FinalChi2_Large_40.root");


  NueFitter  nf2;
  nf2.AddFile("Grid_1e4_MinNum_Fixed.root");
  nf2.AddFile("Grid_1e4_Set2_MinNum.root");
  nf2.AddFile("Grid_1e4_Set3_MinNum.root");
  nf2.AddFile("Grid_1e4_Set4_MinNum.root");
                                                                                
  nf2.AddContour(68);
  nf2.AddContour(90);
  nf2.PerformFit(42, "FinalChi2_Large_42.root");


  NueFitter  nf3;
  nf3.AddFile("Grid_1e4_MinNum_Fixed.root");
  nf3.AddFile("Grid_1e4_Set2_MinNum.root");
  nf3.AddFile("Grid_1e4_Set3_MinNum.root");
  nf3.AddFile("Grid_1e4_Set4_MinNum.root");

  nf3.AddContour(68);
  nf3.AddContour(90);
  nf3.PerformFit(46, "FinalChi2_Large_46.root");


  NueFitter  nf4;
  nf4.AddFile("Grid_1e4_MinNum_Fixed.root");
  nf4.AddFile("Grid_1e4_Set2_MinNum.root");
  nf4.AddFile("Grid_1e4_Set3_MinNum.root");
  nf4.AddFile("Grid_1e4_Set4_MinNum.root");

  nf4.AddContour(68);
  nf4.AddContour(90);
  nf4.PerformFit(57, "FinalChi2_Large_57.root");


  NueFitter  nf5;
  nf5.AddFile("Grid_1e4_MinNum_Fixed.root");
  nf5.AddFile("Grid_1e4_Set2_MinNum.root");
  nf5.AddFile("Grid_1e4_Set3_MinNum.root");
  nf5.AddFile("Grid_1e4_Set4_MinNum.root");

  nf5.AddContour(68);
  nf5.AddContour(90);
  nf5.PerformFit(62, "FinalChi2_Large_62.root");

  NueFitter  nf6;
  nf6.AddFile("Grid_1e4_MinNum_Fixed.root");
  nf6.AddFile("Grid_1e4_Set2_MinNum.root");
  nf6.AddFile("Grid_1e4_Set3_MinNum.root");
  nf6.AddFile("Grid_1e4_Set4_MinNum.root");

  nf6.AddContour(68);
  nf6.AddContour(90);
  nf6.PerformFit(51, "FinalChi2_Large_51.root");


}

