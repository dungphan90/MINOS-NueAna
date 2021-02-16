void LoadLibs(){




  gSystem->Load("libMCNtuple.so");
  gSystem->Load("libCandNtupleSR.so");
  gSystem->Load("libCandShield.so");
  gSystem->Load("libTruthHelperNtuple.so");
  gSystem->Load("libStandardNtuple.so");
  gSystem->Load("libCandNtupleEM");
  gSystem->Load("libBeamDataNtuple");
  gSystem->Load("libAnalysisNtuples.so");
  gSystem->Load("libAnalysisNtuplesModule.so");
//  gSystem->Load("libNeugenInterface.so");
  gSystem->Load("libMCReweight.so");
  gSystem->Load("libBeamDataAnaSummary.so");
  gSystem->Load("libBeamDataUtil.so");
  gSystem->Load("libDataUtil.so");
  gSystem->Load("libSpillTiming.so");

  gSystem -> Load("libPhysicsNtuple.so");
  gSystem -> Load("libPhysicsNtupleFill.so");
  gSystem -> Load("libPhysicsNtupleSelect.so");
  gSystem -> Load("libPhysicsNtupleStore.so");
  gSystem -> Load("libPhysicsNtuplekNNAlg.so");
  gSystem->Load("libMad.so");        //needed for DP PID 

  gSystem->Load("libMuonRemoval");   //needed for muon removal info

  gSystem->Load("libAtNuEvent.so");  //needed for MOISolution code
  gSystem->Load("libMuELoss.so");    //needed for MOISolution code
  gSystem->Load("libRunQuality.so");
  gSystem->Load("libCandMorgue.so"); //needed for MOISolution code
  gSystem->Load("libAtNuOutput.so"); //needed for MOISolution code

  gSystem->Load("libSwimmer.so");
  gSystem->Load("libCandTrackSR.so");
  gSystem->Load("libCandFitTrackSR.so");
  gSystem->Load("libCandSubShowerSR.so");
  gSystem->Load("libCandShowerSR.so");
  gSystem->Load("libCalDetSI.so");
  gSystem->Load("libCalDetPID.so");
  gSystem->Load("libCalDetDST.so");

  gSystem->Load("libVertexFinder");
  gSystem->Load("libNtpVtxFinder");
  gSystem->Load("libNueAnaTools.so");
  gSystem->Load("libNueAna.so");
  gSystem->Load("libNueAnaModule.so");
  

  gSystem->Load("libNueAnaParticleManaged.so");

  gSystem->Load("libNueAnaParticleFinder.so");
 


 
     gSystem->Load("$ROOTSYS/lib/libMLP.so");

	gSystem->Load("libNueAnaParticleAna.so");



	gSystem->Load("libNueAnaParticleSystematic.so");

}


void LoadDisplayLibs() {
  LoadLibs();
  
  gSystem->Load("libMidadUtil.so");
  gSystem->Load("libMidad.so");
  gSystem->Load("libMidadMultiPage.so");
  
  gSystem->Load("libCandShowerEM.so");    //for predicting EM shower shape
  gSystem->Load("libCandFitShowerEM.so"); //for predicting EM shower shape
//  gSystem->Load("libNueAnaDisplay.so");


gSystem->Load("libParticleDisplay.so");

}
