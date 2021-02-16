void LoadLibs(){

  gSystem->Load("libMCNtuple.so");
  gSystem->Load("libCandNtupleSR.so");
  gSystem->Load("libCandShield.so");
  gSystem->Load("libTruthHelperNtuple.so");
  gSystem->Load("libStandardNtuple.so");
  gSystem->Load("libCandNtupleEM");
  gSystem->Load("libBeamDataNtuple");
  gSystem->Load("libAnalysisNtuples.so");
  gSystem->Load("libMCReweight.so");
  gSystem->Load("libMad.so");        //needed for DP PID 
  gSystem->Load("libAnalysisNtuplesModule.so");
  gSystem->Load("libNeugenInterface.so");
  gSystem->Load("libBeamDataAnaSummary.so");
  gSystem->Load("libBeamDataUtil.so");
  gSystem->Load("libDataUtil.so");
  gSystem->Load("libOscProbMatrixDecomp.so");
  gSystem->Load("libOscProb.so");
  gSystem->Load("libSpillTiming.so");

  gSystem->Load("libMuonRemoval");   //needed for muon removal info

  gSystem -> Load("libPhysicsNtuple.so");
  gSystem -> Load("libPhysicsNtupleFill.so");
  gSystem -> Load("libPhysicsNtupleSelect.so");
  gSystem -> Load("libPhysicsNtupleStore.so");
  gSystem -> Load("libPhysicsNtuplekNNAlg.so");

  gSystem->Load("libAtNuEvent.so");  //needed for MOISolution code
  gSystem->Load("libMuELoss.so");    //needed for MOISolution code
  gSystem->Load("libRunQuality.so");
  gSystem->Load("libCandMorgue.so"); //needed for MOISolution code
  gSystem->Load("libAtNuOutput.so"); //needed for MOISolution code

  gSystem->Load("libMinuit2");
  gSystem->Load("libNtpFitSA.so");
  gSystem->Load("libNtupleUtils.so");


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




        //needed for steves reconstruction
  gSystem->Load("libMLP.so");
  gSystem->Load("libNueAnaParticleManaged.so");
  gSystem->Load("libNueAnaParticleFinder.so");
  gSystem->Load("libNueAnaTools.so");
  gSystem->Load("libNueAnaParticleAna.so");

        //end steves reconstruction

  gSystem->Load("libNueAna.so");
  gSystem->Load("libNueAnaModule.so");
  gSystem->Load("libXTalkFilter.so");
  gSystem->Load("libMCNNAnalysis.so");
  gSystem->Load("libNueAnaParticlePID.so");
  gSystem->Load("libNueAnaParticlePIDAnalysis.so");

  //Start For LisaFiles
  gSystem->Load("libNueAnaMultiBinAna.so");
  gSystem->Load("libNueAnaExtrapolation.so");
  gSystem->Load("libNueAnaTools.so");
  //End For LisaFiles


}

void LoadDisplayLibs() {
  LoadLibs();
  
  gSystem->Load("libMidadUtil.so");
  gSystem->Load("libMidad.so");
  gSystem->Load("libMidadMultiPage.so");
  
  gSystem->Load("libCandShowerEM.so");    //for predicting EM shower shape
  gSystem->Load("libCandFitShowerEM.so"); //for predicting EM shower shape
  gSystem->Load("libNueAnaDisplay.so");

  gSystem->Load("libNueAnaParticleDisplay.so");


}
