/////////////////////////////////////////////////////////////////////////////////////////
//
// Description:
//  This macro is suppose to serve as an example of how to run over AnaNue files.
//
// It tries to encapsulate many common activities to avoid rewriting chunks of code
// (Hopefully the encapsulation doesn't add confusion)
//
// Some details about AnaNue files:
//  AnaNue files are ROOT files that are made by the NUE Working group to perform an analysis.
//  AnaNue files are referred to a PAN (Physics Analysis Ntuple) files.  In principle, each 
//  analysis group has their own set of PAN files which contain variables that are useful for
//  that analysis.  AnaNue files contain two trees:
//    1) pottree: The pottree stores the total number of Protons On Target (POT) that correspond
//       to the collection of events stored in the  AnaNue file.  This is used to normalize the 
//       number of events we observe to a common reference.  The pottree has a branch of NuePOT
//       objects.
//    2) ana_nue: The ana_nue tree has an entry for every event (eg neutrino interaction).  The
//       ana_nue tree has a branch of NueRecord objects.  The NueRecord stores all the information
//       for the event that the NUE group thought might be necessary to know for an analysis.
//
//  AnaNue files are generated from SNTP files using the NueAna software package.  
//  SNTP files are the standard Ntuple files that are created by the MINOS reconstruction code.
//  NueAna is the software package created and managed by the NUE Working group that is needed
//  to both create and read AnaNue files.
//
// Some details about this ROOT macro:
//  This ROOT macro is actually a ?loon? macro.  You can think of loon as a super version of root
//  that understands the MINOS software environment.  A basic rule of thumb: if you are going to 
//  use the minos software packages, then you should use loon instead of root.  In principle, you
//  should be able to do anything with loon that you can do with root.
//
//  This macros is suppose to serve as an example of how to run over AnaNue files.  However, instead
//  of going through all the low level details of chaining AnaNue files together and looping over 
//  entries in a tree, it tries to encapsulate many of those common activities through the 
//  AnaNueProcessor class so that the user does not need rewrite these chunks of code over and over
//  again.  Hence, this macro demonstrates how to use the AnaNueProcessor which encapsulates many of
//  the common bits of code that are found in many of the nue analysis macros which handle the opening
//  and reading of AnaNue files.  Note, it is not necessary to use an AnaNueProcessor object to process
//  AnaNue files.  You could just create a TTree, set the branch address to a NueRecord, and start looping
//  over Entries yourself.  Do whatever you find most comfortable.  
//
//  The point of this macro is to make histograms of various quantities of interest with the final intent
//  of producing some specific set of measurements or plots.  I personally have found that it is more 
//  efficient to do most analyses in two phases.  In the first phase, a macro is used to run over the events
//  in a AnaNue file.  In this first phase, histograms are created of various quantities and are written out
//  to a ROOT file.  In the second phase, a second ROOT macro is used to open the file with histograms and 
//  manipulate them to produce various measurements or plots.  This macro includes some histograms that are
//  typically used in analysis.
//
//  When doing an analysis, I have found that I typically make the same histograms ~50 times, but change just
//  a single condition for filling it.  For example, I'll make a plot of the true neutrino energy, but I'll 
//  do it for numu CC events, NC events, BeamNue CC events, etc and I'll make those exact same histograms after
//  various cut values.  In order to reduce the amount of coding that I needed to do, I decided to encapsulate
//  the histogram creation, filling, and writing process, in order to streamline the production of these  
//  similarly structured sets and subsets of histograms.  This is accomplished through the Histograms class
//  which is defined in this macro and inherits from the NueUtilities::Histograms.  This object will write the
//  histograms out to a ROOT file and organize the histograms into folders.  Note, this might produce more 
//  histograms than you will inevitably be interested in, but I typically find histograms to be  objects that
//  carry little cost.  It is usually better to have a variety of histograms on hand during the second phase of
//  analysis that you will ignore, than to find out that you are missing a histogram that you may need later.    
//
//
// Task:  
//       Loop over events and make histograms for:
//       CC, NC, BeamNue, CC+NC+BeamNue, Tau, and SignalNue events
//       Note for the some histograms will not be filled and written out depending
//       on the input files used.  For instance, the ND input files will not produce
//       SignalNue histograms as output
//
// Inputs:
//        The name of an ASCII file containing a list of AnaNue files
//        The file names are delimited by a newline character, '\n'.
//        The files can be Data or MC for either the Far or Near detector.
//        However, the files in the list cannot be of mixed type. That means:
//        They only contain Run1 ND Data
//        They only contain Run2 ND Data
//        They only contain Run3 ND Data
//
//        They only contain Run1 ND MC
//        They only contain Run2 ND MC
//        They only contain Run3 ND MC
//
//        They only contain Run1 FD AnaNue-f214* (nue signal) MC files 
//        They only contain Run1 FD AnaNue-f213* (tau background) MC files 
//        They only contain Run1 FD AnaNue-f210* (beam: NC, CC, and beam nue backgrounds) MC files 
//
//        They only contain Run2 FD AnaNue-f214* (nue signal) MC files 
//        They only contain Run2 FD AnaNue-f213* (tau background) MC files 
//        They only contain Run2 FD AnaNue-f210* (beam: NC, CC, and beam nue backgrounds) MC files 
//
//        They only contain Run3 FD AnaNue-f214* (nue signal) MC files 
//        They only contain Run3 FD AnaNue-f213* (tau background) MC files 
//        They only contain Run3 FD AnaNue-f210* (beam: NC, CC, and beam nue backgrounds) MC files 
//
//        They only contain FD data (FD data is special, you can run over all Runs at once)
//
// Outputs:
//        A ROOT file containing the histograms produced by this macro
//
// How to run it:
//      #Setup the nue analysis software
//        <YourPrompt> setup_nue            
//        <YourPrompt> cd <YourWorkingArea> 
//        <YourPrompt> cp $SRT_PRIVATE_CONTEXT/NueAna/macros/nue_analysis_example.C ./ 
//      #Make any edits you want with whatever editor you are most comfortable with
//        <YourPrompt> xemacs nue_analysis_example.C
//      #Compile the LOON macro
//        <YourPrompt> CompileLoonMacro nue_analysis_example.C 
//      #Make/Get a file list of the AnaNue files to analyze
//      #How about making a file list of Run1 ND Data
//        <YourPrompt> find /minos/data/analysis/nue/2ndAnalysis/AnaNue_Files/BeforeLEM/Full/Trimmed_Fid/Near/Data/L010185N/Standard/Run1/ -name 'AnaNue-*.root' > near_data_run1.dat
//      #Run the macro
//        <YourPrompt> loon -bq 'nue_analysis_example.C+("near_data_run1.dat","near_data_run1.root")' >& log_near_data_run1.txt
//      #If you want the job to run in the background, then instead of doing the above, you should do:
//      #<YourPrompt> nohup loon -bq 'nue_analysis_example.C+("near_data_run1.dat","near_data_run1.root")' >& log_near_data_run1.txt &
//
//      #Below is a recap (without the annoying comments) of how you could run this macro
//        <YourPrompt> setup_nue            
//        <YourPrompt> cd <YourWorkingArea> 
//        <YourPrompt> cp $SRT_PRIVATE_CONTEXT/NueAna/macros/nue_analysis_example.C ./ 
//        <YourPrompt> xemacs nue_analysis_example.C
//        <YourPrompt> /minos/app/nue/Scripts/Utilities/CompileLoonMacro nue_analysis_example.C 
//        <YourPrompt> find /minos/data/analysis/nue/2ndAnalysis/AnaNue_Files/BeforeLEM/Full/Trimmed_Fid/Near/Data/L010185N/Standard/Run1/ -name 'AnaNue-*.root' > near_data_run1.dat
//        <YourPrompt> loon -bq 'nue_analysis_example.C+("near_data_run1.dat","near_data_run1.root")' >& log_near_data_run1.txt
//
//
// Author: Gregory Pawloski 
// Created: August, 2010
////////////////////////////////////////////////////////////////////////








//------------------------------------------------------------------------------
//include minossoft related headers
#include "NueAna/NueUtilities.h"  //Header contains functions that are commonly needed
                                  //in macros to process AnaNue files.  For example:
                                  //   chaining files and getting a POT counting
                                  //Contains AnaNueProcessor which tries to encapsulate
                                  //these functions
                                  //
                                  //Header also happens to contains other headers that
                                  //are needed for many of the common NueAna functions

#include "NueAna/NueAnaTools/NueConvention.h" //Header contains function to recalculate E-scale

//------------------------------------------------------------------------------


using namespace std;  //so you don't have to put std:: in front of things like cout


//------------------------------------------------------------------------------
//
//class: Histograms
//
//A class that is used to handle histograms which includes:
//  Opening a ROOT file to store them in
//  Initializing whatever histograms that you are interested in
//  Filling the histograms
//  Writing the histograms
//  Closing the file 
//
//This is an attempt to make the main function more readable
//Also attempt to avoid rewriting lines of code for multiple histograms that are very similar
//Hopefully avoids copy and paste errors in that case
//
//In this object you make sets of histograms which correspond to events passing a cut.
//For events passing this cut, you then have a subset of histograms based on the event type
//Such as whether it is CC or NC
//In a subset, you then have specific histograms for a specified quantity
//The exact histograms are defined at the end of the macro
//
//Note don't make copies of an object of this class (why would you want to)
//Only pass by reference if you find that you have too
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// It is unlikely that you will need to change anything from here to ....
//------------------------------------------------------------------------------

class Histograms : public NueUtilities::Histograms
{
 public:
  Histograms(){};
  ~Histograms(){};
  Histograms(TString outputFile, TString option){ Open(outputFile, option); };

  inline void Open(TString outputFile, TString option){ OpenFile(outputFile, option); };
  void Fill(TString setName, NueUtilities::AnaNueProcessor & anp, double eventWeight);
  void Write(NueUtilities::AnaNueProcessor & anp);
  inline void Close(){ CloseFile(); };


 private:
  //Don't make a copy of this class! Just pass by reference if needed!
  //Here the assignment and copy constructors are made private to prevent you from doing that
  Histograms(const Histograms&) : NueUtilities::Histograms() {};
  void operator=(const Histograms&);
  void SubsetTask(TString action, TString setName, TString subSetName,
                  void (Histograms::*pHistosInSubset)(TString, TString, TString, NueUtilities::AnaNueProcessor&, double),
                  NueUtilities::AnaNueProcessor & anp, bool passesRuleToFill, bool passesRuleToWrite, double eventWeight);


  void HistoSubsets(TString action, TString setName, NueUtilities::AnaNueProcessor & anp, double eventWeight=0);
  void HistosInSubset(TString action, TString setName, TString subSetName, NueUtilities::AnaNueProcessor & anp, double eventWeight=0);
};


void Histograms::Fill(TString setName, NueUtilities::AnaNueProcessor & anp, double eventWeight)
{
  if( !DoesSetExist(setName) ) HistoSubsets("make", setName, anp, eventWeight);
  HistoSubsets("fill", setName, anp, eventWeight); 
}

void Histograms::Write(NueUtilities::AnaNueProcessor & anp)
{
  for( unsigned int set=0; set < NumSets(); ++set)
  {
   cout << "\n\n Writing Histogram Set: " << GetSetName(set) << "\n"<< endl;
   HistoSubsets("write", GetSetName(set), anp); 
  }
}

void Histograms::SubsetTask(TString action, TString setName, TString subSetName,
                             void (Histograms::* pHistosInSubset)(TString, TString, TString, NueUtilities::AnaNueProcessor&, double),
                             NueUtilities::AnaNueProcessor & anp, bool passesRuleToFill, bool passesRuleToWrite, double eventWeight)
{
  if(action == "make")                            (this->*pHistosInSubset) (action, setName, subSetName, anp, eventWeight);
  else if(action == "fill" && passesRuleToFill)   (this->*pHistosInSubset) (action, setName, subSetName, anp, eventWeight);
  else if(action == "write" && passesRuleToWrite) (this->*pHistosInSubset) (action, setName, subSetName, anp, eventWeight);
  else if(action != "make" && action != "fill" && action != "write")
  {
   cout << "\n\nERROR: Unknown SubsetTask = " << action << "\n"<<endl;
   exit(1);
  }
}

//------------------------------------------------------------------------------
//  .... to here
//------------------------------------------------------------------------------






//------------------------------------------------------------------------------
//The main function which should have the same name as the macro without the .C
// inputFileList: The name of an ASCII file containing a list of AnaNue files
//                The file names are delimited by a newline character ('\n')
//
// outputFile: The name of the output ROOT file to store your histograms, etc.
//------------------------------------------------------------------------------
void nue_analysis_example(string inputFileList, TString outputFile) 
{

 //
 //Step 1:  Open the output ROOT file to store histograms or whatever
 //         Initialize the histograms that you are interested in here
 //         This will vary from task to task
 //        
 //         You could do:
 //          TFile file(outputFile,"RECREATE");
 //          TH1D* h = new TH1D("h","",numbins,lowVal,highVal);
 //
 //         But I encapsulate it all with the Histograms object.
  //

  Histograms histos(outputFile,"RECREATE");


 //
 //Step 2: Perform one time configurations before looping over events
 //

  //Setup oscillation parameters for Far Detector MC
  //Comment these out to get no oscillations
   NueStandard::SetDefaultOscParam();       //Choose this for most studies 
  //NueStandard::SetDefaultOscParamNoNue();  //Choose this for FD background only predictions, FD sidebands


 //
 //Step 3: Loop over events in the input AnaNue files
 //
  NueUtilities::AnaNueProcessor anp(inputFileList);  //Create chain of ana_nue trees, pottrees in the files
                                                     //There is too much overhead with copying this object
                                                     //only pass the anp object by reference

  //Print out the amount of POT that will be processed
  cout << "\n\n Processing " << anp.GetTotalPOT()/100000000.0 << " x 10^{20} POT\n" << endl;

  //loop over events
  while( anp.NextEntry() ) //loop over events in ana_nue chain
  {
   //Event information is stored in pointer to a NueRecord that is accessible through anp.nueRecord
   
    anp.PrintProgress(5); //cout progress after increments of 5% of the entries have been analyzed

    
   //
   //Step 3a: Apply Corrections
   //
    //Recalibrate E-scale before you do anything
    NueConvention::NueEnergyCorrection(anp.nueRecord);
 
    //Set event weight 
    double eventWeight = 1.0;
   
    //Apply MC corrections for flux, etc.
    eventWeight *= NueStandard::GetMCWeights(anp.nueRecord);
   
    //normalize POT
      //This normalization is only appropriate if Run 1, 2, and 3 are in separate samples
    eventWeight *= anp.Normalization_RunSeparatedSample();
      
    //Apply oscillations for Far
    //This method of applying oscillations depends on how you do add your samples together
    //This method works if you keep the f210*, f213*, and f214* files separate
    double oscWeight = 1.0;
    //Don't w0orry about oscillations for now
    oscWeight *= anp.GetOscWeight_f210f213f214Separate();


   //
   //Step 3a: Apply Cuts
   //
    //Apply Fiducial Cuts (which includes Data Quality)
    //I use a "continue" because if it fails this cut, then it will fail all others after it
    if(!NueStandard::PassesFarDataTiming(anp.nueRecord)) continue;  //event is within spill window
    if(!NueStandard::PassesSelection(anp.nueRecord, Selection::kFid)) continue;
    //Fill histograms that pass fiducial cuts here
    histos.Fill("AfterFiducialCut", anp, eventWeight*oscWeight);

    //Apply Preselection Cuts (which includes Fiducial and DQ cuts)
    //I use a "continue" because if it fails this cut, then it will fail all others after it
    if(!NueStandard::PassesSelection(anp.nueRecord, Selection::kPre)) continue;
    //Fill histograms that pass preselection cuts here
    histos.Fill("AfterPreselectionCut", anp, eventWeight*oscWeight);

    //Apply PID Cuts (which includes Preselection, Fiducial and DQ cuts)
    //Note I don't use a "continue", becasue I might be interested in looking at other PIDs later
    //and an event that fails one PID might pass another
    if(NueStandard::PassesSelection(anp.nueRecord, Selection::kANN2PE_DAIKON04))
    {
     //Fill histograms that pass ann11 cut here
     histos.Fill("AfterAnn11Cut", anp, eventWeight*oscWeight);
    }

   //Note, if you are using the Histograms class, you can create a new set of histograms
   //for any cut you just made by just calling:
   //  histos.Fill("YourCut", anp, eventWeight*oscWeight);
   //This call will handle both the initialization, filling, and eventual writing of the histograms


  }//End loop over entries




 //
 //Step 4: Post event loop processing
 //

  //Usually this step is better to do in a second ROOT macros that only runs over the histograms that were produced in this macro
  //If the histograms produced in this macro are numerous, varied, and general enough, then the method of using a second macro can help 
  //avoid looping over events again


 //
 //Step 5: Write histograms to the output file
 //

 histos.Write(anp);

 histos.Close();
}






//------------------------------------------------------------------------------
//Define histograms in Histograms class
//
//In this object you make sets of histograms which correspond to events passing a cut.
//For events passing this cut, you then have a subset of histograms based on the event type
//Such as whether it is CC or NC
//In a subset, you then have specific histograms for a specified quantity
//
//The below methods are used to define subsets that are orgainized by type
//and to define the specific histograms in a subset
//
//------------------------------------------------------------------------------




//------------------------------------------------------------------------------
//
// Defintion for Histograms::HistoSubsets
//
// In this method, you define the histogram subsets
// To define a subset you must give:
//   1)a name for the subset
//   2)a collection of histograms that are in the subset
//   3)the condition for an event to be in the subset
//   4)the condition to write out the subset
//
//
//------------------------------------------------------------------------------
void Histograms::HistoSubsets(TString action, TString setName, NueUtilities::AnaNueProcessor & anp, double eventWeight)
{
// New histograms subsets are created by using the SubsetTask method
// For each new subset you want, you need to call SubsetTask
// This is how you use SubsetTask:
//
// SubsetTask(
//            action,setName,   //DO NOT change this line, always pass these variables
//            <name for the subset>,   //change this line, specify the set name (for example: "Nue")
//            &Histograms::HistosInSubset,  //You probabably won't change this line
//                                          //You have the option of defining a different collection of histograms
//                                          //for this subset by defining a new method like HistosInSubset.
//            anp,   //DO NOT change this line, always pass this variable
//            <condition to fill subset>,  //change this line, specify the condition that must be satisfied for the subset to be filled
//                                         //(for example: anp.isSignal())
//            <condition to write subset>,  //change this line, specify the condition that must be satisfied for the subset to be written to file
//                                         //(for example: anp.didSignal())
//            eventWeight);   //DO NOT change this line, always pass this variable


 //The following subsets group nu+nubar events into the same histograms
  //MC Signal Nue CC Events
   SubsetTask(action,setName,"Nue", &Histograms::HistosInSubset, anp, anp.isSignal(), anp.didSignal(), eventWeight);
  //MC Tau CC Events
   SubsetTask(action,setName,"Tau", &Histograms::HistosInSubset, anp, anp.isTau(), anp.didTau(), eventWeight);
  //MC NC Events
   SubsetTask(action,setName,"NC", &Histograms::HistosInSubset, anp, anp.isNC() && !(anp.isTau() || anp.isSignal()), anp.didNC() && !(anp.didTau() || anp.didSignal()), eventWeight);
  //MC numu+numubar CC Events
   SubsetTask(action,setName,"CC", &Histograms::HistosInSubset, anp, anp.isCC() && !(anp.isTau() || anp.isSignal()), anp.didCC() && !(anp.didTau() || anp.didSignal()), eventWeight);
  //MC Beam Nue CC Events
   SubsetTask(action,setName,"BeamNue", &Histograms::HistosInSubset, anp, anp.isBeamNue() && !(anp.isTau() || anp.isSignal()), anp.didBeamNue() && !(anp.didTau() || anp.didSignal()), eventWeight);
  //For data this is all events or for mc this is all NC+CC+BeamNue
   SubsetTask(action,setName,"All", &Histograms::HistosInSubset, anp, true, !(anp.didTau() || anp.didSignal()), eventWeight);

 //The following subsets only have nu events in the histograms
  //MC Signal Nue CC Events
   SubsetTask(action,setName,"Nue_Nu", &Histograms::HistosInSubset, anp, anp.isNu() && anp.isSignal(), anp.didSignal(), eventWeight);
  //MC Tau CC Events
   SubsetTask(action,setName,"Tau_Nu", &Histograms::HistosInSubset, anp, anp.isNu() && anp.isTau(), anp.didTau(), eventWeight);
  //MC NC Events
   SubsetTask(action,setName,"NC_Nu", &Histograms::HistosInSubset, anp, anp.isNu() && anp.isNC() && !(anp.isTau() || anp.isSignal()), anp.didNC() && !(anp.didTau() || anp.didSignal()), eventWeight);
  //MC numu CC Events
   SubsetTask(action,setName,"CC_Nu", &Histograms::HistosInSubset, anp, anp.isNu() && anp.isCC() && !(anp.isTau() || anp.isSignal()), anp.didCC() && !(anp.didTau() || anp.didSignal()), eventWeight);
  //MC Beam Nue CC Events
   SubsetTask(action,setName,"BeamNue_Nu", &Histograms::HistosInSubset, anp, anp.isNu() && anp.isBeamNue() && !(anp.isTau() || anp.isSignal()), anp.didBeamNue() && !(anp.didTau() || anp.didSignal()), eventWeight);

 //The following subsets only have nubar events in the histograms
  //MC Signal Nue CC Events
   SubsetTask(action,setName,"Nue_NuBar", &Histograms::HistosInSubset, anp, anp.isNuBar() && anp.isSignal(), anp.didSignal(), eventWeight);
  //MC Tau CC Events
   SubsetTask(action,setName,"Tau_NuBar", &Histograms::HistosInSubset, anp, anp.isNuBar() && anp.isTau(), anp.didTau(), eventWeight);
  //MC NC Events
   SubsetTask(action,setName,"NC_NuBar", &Histograms::HistosInSubset, anp, anp.isNuBar() && anp.isNC() && !(anp.isTau() || anp.isSignal()), anp.didNC() && !(anp.didTau() || anp.didSignal()), eventWeight);
  //MC numubar CC Events
   SubsetTask(action,setName,"CC_NuBar", &Histograms::HistosInSubset, anp, anp.isNuBar() && anp.isCC() && !(anp.isTau() || anp.isSignal()), anp.didCC() && !(anp.didTau() || anp.didSignal()), eventWeight);
  //MC Beam Nue CC Events
   SubsetTask(action,setName,"BeamNue_NuBar", &Histograms::HistosInSubset, anp, anp.isNuBar() && anp.isBeamNue() && !(anp.isTau() || anp.isSignal()), anp.didBeamNue() && !(anp.didTau() || anp.didSignal()), eventWeight);

}




//------------------------------------------------------------------------------
//
// Defintion for Histograms::HistosInSubset
//
// In this method, you define the histograms that are in a subset
// To define a histogram you must give:
//   1)a name for the histogram (The class will handle making it unique for various sets and subsets)
//   2)the binning for the histogram
//   3)the quantity to store in the histogram
//
//
//------------------------------------------------------------------------------
void Histograms::HistosInSubset(TString action, TString setName, TString subSetName, NueUtilities::AnaNueProcessor & anp, double eventWeight)
{
// New histograms are created by using the HistoTask method
// For each new histogram you want, you need to call HistoTask
// This is how you use HistoTask:
//
//     HistoTask(
//               action,setName,subSetName,   //DO NOT change this line, always pass these variables
//               <histogram name>,  //change this line, specify the histogram name (for example: "trueE") 
//               numXBins, lowXVal, highXval, //change this line, specify the binning for the histogram
//                                            //This is for a 1D histogram for a 2D you will need
//                                            //to specify numYBins, lowYVal, highYval, after the X values
//                                            //for a 3D histogram, you will have to specify numZBins,
//                                            //lowZVal, highZval, after the Y values
//               anp.GetTrueE(), //change this line, specify the quantity to fill (for example: anp.GetTrueE())
//               eventWeight);   //DO NOT change this line, always pass this variable
//


    //
    //Some 1D histos
    //
    //true E distribution
     HistoTask(action,setName,subSetName,"trueE", 80, 0, 20, anp.GetTrueE(), eventWeight);
    //reco E distribution
     HistoTask(action,setName,subSetName,"recoE", 9, 0, 9, anp.GetRecoE(), eventWeight);
    //ann11 distribution
     HistoTask(action,setName,subSetName,"ann11", 20, 0, 1, NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN2PE_DAIKON04), eventWeight);

    //ann11 variables
     HistoTask(action,setName,subSetName,"ann_par_a", 50,0,8, anp.nueRecord->shwfit.par_a, eventWeight);
     HistoTask(action,setName,subSetName,"ann_par_b", 50,0,2, anp.nueRecord->shwfit.par_b, eventWeight);
     HistoTask(action,setName,subSetName,"ann_uv_molrad", 50,0,20, anp.nueRecord->shwfit.uv_molrad_peak_9s_2pe_dw, eventWeight);
     HistoTask(action,setName,subSetName,"ann_uv_rms", 50,0,8, anp.nueRecord->shwfit.uv_rms_9s_2pe_dw, eventWeight);
     HistoTask(action,setName,subSetName,"ann_mstvars", 50,0,400, anp.nueRecord->mstvars.e4w+anp.nueRecord->mstvars.o4w, eventWeight);
     HistoTask(action,setName,subSetName,"ann_fract_2_planes", 50,0,1, anp.nueRecord->fracvars.fract_2_planes, eventWeight);
     HistoTask(action,setName,subSetName,"ann_fract_4_planes", 50,0,1, anp.nueRecord->fracvars.fract_4_planes, eventWeight);
     HistoTask(action,setName,subSetName,"ann_fract_6_planes", 50,0,1, anp.nueRecord->fracvars.fract_6_planes, eventWeight);
     HistoTask(action,setName,subSetName,"ann_fract_8_counters", 50,0,1, anp.nueRecord->fracvars.fract_8_counters, eventWeight);
     HistoTask(action,setName,subSetName,"ann_fract_road", 50,0,1, anp.nueRecord->fracvars.fract_road, eventWeight);
     HistoTask(action,setName,subSetName,"ann_LongE", 50,0,200, anp.nueRecord->shwfit.LongE, eventWeight);

    //
    //Some 2D histos
    //
    //(true E, reco E) distribution
     HistoTask(action,setName,subSetName,"trueE_recoE", 80, 0, 20, 9, 0, 9, anp.GetTrueE(), anp.GetRecoE(), eventWeight);
    //(true E, ann11) distribution
     HistoTask(action,setName,subSetName,"trueE_ann11", 80, 0, 20, 20, 0, 1, anp.GetTrueE(), NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN2PE_DAIKON04), eventWeight);
    //(reco E, ann11) distribution
     HistoTask(action,setName,subSetName,"recoE_ann11", 9, 0, 9, 20, 0, 1, anp.GetRecoE(), NueStandard::GetPIDValue(anp.nueRecord, Selection::kANN2PE_DAIKON04), eventWeight);
}








