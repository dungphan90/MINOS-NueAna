///////////////////////////////////////////////
////
////            XTalkFilter
////
////    This module will attempt to remove xtalk hits and put
////    them back in the strips where they belong
////
////
////    Code includes code by from XTalkRemover from Pedro Ochoa,
////    April 26 2007
////
////    Steven Cavanaugh, April 15, 2009
////
///////////////////////////////////////////////

#ifndef XTalkFilter_H
#define XTalkFilter_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include <vector>
#include "TTree.h"

#include "StandardNtuple/NtpStRecord.h"
#include "Plex/PlexHandle.h"
#include "Plex/PlexPixelSpotId.h"
#include "NueAna/NueAnalysisCuts.h"

#include "Pixel.h"
#include "PixelSpot.h"
#include "Conventions/Detector.h"



using namespace std;

class PMT;

class TH1F;
class TFile;
class NtpStRecord;

const Int_t MAX_NUMBER_OF_PMTS = 1800;

class XTalkFilter : public JobCModule
{

public:

  XTalkFilter();
  ~XTalkFilter();

public:

// Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);
  void EndJob();
  void BeginJob();
  void Config(const Registry&);
  const Registry& DefaultConfig() const;

        void Initialize();

        PMT * GetPMT(std::vector<PMT> &v, int id, int create=1);

        int GetHops(PMT *pmt, Pixel*startPixel);

        int GetHops(PMT *pmt, std::vector<int> &pixels, int pixelID,
            Pixel* startPixel,std::vector<int> visited);

        static int printit;
	void DoTrueFilter(NtpStRecord * record, int ievt);


private:

  PlexHandle *fPlexHandle;
  Int_t _pixelToVaChannel[16];
  Float_t _pixelSpotXTalkMap[8][9]; 

  Int_t CalculatePMTmuxID(PlexPixelSpotId &ppsid);
  int method;
  
  std::vector<PMT> east;
  std::vector<PMT> west;
  
  std::vector<std::vector<std::pair<int, int> > >pixelmaps;
  
  Detector::Detector_t lastpixelmap;
  void GetPixelMaps();


  void SaveAdjustment(std::vector<PMT> & pmts, int side,
                      NtpStRecord * record);
 
  void CleanXTalk(std::vector<PMT> & pmts, NtpStRecord * record);
  
// std::vector<std::pair<int,int > > GetXTalkContributors(int pixel);

   std::vector<std::pair<int,int > > GetXTalkPixelsFar(int pixel);

   std::vector<std::pair<int,int > > GetXTalkPixelsNear(int pixel);
 

   Detector::Detector_t detector;
 
};
#endif // COMPAREALL_H
////////////////////////////////////////////////////////////////////////
