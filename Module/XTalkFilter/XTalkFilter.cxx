#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include "TH1F.h"
#include "TFile.h"

#include "NueAna/NueRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro

#include "AnalysisNtuples/ANtpDefaultValue.h"                                   
#include "TClass.h"
#include "TDataType.h"
#include "TDataMember.h"
#include "TRealData.h"
#include "XTalkFilter.h"
#include "StandardNtuple/NtpStRecord.h"

#include "CandNtupleSR/NtpSRRecord.h"

#include "TClonesArray.h"

#include "CandNtupleSR/NtpSRStrip.h"  
#include "CandNtupleSR/NtpSREvent.h" 
#include "CandNtupleSR/NtpSRShower.h" 
#include "CandNtupleSR/NtpSRTrack.h" 

#include "TObject.h"

#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"

#include "PMT.h"

JOBMODULE(XTalkFilter, "XTalkFilter",
          "Filter out strips below a given PE");
          
CVSID("$Id: XTalkFilter.cxx,v 1.9 2009/06/23 22:50:10 scavan Exp $");

int XTalkFilter::printit=-1;

// helper functions 
// to get an array of pairs of xtalk pixel and index into
// pixelspot xtalk map for a given true pixel


        
std::vector<std::pair<int,int > > XTalkFilter::GetXTalkPixelsNear(int pixel)
{     
        std::vector<std::pair<int,int> > ret;
                


  if(pixel==0)
  {
        ret.push_back(std::pair<int,int>(1,5));
        ret.push_back(std::pair<int,int>(8,7));
        ret.push_back(std::pair<int,int>(9,8));
  }else if(pixel==7){
        ret.push_back(std::pair<int,int>(6,3));
        ret.push_back(std::pair<int,int>(14,6));
        ret.push_back(std::pair<int,int>(15,7));
  }else if(pixel==56){
        ret.push_back(std::pair<int,int>(48,1));
        ret.push_back(std::pair<int,int>(49,2));
        ret.push_back(std::pair<int,int>(57,5));
  }else if(pixel==63){
        ret.push_back(std::pair<int,int>(62,3));
        ret.push_back(std::pair<int,int>(54,0));
        ret.push_back(std::pair<int,int>(55,1));
  }

  //first row
  else if(pixel>0 && pixel<7)
  {  
        ret.push_back(std::pair<int,int>(pixel-1,3));
        ret.push_back(std::pair<int,int>(pixel+1,5));
        ret.push_back(std::pair<int,int>(pixel+7,6));
        ret.push_back(std::pair<int,int>(pixel+8,7));
        ret.push_back(std::pair<int,int>(pixel+9,8));
  }

  //last orw
  else if(pixel>56 && pixel<63)
  { 
        ret.push_back(std::pair<int,int>(pixel-1,3));
        ret.push_back(std::pair<int,int>(pixel+1,5));
        ret.push_back(std::pair<int,int>(pixel-7,2));
        ret.push_back(std::pair<int,int>(pixel-8,1));
        ret.push_back(std::pair<int,int>(pixel-9,0));
  }   
    
  //left column

  else if(pixel%8==0 && pixel>0 && pixel<56)
  { 
        ret.push_back(std::pair<int,int>(pixel-8,1));
        ret.push_back(std::pair<int,int>(pixel-7,2));
        ret.push_back(std::pair<int,int>(pixel+1,5));
        ret.push_back(std::pair<int,int>(pixel+8,7));
        ret.push_back(std::pair<int,int>(pixel+9,8));
  }   
    
  //right column
  else if((pixel+1)%8==0 && pixel>7 && pixel<63)
  {
        ret.push_back(std::pair<int,int>(pixel-8,1));
        ret.push_back(std::pair<int,int>(pixel-9,0));
        ret.push_back(std::pair<int,int>(pixel-1,3));
        ret.push_back(std::pair<int,int>(pixel+7,6));
        ret.push_back(std::pair<int,int>(pixel+8,7));
  }

  //interior pixel
  else{
        ret.push_back(std::pair<int,int>(pixel-9,0));
        ret.push_back(std::pair<int,int>(pixel-8,1));
        ret.push_back(std::pair<int,int>(pixel-7,2));
        ret.push_back(std::pair<int,int>(pixel-1,3));
        ret.push_back(std::pair<int,int>(pixel+1,5));
        ret.push_back(std::pair<int,int>(pixel+7,6));
        ret.push_back(std::pair<int,int>(pixel+8,7));
        ret.push_back(std::pair<int,int>(pixel+9,8));
  }  
        return ret;
}       
    




std::vector<std::pair<int,int > > XTalkFilter::GetXTalkPixelsFar(int pixel)
{
        std::vector<std::pair<int,int> > ret;

  switch(pixel){
    case 0:
        ret.push_back(std::pair<int,int>(1,5));
        ret.push_back(std::pair<int,int>(4,7));
        ret.push_back(std::pair<int,int>(5,8));
      break;
    case 1:
        ret.push_back(std::pair<int,int>(0,3));
        ret.push_back(std::pair<int,int>(2,5));
        ret.push_back(std::pair<int,int>(4,6));
        ret.push_back(std::pair<int,int>(5,7));
        ret.push_back(std::pair<int,int>(6,8));
      break;
    case 2:
                ret.push_back(std::pair<int,int>(1,3));
        ret.push_back(std::pair<int,int>(3,5));
        ret.push_back(std::pair<int,int>(5,6));
        ret.push_back(std::pair<int,int>(6,7));
        ret.push_back(std::pair<int,int>(7,8));
      break;
    case 3:
                ret.push_back(std::pair<int,int>(2,3));
        ret.push_back(std::pair<int,int>(6,6));
        ret.push_back(std::pair<int,int>(7,7));
      break;
    case 4:
                ret.push_back(std::pair<int,int>(0,1));
        ret.push_back(std::pair<int,int>(1,2));
        ret.push_back(std::pair<int,int>(5,5));
        ret.push_back(std::pair<int,int>(8,7));
        ret.push_back(std::pair<int,int>(9,8));
      break;
    case 5:
                ret.push_back(std::pair<int,int>(0,0));
        ret.push_back(std::pair<int,int>(1,1));
        ret.push_back(std::pair<int,int>(2,2));
        ret.push_back(std::pair<int,int>(4,3));
        ret.push_back(std::pair<int,int>(6,5));
        ret.push_back(std::pair<int,int>(8,6));
        ret.push_back(std::pair<int,int>(9,7));
        ret.push_back(std::pair<int,int>(10,8));
      break;
    case 6:
                ret.push_back(std::pair<int,int>(1,0));
        ret.push_back(std::pair<int,int>(2,1));
        ret.push_back(std::pair<int,int>(3,2));
        ret.push_back(std::pair<int,int>(5,3));
        ret.push_back(std::pair<int,int>(7,5));
        ret.push_back(std::pair<int,int>(9,6));
        ret.push_back(std::pair<int,int>(10,7));
        ret.push_back(std::pair<int,int>(11,8));
        break;
    case 7:
                ret.push_back(std::pair<int,int>(2,0));
        ret.push_back(std::pair<int,int>(3,1));
        ret.push_back(std::pair<int,int>(6,3));
        ret.push_back(std::pair<int,int>(10,6));
        ret.push_back(std::pair<int,int>(11,7));
      break;
    case 8:
                ret.push_back(std::pair<int,int>(4,1));
        ret.push_back(std::pair<int,int>(5,2));
        ret.push_back(std::pair<int,int>(9,5));
        ret.push_back(std::pair<int,int>(12,7));
        ret.push_back(std::pair<int,int>(13,8));
      break;
    case 9:
                ret.push_back(std::pair<int,int>(4,0));
        ret.push_back(std::pair<int,int>(5,1));
        ret.push_back(std::pair<int,int>(6,2));
        ret.push_back(std::pair<int,int>(8,3));
        ret.push_back(std::pair<int,int>(10,5));
        ret.push_back(std::pair<int,int>(12,6));
        ret.push_back(std::pair<int,int>(13,7));
        ret.push_back(std::pair<int,int>(14,8));
      break;
    case 10:
                ret.push_back(std::pair<int,int>(5,0));
        ret.push_back(std::pair<int,int>(6,1));
        ret.push_back(std::pair<int,int>(7,2));
        ret.push_back(std::pair<int,int>(9,3));
        ret.push_back(std::pair<int,int>(11,5));
        ret.push_back(std::pair<int,int>(13,6));
        ret.push_back(std::pair<int,int>(14,7));
        ret.push_back(std::pair<int,int>(15,8));
      break;
    case 11:
                ret.push_back(std::pair<int,int>(6,0));
        ret.push_back(std::pair<int,int>(7,1));
        ret.push_back(std::pair<int,int>(10,3));
        ret.push_back(std::pair<int,int>(14,6));
        ret.push_back(std::pair<int,int>(15,7));
      break;
    case 12:
                ret.push_back(std::pair<int,int>(8,1));
        ret.push_back(std::pair<int,int>(9,2));
        ret.push_back(std::pair<int,int>(13,5));
      break;
    case 13:
                ret.push_back(std::pair<int,int>(8,0));
        ret.push_back(std::pair<int,int>(9,1));
        ret.push_back(std::pair<int,int>(10,2));
        ret.push_back(std::pair<int,int>(12,3));
        ret.push_back(std::pair<int,int>(14,5));
      break;
    case 14:
                ret.push_back(std::pair<int,int>(9,0));
        ret.push_back(std::pair<int,int>(10,1));
        ret.push_back(std::pair<int,int>(11,2));
        ret.push_back(std::pair<int,int>(13,3));
        ret.push_back(std::pair<int,int>(15,5));
      break;
    case 15:
                ret.push_back(std::pair<int,int>(10,0));
        ret.push_back(std::pair<int,int>(11,1));
        ret.push_back(std::pair<int,int>(14,3));
      break;
    default:
      break;
        }

        return ret;
}

/////////////////////////////////////////////////////////

void XTalkFilter::Initialize()
{

// Necessary arrays (obtained from AltDeMux)
  _pixelToVaChannel[0]  = 3;
  _pixelToVaChannel[1]  = 5;
  _pixelToVaChannel[2]  = 14;
  _pixelToVaChannel[3]  = 16;
  _pixelToVaChannel[4]  = 7;
  _pixelToVaChannel[5]  = 9;
  _pixelToVaChannel[6]  = 10;
  _pixelToVaChannel[7]  = 12;
  _pixelToVaChannel[8]  = 11;
  _pixelToVaChannel[9]  = 13;
  _pixelToVaChannel[10] = 6;
  _pixelToVaChannel[11] = 8;
  _pixelToVaChannel[12] = 17;
  _pixelToVaChannel[13] = 15;
  _pixelToVaChannel[14] = 2;
  _pixelToVaChannel[15] = 4;


//not actually using _pixelSpotXTalkMap now...
/*

  //-->spot 0
  _pixelSpotXTalkMap[2][0] = 0.007; 
  _pixelSpotXTalkMap[2][1] = 0.025;
  _pixelSpotXTalkMap[2][2] = 0.002;
  _pixelSpotXTalkMap[2][3] = 0.025;
  _pixelSpotXTalkMap[2][5] = 0.003;
  _pixelSpotXTalkMap[2][6] = 0.001;
  _pixelSpotXTalkMap[2][7] = 0.003;
  _pixelSpotXTalkMap[2][8] = 0.002;
  
  //-->spot 1
  _pixelSpotXTalkMap[1][0] = 0.003; 
  _pixelSpotXTalkMap[1][1] = 0.025;
  _pixelSpotXTalkMap[1][2] = 0.003;
  _pixelSpotXTalkMap[1][3] = 0.006;
  _pixelSpotXTalkMap[1][5] = 0.006;
  _pixelSpotXTalkMap[1][6] = 0.001;
  _pixelSpotXTalkMap[1][7] = 0.004;
  _pixelSpotXTalkMap[1][8] = 0.001;
  

  //-->spot 2
  _pixelSpotXTalkMap[0][0] = 0.002; 
  _pixelSpotXTalkMap[0][1] = 0.025;
  _pixelSpotXTalkMap[0][2] = 0.007;
  _pixelSpotXTalkMap[0][3] = 0.003;
  _pixelSpotXTalkMap[0][5] = 0.025;
  _pixelSpotXTalkMap[0][6] = 0.001;
  _pixelSpotXTalkMap[0][7] = 0.003;
  _pixelSpotXTalkMap[0][8] = 0.002;
  
  //-->spot 3
  _pixelSpotXTalkMap[4][0] = 0.002; 
  _pixelSpotXTalkMap[4][1] = 0.008;
  _pixelSpotXTalkMap[4][2] = 0.001;
  _pixelSpotXTalkMap[4][3] = 0.011;
  _pixelSpotXTalkMap[4][5] = 0.004;
  _pixelSpotXTalkMap[4][6] = 0.002;
  _pixelSpotXTalkMap[4][7] = 0.008;
  _pixelSpotXTalkMap[4][8] = 0.001;
  

  //-->spot 4
  _pixelSpotXTalkMap[3][0] = 0.001; 
  _pixelSpotXTalkMap[3][1] = 0.008;
  _pixelSpotXTalkMap[3][2] = 0.002;
  _pixelSpotXTalkMap[3][3] = 0.004;
  _pixelSpotXTalkMap[3][5] = 0.011;
  _pixelSpotXTalkMap[3][6] = 0.001;
  _pixelSpotXTalkMap[3][7] = 0.008;
  _pixelSpotXTalkMap[3][8] = 0.002;
  
  //-->spot 5
  _pixelSpotXTalkMap[7][0] = 0.002; 
  _pixelSpotXTalkMap[7][1] = 0.003;
  _pixelSpotXTalkMap[7][2] = 0.001;
  _pixelSpotXTalkMap[7][3] = 0.025;
  _pixelSpotXTalkMap[7][5] = 0.003;
  _pixelSpotXTalkMap[7][6] = 0.007;
  _pixelSpotXTalkMap[7][7] = 0.025;
  _pixelSpotXTalkMap[7][8] = 0.002;
  
  //-->spot 6
  _pixelSpotXTalkMap[6][0] = 0.001; 
  _pixelSpotXTalkMap[6][1] = 0.004;
  _pixelSpotXTalkMap[6][2] = 0.001;
  _pixelSpotXTalkMap[6][3] = 0.006;
  _pixelSpotXTalkMap[6][5] = 0.006;
  _pixelSpotXTalkMap[6][6] = 0.003;
  _pixelSpotXTalkMap[6][7] = 0.025;
  _pixelSpotXTalkMap[6][8] = 0.003;
  
  //-->spot 7
  _pixelSpotXTalkMap[5][0] = 0.001; 
  _pixelSpotXTalkMap[5][1] = 0.003;
  _pixelSpotXTalkMap[5][2] = 0.002;
  _pixelSpotXTalkMap[5][3] = 0.003;
  _pixelSpotXTalkMap[5][5] = 0.025;
  _pixelSpotXTalkMap[5][6] = 0.002;
  _pixelSpotXTalkMap[5][7] = 0.025;
  _pixelSpotXTalkMap[5][8] = 0.007;

*/  
 


}

PMT * XTalkFilter::GetPMT(std::vector<PMT> &v, int id, int create)
{
        for(unsigned int i=0;i<v.size();i++)
        {
                if(v[i].GetId()==id)return &(v[i]);
        }
        
        if(!create)return 0;
        
        //if its not there... create it
        PMT newpmt(id);
        v.push_back(newpmt);
        return & (v[v.size()-1]);

}

void XTalkFilter::DoTrueFilter(NtpStRecord * record, int ievt){

// clear the pmt vectors
        east.clear();
        west.clear();

	detector = record->GetHeader().GetVldContext().GetDetector();

	//make sure we are using the maps for the correct detector
	GetPixelMaps();

   NtpSREvent* event = dynamic_cast<NtpSREvent*>(record->evt->At(ievt));  
   int nstrips=event->nstrip;

// Getting a handle on the plex
   fPlexHandle = new PlexHandle(record->GetHeader().GetVldContext());

   if (record!=NULL){
    
// -->strips loop
     for(int istp = 0; istp< nstrips; istp++){

// Getting east and west strip-ends and PMT pixels
     const NtpSRStrip* strip = dynamic_cast<const NtpSRStrip*>
                          (record->stp->At(event->stp[istp]));

     if (!strip)
                {
       printf("missing strip!!! idx %d total %d\n",event->stp[istp],
                                         record->stp->GetEntries());
                
                }

      PlexStripEndId pseid_east(detector,
                     strip->plane,strip->strip,StripEnd::kEast);
      PlexStripEndId pseid_west(detector,
		     strip->plane,strip->strip,StripEnd::kWest);

      PlexPixelSpotId ppsid_east = fPlexHandle->GetPixelSpotId(pseid_east); 
      PlexPixelSpotId ppsid_west = fPlexHandle->GetPixelSpotId(pseid_west);

      Int_t PMTmuxID_east = CalculatePMTmuxID(ppsid_east);
      Int_t PMTmuxID_west = CalculatePMTmuxID(ppsid_west);

      Float_t strip_e = (strip->ph0.sigcor+strip->ph1.sigcor)/60.;

// East
//----
      Float_t q  = strip->ph0.pe;
      Float_t qc = strip->ph0.sigcor/60.;
      if(q!=0){
        if(PMTmuxID_east<MAX_NUMBER_OF_PMTS && !pseid_east.IsVetoShield()){
 
// if(qc/q<0.1 || qc/q > 10.0){
//   qc=q;
// }
// UpdatePixelMap(PMTmuxID_east,ppsid_east.GetPixel(),qc);
// UpdateXTalkMap(PMTmuxID_east,ppsid_east.GetPixel(), ppsid_east.GetSpot(),qc);         
                PMT * pmte = GetPMT(east,PMTmuxID_east);
                pmte->GetPixel(ppsid_east.GetPixel())->GetPixelSpot(ppsid_east.GetSpot())->AddStrip(event->stp[istp], qc, strip_e);
 
        } else {
          cout << "Warning in XTalkRemover: Had a PMTmuxID higher than MAX_NUMBER_OF_PMTS and/or a veto shield hit" << endl;
        }//<--not bigger than MAX_NUMBER_OF_PMTS, and not veto shield
      }//<--no charge on east side

      //West
      //----
      q  = strip->ph1.pe;
      qc = strip->ph1.sigcor/60.;
      if(q!=0){
        if(PMTmuxID_west<MAX_NUMBER_OF_PMTS && !pseid_west.IsVetoShield()){       
//        if(qc/q<0.1 || qc/q > 10.0){
//          qc=q;
//        }
//        UpdatePixelMap(PMTmuxID_west,ppsid_west.GetPixel(),qc);
//        UpdateXTalkMap(PMTmuxID_west,ppsid_west.GetPixel(), ppsid_west.GetSpot(),qc);
                PMT * pmtw = GetPMT(west,PMTmuxID_west);
                pmtw->GetPixel(ppsid_west.GetPixel())->GetPixelSpot(ppsid_west.GetSpot())->AddStrip(event->stp[istp], qc, strip_e);

        } else {
          cout << "Warning in XTalkRemover: Had a PMTmuxID higher than MAX_NUMBER_OF_PMTS and/or a veto shield hit" << endl;
        }//<--not bigger than MAX_NUMBER_OF_PMTS, and not veto shield
      }//<--no charge on west side
      
    }//strips loop
    
  } else {//<--not empty record
    cout << "Warning: Passed empty record to MakePixelMap in XTalkRemover" << endl;
  }


/*
        printf("--------------------------\nEast PMTs\n");
        for(unsigned int i=0;i<east.size();i++)
        {
                east[i].Dump();
        }
        printf("--------------------------\n");



        printf("--------------------------\nWest PMTs\n");
        for(unsigned int i=0;i<west.size();i++)
        {
                west[i].Dump();
        }
        printf("--------------------------\n");
*/      
        
       if(printit)printf("CleanXTalk east\n"); 
        CleanXTalk(east,record);
       if(printit)printf("CleanXTalk west\n"); 
        CleanXTalk(west,record);        
        
        SaveAdjustment(east,0,record);
        SaveAdjustment(west,1,record);
        
}

void XTalkFilter::SaveAdjustment(std::vector<PMT> & pmts, int side, NtpStRecord * record)
{

// iterate over pmts that have changed spots
        for(unsigned int i=0;i<pmts.size();i++)
        {
                PMT* pmt = &pmts[i];
                if(!pmt->GetStatus())continue; //not changed

                std::vector<int>pixels=pmt->GetPixelVector();
                for(unsigned int j=0;j<pixels.size();j++)
                {
                        Pixel* p = pmt->GetPixel(pixels[j],0);
                        if(!p)continue;
                        if(!p->GetStatus())continue;
                                
                        std::vector<int>  psv = p->GetPixelSpotVector();
                        for(unsigned int l=0;l<psv.size();l++)
                        {
                                PixelSpot * ps=p->GetPixelSpot(psv[l],0);
                                if(!ps)continue;
                                if(!ps->GetStatus())continue;
                                        

				double adjfrac=0;
				double totalstrippe=0;
                                std::vector<int> strips = ps->GetStrips();
                                for(unsigned int m=0;m<strips.size();m++)
                                {                             
					NtpSRStrip* strip = (NtpSRStrip*)(record->stp->At(strips[m]));  
                                	totalstrippe += side==0?strip->ph0.sigcor/60.:strip->ph1.sigcor/60.;
				}

                                double thispe = ps->GetE();
				adjfrac = totalstrippe>0?thispe / totalstrippe : 0;
				if(adjfrac<1e-5)adjfrac=0; //don't want small energy... don't want <0 energy
			


if(printit)printf("adj pmt %d p %d ps %d ee %f adjfrac %f stripe %f\n",pmt->GetId(), p->GetId(), ps->GetId(), thispe, adjfrac, totalstrippe);
				
                                for(unsigned int m=0;m<strips.size();m++)
                                {                               


        
                                NtpSRStrip* strip = (NtpSRStrip*)(record->stp->At(strips[m]));
                                
                                if(!strip)
                                {
                                        printf("CANT FIND STRIP AT IDX %d\n",strips[m]);
                                        continue;
                                        
                                }
                                
                                    
                                   
                      
	
        
					//if(printit)printf("adjfrac %f  from old pe %f new pe %f\n",adjfrac,strippe,thispe);
                                        if(side==0)
                                        {
                                                if(printit)printf("adjusting frac %f strip %d pe %f to %f sc %f to %f\n",adjfrac,strips[m],strip->ph0.pe,strip->ph0.pe*adjfrac,strip->ph0.sigcor,strip->ph0.sigcor*adjfrac);
                                        
                                                strip->ph0.pe*=adjfrac;
						strip->ph0.siglin*=adjfrac;
                                                strip->ph0.sigcor*=adjfrac;
						strip->ph0.raw*=adjfrac;
                                        }else{
                                                if(printit)printf("adjusting frac %f strip %d pe %f to %f sc %f to %f\n",adjfrac,strips[m],strip->ph1.pe,strip->ph1.pe*adjfrac,strip->ph1.sigcor,strip->ph1.sigcor*adjfrac);
                                                
                                                strip->ph1.pe*=adjfrac;
						strip->ph1.siglin*=adjfrac;
                                                strip->ph1.sigcor*=adjfrac;
						strip->ph1.raw*=adjfrac;
                                        }
                                        
                                        
                                }
                        }
                }
        }       


}

int XTalkFilter::GetHops(PMT *pmt, Pixel* startPixel)
{
        std::vector<int>pixels=pmt->GetPixelVector();
        std::vector<int>visited;
        int nhops= GetHops(pmt, pixels,-1,startPixel,visited);
        if(nhops>=1000)nhops=-1; //not possible to get to mixed pixel!
        return nhops;
}


//visited vector must not be passed by reference...
int XTalkFilter::GetHops(PMT *pmt, std::vector<int> &pixels, int pixelID, Pixel* startPixel,std::vector<int> visited)
{

        if(pixelID==startPixel->GetId())return 1000; //circular loop... check for this value before using hops!.
        if(pixelID<0)pixelID=startPixel->GetId();
        visited.push_back(pixelID);

        Pixel* p = pmt->GetPixel(pixelID,0);

        if(p->type==1)return 0;
        
        int minhops=1000;
        std::vector<std::pair<int,int> > * pmap=&pixelmaps[pixelID];
        for(unsigned int k=0;k<pmap->size();k++)
        {
                int pid = (*pmap)[k].first;
                Pixel * q = 0;
                
                //don't infinitely loop!
                int alreadydid=0;
                for(unsigned int j=0;j<visited.size();j++)
                if(visited[j]==pid){alreadydid=1;break;}
                if(alreadydid)continue;
                
                q=pmt->GetPixel(pid,0);
                if(!q)continue;
                
                int hops = GetHops(pmt,pixels,q->GetId(),startPixel,visited);
//              printf("hops %d start %d to %d from %d\n",hops,startPixel->GetId(), q->GetId(),pixelID);
                if(hops<minhops)minhops=hops;   
        }
        return minhops+1;
}

void XTalkFilter::CleanXTalk(std::vector<PMT> & pmts, NtpStRecord * /*record*/)
{
        //int print=1;
        
        double xtalke=0;
        double totale=0;
        
        int xtalkstrips=0;
        int xtalkpmts=0;


        //double xtalkthreshold=1.0;  //xtalk if meus in this side less than this number
        //double xtalksinglethreshold=2.0;  //xtalk if all strip meus on this side and meus less than this number
        //double contribFrac=0.1; //less than x % of possible contributing energy indicated crosstalk?


        for(unsigned int i=0;i<pmts.size();i++)
        {
                PMT* pmt = &pmts[i];
                if(pmt->GetNStrips()<2)continue;//no need to look for XTalk here...
        
                xtalkpmts++;
                xtalkstrips+=pmt->GetNStrips();
                
           //     if(printit)printf("\n/////////////////////////////////\nPMT id %d\n",pmt->GetId());



                
                std::vector<int>pixels=pmt->GetPixelVector();

                //iterate over xtalk pixels, and add back in this pixels energy to those pixels by fraction of what we expect....
                //perhaps than just a below a certain energy xtalk check... also require energy of pixel to be larger than some frac of energy of pixels that could have made this xtalk
                        
                        

                        
                        //make a multimap of xtalks to consider (energy, pixel idx)
                 //       multimap<float, int> xtalkcandidate;


                //attempt to classify the pixel as crosstalk only or mixed
                
		//Far det
		//crosstalk is defined as all strip energy only on this side of the detector
                //mixed is defined as energy on both sides of the detector (ie may contain xtalk)

		if(detector==Detector::kFar)
		{
                	for(unsigned int j=0;j<pixels.size();j++)
                	{
                        	Pixel* p = pmt->GetPixel(pixels[j],0);
                        	if(!p)continue;
                        
                        
                        	if(fabs(p->GetE() - p->GetStripE() )<1e-5)
                        	{
                                	//is xtalk
                                	p->type=0;
                        	}else{
                       	        	//is mixed
                       	         	p->type=1;
                 	       	}
                	}       
		}


		//Near det
		//since it is single sided... define an xtalk pixel as one that is surrounded by other pixels with an energy less than some fraction of the summed energy of the surrounding pixels
	
		double nearthresh=0.1;

		if(detector==Detector::kNear)
		{
                        for(unsigned int j=0;j<pixels.size();j++)
                        {
                                Pixel* p = pmt->GetPixel(pixels[j],0);
                                if(!p)continue;
				double neighborE=0;



	                        std::vector<std::pair<int,int> > * pmap=&pixelmaps[p->GetId()];
                       	 	for(unsigned int k=0;k<pmap->size();k++)
                        	{
                                	int pid = (*pmap)[k].first;

					Pixel *pn = pmt->GetPixel(pid,0);
					if(!pn)continue;					
					neighborE+=pn->GetE();
				}

				if(neighborE)
				{
					if((p->GetE()/neighborE)<nearthresh)p->type=0;
					else p->type=1;
				}else{
					p->type=1;
				}
			}
		}


/*
                for(unsigned int j=0;j<pixels.size();j++)
                {
                        Pixel* p = pmt->GetPixel(pixels[j],0);
                        if(!p)continue;
                        
                        double e=p->GetE();
                        if(e==0)continue;
                        
                        totale+=e;
                        
                //      if(e<xtalkthreshold || (e<xtalksinglethreshold && e == p->GetStripE()))
                //              xtalkcandidate.insert(std::pair<float,int>(e,p->GetId()));

                        //see what the energy of the pixels that could have contribued to xtalk here are
                        std::vector<std::pair<int,int> > * pmap=&pixelmaps[p->GetId()];

                
                        double contribE=0;
                        for(unsigned int k=0;k<pmap->size();k++)
                        {
                                int pid = (*pmap)[k].first;
                                Pixel * q = pmt->GetPixel(pid,0);
                                if(!q)continue;
                                contribE+=q->GetE();
                        }
                        
                
                        //single sided energy strip
                        if( (e<xtalksinglethreshold && e == p->GetStripE()))
                                xtalkcandidate.insert(std::pair<float,int>(e,p->GetId()));
                
                        //above some fraction of contributor energy
                        else if(e<contribE * contribFrac)
                                xtalkcandidate.insert(std::pair<float,int>(e,p->GetId()));
                
                
                }

*/
                std::multimap<int,int>hops;
                                                
                for(unsigned int j=0;j<pixels.size();j++)
                {
                        Pixel* p = pmt->GetPixel(pixels[j],0);
                        if(!p)continue;
                       totale+=p->GetE(); 
                        int nhops = GetHops(pmt,p); //how far to a mixed pixel?
                        p->hopstomixed=nhops;
                        hops.insert(std::pair<int,int>(nhops,pixels[j]));       
          //              if(printit)printf("minhops pixeld id %d hops %d\n",pixels[j],nhops);
                }       
                        
        
                        
                
///for printing out what is currently there....

//              for(unsigned int j=0;j<pixels.size();j++)
//              {
//                      Pixel* p = pmt->GetPixel(pixels[j],0);

//order of hops....





                std::multimap<int,int>::reverse_iterator ith;
                for(ith=hops.rbegin();ith!=hops.rend();ith++)
                {
                        Pixel *p = pmt->GetPixel(ith->second,0);
                        if(!p)continue;
                        
                        
                        double e=p->GetE();
                        if(e==0)continue;

                        
                        
                       // if(printit)printf("\t pixel id %d addr %p with e %f strips(%d) e %f type %d hops %d\n",p->GetId(),(void*)p,p->GetE(),p->GetNStrips(),p->GetStripE(),p->type,p->hopstomixed);
                        
                        //see what the energy of the pixels that could have contribued to xtalk here are
                        std::vector<std::pair<int,int> > * pmap=&pixelmaps[p->GetId()];
                        
                        std::vector<std::pair<int, int> > toadj; //pixel, pixelspot of possible contributers
                        std::vector<double> contribe; //current e of possible contributors
                        
                        double sum_e=0;
                        for(unsigned int k=0;k<pmap->size();k++)
                        {
                                int pid = (*pmap)[k].first;

                         //       if(printit)printf("\t\t xtalks to pixel %d\n",pid);

                                Pixel * q = pmt->GetPixel(pid,0);
                                if(!q)continue;
                                
                                //int dir = (*pmap)[k].second;
                                
                                std::vector<int>  psv = q->GetPixelSpotVector();

                                
                                for(unsigned int l=0;l<psv.size();l++)
                                {
                                        PixelSpot * ps=q->GetPixelSpot(psv[l],0);
                                        if(!ps)continue;
                                
                                
                                        double xe = ps->GetE() ;//* _pixelSpotXTalkMap[ps->GetId()-1][8-dir];//reverse the dir
                                
                                        toadj.push_back(std::pair<int,int>(pid,ps->GetId()));
                                        contribe.push_back(xe);
                                        sum_e+=xe;

                           //             if(printit)printf("\t\t\t pixelspot %d addr %p dir %d with xe %f    e %f  psxtm %f stripE %f\n",ps->GetId(),(void*)ps,8-dir,xe,ps->GetE(), _pixelSpotXTalkMap[ps->GetId()-1][8-dir],ps->GetStripE());
                                }
                        }
                        
                                        
                }



///////////////////

                
                


                //for(unsigned int j=0;j<pixels->size();j++)
                //{
                //      Pixel* p = &(*pixels)[j];

//              multimap<float, int>::iterator it;
//              for(it=xtalkcandidate.begin();it!=xtalkcandidate.end();it++)
//              {
//                      Pixel * p = pmt->GetPixel(it->second,0);


                for(ith=hops.rbegin();ith!=hops.rend();ith++)
                {
                        Pixel *p = pmt->GetPixel(ith->second,0);
                        if(!p)continue;
                        
                        double pixelE=p->GetE();
                        if(pixelE<=0)continue;
                        
                        if(p->hopstomixed<0)continue;//it doesn't connect to other pixels...
                        if(p->type==1)continue;//its not xtalk only... so we wouldn't know what to share...!
                        
                    //    if(printit)printf("considering pixel id %d addr %p with e %f\n",p->GetId(),(void*)p,p->GetE());
                        
                        //see what the energy of the pixels that could have contribued to xtalk here are
                        std::vector<std::pair<int,int> > * pmap=&pixelmaps[p->GetId()];
                        
                        std::vector<std::pair<int, int> > toadj; //pixel, pixelspot of possible contributers
                        std::vector<double> contribe; //current e of possible contributors
                        
                        double sum_e=0;
                        for(unsigned int k=0;k<pmap->size();k++)
                        {
                                int pid = (*pmap)[k].first;


                                Pixel * q = pmt->GetPixel(pid,0);
                                if(!q)continue;
                                
                                double qe=q->GetE();
                                if(qe<1e-5)continue; //don't pass xtalk back to now empty xtalk hit
                                
                                //only pass up the chain
                                if(q->hopstomixed>p->hopstomixed)continue;
                                
                                
                //              if((qe<xtalksinglethreshold && qe == q->GetStripE() ) )continue;

                      //          if(printit)printf("\t\t xtalks to pixel %d\n",pid);

                                
                                //int dir = (*pmap)[k].second;
                                
                                std::vector<int>  psv = q->GetPixelSpotVector();

                                
                                for(unsigned int l=0;l<psv.size();l++)
                                {
                                        PixelSpot*ps = q->GetPixelSpot(psv[l],0);
                                        if(!ps)continue;
                                        
                                        if(ps->GetStatus())continue;
                                
                                        double xe = ps->GetE();//*_pixelSpotXTalkMap[ps->GetId()-1][8-dir];//reverse the dir
                                
                                        toadj.push_back(std::pair<int,int>(pid,ps->GetId()));
                                        contribe.push_back(xe);
                                        sum_e+=xe;

                        //                if(printit)printf("\t\t\t pixelspot %d addr %p dir %d with xe %f    e %f  psxtm %f stripE %f\n",ps->GetId(),(void*)ps,8-dir,xe,ps->GetE(), _pixelSpotXTalkMap[ps->GetId()-1][8-dir],ps->GetStripE());
                                }
                        }
                        
                        if(sum_e>0)
                        for(unsigned int k=0;k<toadj.size();k++)
                        {
                                
                                double toadd=contribe[k]/sum_e*pixelE;
                                p->TakeE(toadd);
                                Pixel * q = pmt->GetPixel(toadj[k].first,0);
                                Pixel oldp = *q;
                                
                                if(!q)continue;
                                PixelSpot *ps = q->GetPixelSpot(toadj[k].second,0);
                                if(!ps)continue;
                                ps->AddE(toadd);
                          //      if(printit)printf("q addr %p\n",(void*)q);
                                
                            //    if(printit)printf("adjusting pixel id %d by %f  e %f %f stripe %f %f type %d\n",q->GetId(),toadd,oldp.GetE(),q->GetE(),oldp.GetStripE(),q->GetStripE(),q->type);
                                xtalke+=toadd;
                        }
                        
                                        
                }
                
                

///for printing out what is currently there....
/*
                for(ith=hops.rbegin();ith!=hops.rend();ith++)
                {
                        Pixel *p = pmt->GetPixel(ith->second,0);
                        if(!p)continue;

                        
                        double e=p->GetE();
                        if(e==0)continue;
                        
                        
                        if(printit)printf("\t pixel id %d with e %f strips(%d) e %f type %d hops %d\n",p->GetId(),p->GetE(),p->GetNStrips(),p->GetStripE(),p->type,p->hopstomixed);
                        
                        //see what the energy of the pixels that could have contribued to xtalk here are
                        std::vector<std::pair<int,int> > * pmap=&pixelmaps[p->GetId()];
                        
                        
                        double sum_e=0;
                        for(unsigned int k=0;k<pmap->size();k++)
                        {
                                int pid = (*pmap)[k].first;

                                if(printit)printf("\t\t xtalks to pixel %d\n",pid);

                                Pixel * q = pmt->GetPixel(pid,0);
                                if(!q)continue;
                                
                                int dir = (*pmap)[k].second;
                                
                                std::vector<int> psv = q->GetPixelSpotVector();

                                
                                for(unsigned int l=0;l<psv.size();l++)
                                {
                                        PixelSpot *ps = q->GetPixelSpot(psv[l],0);
                                        if(!ps)continue;
                                        
  //                                      double xe = ps->GetE()*_pixelSpotXTalkMap[ps->GetId()-1][8-dir];//reverse the dir

    //                                    sum_e+=xe;

                                        if(printit)printf("\t\t\t pixelspot %d addr %p dir %d with xe %f    e %f  psxtm %f stripE %f\n",ps->GetId(),(void*)ps,8-dir,xe,ps->GetE(), _pixelSpotXTalkMap[ps->GetId()-1][8-dir],ps->GetStripE());
                                }
                        }
                }
*/
///////////////////
        }

        if(printit)printf("found %d strips in %d pmts with total %f energy and %f redistributed\n",xtalkstrips, xtalkpmts,totale,xtalke);

}




//......................................................................

XTalkFilter::XTalkFilter()
{
        if(printit<0)
        {
                if(MsgService::Instance()->IsActive("XTalkFilter",Msg::kDebug))printit=1;
                else printit=0;
        }

	detector = Detector::kUnknown;
        lastpixelmap = Detector::kUnknown;
	method=0;
        Initialize();
}


void XTalkFilter::GetPixelMaps()
{
	if(lastpixelmap != detector)
	{
		pixelmaps.clear();
		if(detector==Detector::kFar)
		{
			for(int i=0;i<16;i++)
        		{
                		pixelmaps.push_back(GetXTalkPixelsFar(i));
        		}
		}
		if(detector==Detector::kNear)
                {
                        for(int i=0;i<64;i++)
                        {
                                pixelmaps.push_back(GetXTalkPixelsNear(i));
                        }
                }


		lastpixelmap=detector;

	}


}


//......................................................................

XTalkFilter::~XTalkFilter()
{}

//......................................................................
void XTalkFilter::BeginJob()
{

}


void XTalkFilter::EndJob()
{
}

JobCResult XTalkFilter::Reco(MomNavigator* mom)
{
   bool foundST=false;
                                                                               
//   VldContext vc;
   //first ask mom for a NtpStRecord

//iterate over all NtpStRecords in MOM

std::vector<TObject*> records = mom->GetFragmentList("NtpStRecord");


for(unsigned int i=0;i<records.size();i++)
{

   NtpStRecord *str=dynamic_cast<NtpStRecord *>(records[i]); //(mom->GetFragment("NtpStRecord","Primary"));
   if(str){
     foundST=true;
  //   vc=str->GetHeader().GetVldContext();
   }else{
        continue;
   }
   //cout<<"Running over "<<str->GetName()<<endl;

        NtpSREventSummary * evthdr = dynamic_cast<NtpSREventSummary *>(&str->evthdr);

        for(int ievt=0;ievt<evthdr->nevent;ievt++)
        {
        
        
/*              for(int j=0;j<str->stp->GetEntries();j++)
                {
                        NtpSRStrip* strip = (NtpSRStrip*)(str->stp->At(j));
                        strip->ph0.pe=5;
                        strip->ph0.sigcor=5;
                    strip->ph0.siglin=5;
                    strip->ph0.raw=5;

                        strip->ph1.pe=5;
                        strip->ph1.sigcor=5;
                    strip->ph1.siglin=5;
                    strip->ph1.raw=5;

                }continue;              
*/      
                if(method==0)DoTrueFilter(str,ievt);    
        }
}
   
   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}




const Registry& XTalkFilter::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("XTalkFilter",Msg::kDebug)<<"In Trimmer::DefaultConfig"<<endl;

  
  static Registry r;
  std::string name = this->JobCModule::GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  r.UnLockValues();
  r.Set("Method",0);         

  r.LockValues();             
               
                                                                              
  return r;
}

void XTalkFilter::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("XTalkFilter",Msg::kDebug)<<"In Trimmer::Config"<<endl;
                                                                               
  bool islocked = this->GetConfig().ValuesLocked();
  if (islocked) this->GetConfig().UnLockValues();
 
  int valint;
  if(r.Get("Method", valint)) 
  {
    method=valint;
  }

  if (islocked) this->GetConfig().LockValues();

        MSG("XTalkFilter",Msg::kDebug)<< "Using method "<<method<<endl;

}



//----------------------------------------------------------------------------
// PMTmuxID is the internal number scheme used by Mark Thomson in the AltDeMux
//---------------------------------------------------------------------------- 
Int_t XTalkFilter::CalculatePMTmuxID(PlexPixelSpotId &ppsid){

  Int_t elvalue=0;
  elvalue = ppsid.GetTube() + 3*ppsid.GetInRack() + 24*ppsid.GetRackBay();
  
  Char_t level = ppsid.GetRackLevel(); 
  Char_t ew = ppsid.GetEastWest(); 
  
  if(level=='U'){
    elvalue += 450;
  }
  if(ew =='E'){
    elvalue += 900;
  }
  
  return elvalue;

}
