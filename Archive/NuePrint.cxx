
/**
 *
 * $Id: NuePrint.cxx,v 1.2 2008/11/19 18:22:51 rhatcher Exp $
 *
 * \class NuePrint
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: M. Sanchez
 *
 * Created on: Thu Apr 21 19:52:28 2005
 *
 */

#include "TMath.h"
#include "NueAna/NuePrint.h"
#include "NueAna/NueRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "TClass.h"
#include "TDataType.h"
#include "TDataMember.h"
#include "TRealData.h"
#include "ClassType.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "Util/UtilString.h"
#include "NueAna/NueRWHelpers.h"

JOBMODULE(NuePrint, "NuePrint",
          "Prints ana_nue ntuples for weka/sas analysis");
CVSID("$Id:");

NuePrint::NuePrint():
    counter(0),
    passcounter(0),
    outFilename(""),
    outFormat(""),
    excludeVars(""),
    kHiPlaneTrackCut(25),
    kPhProngCut(-1),
    kMeuEnergyCut(-1),
    kLoPlaneEventCut(-1),
    kHiTrackLikeCut(-1),
    kLoPhNStripCut(-1),
    kLoPhNPlaneCut(-1),
    kHiEnergyCut(-1),
    kLoEnergyCut(-1),
    kHiEnergyShowerCut(-1),
    kLoEnergyShowerCut(-1),
    kSigClass(2),
    kBgNCClass(0),
    kBgNumuClass(1),
    kBgNutauClass(3),
    kBgBnueClass(4),
    kSigResCode(1004),
    kEmFrac(-1),
    outAll(0),
    training(0),
    recalc(0),
    deftrk(0),
    kDM2(0.0025),
    kTheta23(TMath::Pi()/4.0),
    kUE32(0.01),
    incFile("")
{

}

NuePrint::~NuePrint() 
{

}
JobCResult NuePrint::Ana(const MomNavigator* mom)
{


   TObject *obj=0;
   TIter objiter = mom->FragmentIter();

   TList *varList;
   string sepVal, defVal="";
   if(outFormat=="WEKA"){ sepVal=","; defVal="?";}
   else if(outFormat=="SAS"){ sepVal=" "; defVal="-9999.99";}
   else if(outFormat=="SPR"){ sepVal=" "; defVal="-9999.99";}

   while((obj=objiter.Next())){


      NueRecord *nr = static_cast<NueRecord *>(obj);


      if(nr){
	 MSG("NuePrint",Msg::kDebug)<<"Found a NueRecord in MOM"<<endl;
      }
      else{
	 MSG("NuePrint",Msg::kError)<<"Didn't find a NueRecord in MOM"<<endl;
	 counter++;
	 continue;
      }
      MSG("NuePrint",Msg::kDebug)<<"Found a NueRecord "<<nr<<endl;

       if(!counter){
           varList=GetVarList(nr);
           PrintHeader(varList);
 
       }
      counter++;

      MSG("NuePrint",Msg::kDebug) << nr->GetHeader().GetRun() << " " 
                                  << nr->GetHeader().GetSnarl() << " " 
                                  << nr->GetHeader().GetEventNo() << " " 
                                  << nr->GetHeader().GetEvents() 
                                 << " " << outAll << endl;
      

      // Add any desired cuts here
      // Cut for no reco event
      if (nr->GetHeader().GetEventNo()<0&&!outAll) continue;
      // Cut for Fiducial volume
      if (nr->anainfo.inFiducialVolume != 1&&!outAll) continue;
      // Cut for Full Containment, note difference in Near/Far
      if (nr->GetHeader().GetVldContext().GetDetector()== Detector::kFar
          && nr->anainfo.isFullyContained != 1&&!outAll) continue;
      if (nr->GetHeader().GetVldContext().GetDetector()== Detector::kNear
          && nr->anainfo.isFullyContained != 1 
          && nr->anainfo.isFullyContained != -2&&!outAll) continue;
      // Cut for max number of planes in Track
      if (nr->srtrack.planes>=kHiPlaneTrackCut&&kHiPlaneTrackCut>=0
          &&!outAll) continue;
      // Cut on max total event energy in Meu
      if(nr->srevent.phMeu>kMeuEnergyCut&&kMeuEnergyCut>=0&&!outAll) continue;
      // Cut on min Total pulse height per prong (sigcor)
      if(TMath::Max(nr->srtrack.pulseHeight,nr->srshower.pulseHeight)
	 <kPhProngCut&&kPhProngCut>=0&&!outAll) continue;
      // Cut on max number of trklike planes in Track
      if (nr->srtrack.trklikePlanes>=kHiTrackLikeCut&&kHiTrackLikeCut>=0
          &&!outAll) continue;
      // Cut on min number of planes in Event
      if (nr->srevent.planes<kLoPlaneEventCut&&kLoPlaneEventCut>=0
          &&!outAll) continue;
      // Cut on max total event energy in Meu
      if (nr->srevent.phMeu>=kHiEnergyCut
	  &&kHiEnergyCut>=0&&!outAll) continue;
      // Cut on min total event energy in Meu
      if (nr->srevent.phMeu<kLoEnergyCut
	  &&kLoEnergyCut>=0&&!outAll) continue;
      // Cut on max shower energy (in gev) 
      if (nr->srevent.phNueGeV>=kHiEnergyShowerCut
	  &&kHiEnergyShowerCut>=0&&!outAll) continue;
      // Cut on min shower energy (in gev) 
      if (nr->srevent.phNueGeV<kLoEnergyShowerCut
	  &&kLoEnergyShowerCut>=0&&!outAll) continue;
      // Cut on min number of strips with above a threshold ph
      if (nr->shwfit.hiPhStripCount<kLoPhNStripCut
	  &&kLoPhNStripCut>=0&&!outAll) continue;
      // Cut on min number of planes with above a threshold ph
      if (nr->shwfit.hiPhPlaneCount<kLoPhNPlaneCut
	  &&kLoPhNPlaneCut>=0&&!outAll) continue;
      // Cuts for signal/background selection
      if(training==1&&(nr->mctrue.fNueClass==kSigClass
		       ||nr->mctrue.fNueClass==kBgNCClass 
		       ||nr->mctrue.fNueClass==kBgNumuClass
		       ||nr->mctrue.fNueClass==kBgNutauClass  
		       ||nr->mctrue.fNueClass==kBgBnueClass)) { 
	MSG("NuePrint",Msg::kDebug) << "matched class" <<endl;}
      else if (training==0) { 
	MSG("NuePrint",Msg::kDebug) << "training = 0, nothing to do " <<endl;}
      else{ continue;} 
      if(nr->mctrue.fNueClass==kSigClass
	 &&nr->mctrue.resonanceCode>kSigResCode) continue;
      if(nr->mctrue.fNueClass==kSigClass
	 &&nr->mctrue.emShowerFraction<=kEmFrac) continue;
      // End of cuts section

      PrintValues(nr,sepVal,defVal);
      
      // here is where we are going to print
      //hpar0->Fill(nr->shwfit.par_a);
      //*fileOut << nr->shwfit.par_a << endl;

      passcounter++;

   }
   return JobCResult::kPassed; // kNoDecision, kFailed, etc.


}
void NuePrint::BeginJob()
{
    MSG("NuePrint",Msg::kDebug)<<"In NuePrint::BeginJob"<<endl; 
    MSG("NuePrint",Msg::kInfo)<<"Writing out to file: "
                              << outFilename << endl;

    fileOut =  new ofstream (outFilename.c_str());

    string str=",";

    MSG("NuePrint",Msg::kInfo) << "Excluding Vars: " << excludeVars <<endl;
    if(excludeVars.find(",")==string::npos) {
        excludeVec.push_back(excludeVars);
    } else { UtilString::StringTok(excludeVec,excludeVars,str);}

    MSG("NuePrint",Msg::kInfo) << excludeVec.size() 
                               << " var: " ;
    for(int ivar=0; ivar < (int) excludeVec.size(); ivar++){ 
        MSG("NuePrint",Msg::kInfo) << excludeVec[ivar] << " " ; 
    } 
    MSG("NuePrint",Msg::kInfo) << endl; 

    if(incFile.size()) FillIncludeVec();

    MSG("NuePrint",Msg::kInfo) << "Printing all events: " << outAll << endl; 
    MSG("NuePrint",Msg::kInfo) << "Hiding dummy vars: " << training << endl; 
    MSG("NuePrint",Msg::kInfo) << "Using Cuts " << endl; 
    if(kHiPlaneTrackCut>=0)
    MSG("NuePrint",Msg::kInfo) << "HiPlaneTrackCut: " << kHiPlaneTrackCut << endl; 
    if(kPhProngCut>=0)
    MSG("NuePrint",Msg::kInfo) << "PhProngCut: " << kPhProngCut << endl; 
    if(kMeuEnergyCut>=0)
    MSG("NuePrint",Msg::kInfo) << "MeuEnergyCut: " << kMeuEnergyCut << endl; 
    if(kHiTrackLikeCut>=0)
    MSG("NuePrint",Msg::kInfo) << "HiTrackLikeCut: " << kHiTrackLikeCut << endl; 
    if(kLoPlaneEventCut>=0)
    MSG("NuePrint",Msg::kInfo) << "LoPlaneEventCut: " << kLoPlaneEventCut << endl; 
    if(kHiEnergyCut>=0)
    MSG("NuePrint",Msg::kInfo) << "HiEnergyCut: " << kHiEnergyCut << endl; 
    if(kLoEnergyCut>=0)
    MSG("NuePrint",Msg::kInfo) << "LoEnergyCut: " << kLoEnergyCut << endl; 
    if(kHiEnergyShowerCut>=0)
    MSG("NuePrint",Msg::kInfo) << "HiEnergyShowerCut: " << kHiEnergyShowerCut << endl; 
    if(kLoEnergyShowerCut>=0)
    MSG("NuePrint",Msg::kInfo) << "LoEnergyShowerCut: " << kLoEnergyShowerCut << endl; 
    if(kLoPhNStripCut>=0)
    MSG("NuePrint",Msg::kInfo) << "LoPhNStripCut: " << kLoPhNStripCut << endl; 
    if(kLoPhNPlaneCut>=0)
    MSG("NuePrint",Msg::kInfo) << "LoPhNPlaneCut: " << kLoPhNPlaneCut << endl; 




}
void NuePrint::FillIncludeVec()
{
    std::string dumName, dumtype;
    Float_t dumstart, dumend;
    Int_t dumbins, ivar;
    ifstream ins;
    ins.open(incFile.c_str());

    ivar = 0;

    //read in the file
    while(!ins.eof()) {
        ins>>dumName>>dumstart>>dumend>>dumbins>>dumtype;
        if(!ins.eof()){
            includeVec.push_back(dumName);
//            beg.push_back(dumstart);  end.push_back(dumend);
//            nbins.push_back(dumbins);   gtype.push_back(dumtype);
            ivar++;
        }
    }                              

    MSG("NuePrint",Msg::kInfo) << "Using include file: "<< incFile << endl;     
    MSG("NuePrint",Msg::kInfo) << ivar<<" variables read in"<<endl; 

}
void NuePrint::EndJob()
{
   MSG("NuePrint",Msg::kInfo)<<"Counter "<<counter<<" passcounter "<<passcounter<<endl;
   
}
TList * NuePrint::GetVarList(NueRecord *nr){
    TClass *cl;
    TList *vlist;

    cl=nr->IsA();
    if (!cl->GetListOfRealData()) cl->BuildRealData(nr);
    
    vlist=cl->GetListOfRealData();
    MSG("NuePrint",Msg::kInfo)<<"Total Variables being in NueRecord "
                              << vlist->GetSize()  << endl;         

    return vlist;
}

void NuePrint::PrintHeader(TList * vlist){
    string varName;
    TRealData *rd;

    TIter next(vlist);
    Int_t counter=0;
    MSG("NuePrint",Msg::kDebug)<< "Attributes list: ";

    if(outFormat=="WEKA")
        *fileOut << "@relation nue_output" << endl << endl;

    if(outFormat=="SPR"){
      *fileOut << "# SPR file definition" << endl; 
      *fileOut << "nvar"<< endl;
    }    
    while((rd =dynamic_cast<TRealData*>(next()))){
        varName=rd->GetName();
        if(!IsVarExcluded(varName)){
            MSG("NuePrint",Msg::kDebug)<< varName << ",";
            if(outFormat=="WEKA"){
            *fileOut << "@attribute " << varName << " real" << endl;    
            }else if(outFormat=="SAS"){
                *fileOut << varName << ",";
            }else if(outFormat=="SPR"){
                *fileOut << varName << " ";
            }
            counter++;
	}
    }

    
    if(outFormat=="WEKA"){
        *fileOut << "@attribute fHeader.fRun real" <<endl;
        *fileOut << "@attribute fHeader.fSubRun real" <<endl;
        *fileOut << "@attribute fHeader.fSnarl real" <<endl;
        *fileOut << "@attribute fHeader.fEvtNo real" <<endl;
    }
    else if(outFormat=="SAS"){
     *fileOut << "fHeader.fRun,fHeader.fSubRun,fHeader.fSnarl,fHeader.fEvtNo,";
    }
    else if(outFormat=="SPR"){
     *fileOut << "fHeader.fRun fHeader.fSubRun fHeader.fSnarl fHeader.fEvtNo ";
    }

    if(outFormat=="WEKA"){
        *fileOut << "@attribute class {nc,numu,nue,nutau,bnue,empty}" << endl;
    }else if(outFormat=="SAS"){
      *fileOut << "fHeader.fFileNo" << "," << "weight" << ","<< "class" << endl;
    }else if(outFormat=="SPR"){
        *fileOut << "fHeader.fFileNo" << endl;
    }
    MSG("NuePrint",Msg::kDebug)<< "class" << endl;
    
    if(outFormat=="WEKA"){ counter=counter+5; nvar=counter;}
    if(outFormat=="SAS"){ counter=counter+7; nvar=counter;}
    if(outFormat=="SPR"){ counter=counter+7; nvar=counter;}

    MSG("NuePrint",Msg::kInfo)<< "Total vars being used: " << nvar << endl;
    if(outFormat=="WEKA")
        *fileOut << endl << "@data" << endl;
}

Bool_t NuePrint::IsVarExcluded(std::string &sr){

    if(sr=="mctrue.fNueClass") return true;
    if(sr=="fHeader.fRun") return true;
    if(sr=="fHeader.fSubRun") return true;
    if(sr=="fHeader.fSnarl") return true;
    if(sr=="fHeader.fEvtNo") return true;

    if(sr.find(".")==string::npos) return true;

    for(int ivar=0; ivar < (int) excludeVec.size(); ivar++){
//        cout << ivar << " " << excludeVec[ivar]<< endl;
        if(sr.find(excludeVec[ivar])!=string::npos) return true;
    }

    if(sr.find("*")!=string::npos) return true;
    if(sr.find("fUniqueID")!=string::npos) return true;
    if(sr.find("fBits")!=string::npos) return true;
    if(sr.find("fTempTags")!=string::npos) return true;
    if(sr.find("fData")!=string::npos) return true;
    if(sr.find("fVldContext")!=string::npos) return true;
    if(sr.find("found")!=string::npos) return true;
    
    if(includeVec.size()>0){
        if(IsVarIncluded(sr)) return false;
        else return true;
    }

    return false;
}
Bool_t NuePrint::IsVarIncluded(std::string &sr){

    // Note exclusion precedes inclusion, so if you exclude a branch
    // and then include a var in your file is not going to get printed

    for(int ivar=0; ivar < (int) includeVec.size(); ivar++){
//        cout << ivar << " " << includeVec[ivar]<< endl;
        if(sr == includeVec[ivar]) return true;
    }
    return false;
}
void NuePrint::PrintValues(NueRecord *nr, std::string &sr, std::string &defval){
   TClass *cl;
   TRealData *rd;
   string varName;
   TDataMember *member;
   TDataType *membertype;
   Float_t value;
   Int_t fileNo=0;
   string classStr;
   const char *buffer;
   Int_t valrun=0,valsub=0,valsnr=0,valevt=0;
   Float_t valweight = 0.;
   Int_t ivar=0;
   Int_t pvar=0;

   cl=nr->IsA();      
   TIter  next(cl->GetListOfRealData());

   while ((rd =dynamic_cast<TRealData*>(next()))) {
       member = rd->GetDataMember();
       membertype = member->GetDataType();
       varName=rd->GetName();

       Int_t offset = rd->GetThisOffset();
       char *pointer = (char*)nr  + offset;

       if (!IsVarExcluded(varName)) {
//           if (varName.find("fHeader.fRun")!=string::npos){

//           }
           value=atof(membertype->AsString(pointer));
           ivar++;
           MSG("NuePrint",Msg::kDebug)<<member->GetFullTypeName() << " " << varName << " " << value << endl;

           if(strcmp(member->GetTypeName(),"Float_t")==0
              || strcmp(member->GetTypeName(),"float")==0
              || strcmp(member->GetTypeName(),"Double_t")==0
              || strcmp(member->GetTypeName(),"double")==0){
               if ((ANtpDefVal::IsDefault(value)
                    || ANtpDefVal::IsDefault(static_cast<Double_t> (value))
                    || !TMath::Finite(value) || TMath::IsNaN(value) || value>1e+8 )
                   &&!defval.empty()){
                   MSG("NuePrint",Msg::kDebug) << "Found default value" << endl;
                   if(deftrk&&varName.find("srtrack")!=string::npos
                      &&varName!="srtrack.rangeMomentum"
                      &&varName!="srtrack.sigmaQoverP"){
                       *fileOut << deftrk << sr ; pvar++;
                   }
                   else{*fileOut << defval << sr ; pvar++;}
               }
               else if((value>-1e-8&&value<1e-8)&&!varName.find("timeLength")){
                           value=0.0;
                           *fileOut << value << sr ; pvar++; 
               } else {
		 valweight=1;
                   if(recalc==1 
                      && varName.find("mctrue.fOscProb")!=string::npos){
                       MSG("NuePrint",Msg::kDebug)<<" Value in ntuple " 
                                                  << value << endl;
                       value=NueRWHelpers::Oscillate(nr->mctrue.nuFlavor
                                                     ,nr->mctrue.nonOscNuFlavor
                                                     ,nr->mctrue.nonOscNuEnergy
                                                     ,735, kDM2, kTheta23
                                                     , kUE32);
                       valweight=value;
                       MSG("NuePrint",Msg::kDebug)<<" Value of " << varName 
                                                  <<" after recalc " << value 
                                                  <<" from "
                                                  << nr->mctrue.nuFlavor 
                                                  << " " 
                                                  << nr->mctrue.nonOscNuFlavor 
                                                  << " "
                                                  << nr->mctrue.nonOscNuEnergy 
                                                  << " "
                                                  << endl;
                      
                   }                    
                   if(deftrk&&varName.find("srtrack")!=string::npos
                      &&varName!="srtrack.rangeMomentum"
                      &&varName!="srtrack.sigmaQoverP"){
                       *fileOut << deftrk << sr ; pvar++;
                   }
                   else{*fileOut << value << sr ; pvar++;}
               }
           } else if(strcmp(member->GetTypeName(),"Int_t")==0
                     || strcmp(member->GetTypeName(),"int")==0){
               if ((ANtpDefVal::IsDefault(static_cast<Int_t> (value))
                    || ANtpDefVal::IsDefault(static_cast<Int_t> (value))
                    || !TMath::Finite(value) || TMath::IsNaN(value))
                   &&!defval.empty()){
                   MSG("NuePrint",Msg::kDebug) << "Found default value" << endl;
                   *fileOut << defval << sr ; pvar++;
               }
               else { 
                   *fileOut << value << sr ; pvar++;
               }
           }
           else{
               MSG("NuePrint",Msg::kError)
                   << "Problem: you are not taking into account this type: " 
                   <<"Membertype "<< member->GetFullTypeName() 
                   << " VarName " << varName  
                   << " Value " <<   membertype->AsString(pointer) << endl;  
               
           }
           MSG("NuePrint",Msg::kDebug)
               <<"Membertype "<< member->GetFullTypeName()
               << " VarName " << varName 
               << " Value " <<   membertype->AsString(pointer) << endl;     
       }
       else if(varName.find("NueClass")!=string::npos){
           ivar++;
           buffer=membertype->AsString(pointer);
           classStr=GetClassLabel(*buffer);
       }
       else if(varName.find("fHeader.fRun")!=string::npos){
           ivar++;
           valrun=atoi(membertype->AsString(pointer));
           if (outFormat=="SAS" || outFormat=="SPR"){
               fileNo=valrun%1000;
               ivar++;
           }
       }
       else if(varName.find("fHeader.fSubRun")!=string::npos){
           ivar++;
           valsub=atoi(membertype->AsString(pointer));
       }
       else if(varName.find("fHeader.fSnarl")!=string::npos){
           ivar++;
           valsnr=atoi(membertype->AsString(pointer));
       }
       else if(varName.find("fHeader.fEvtNo")!=string::npos){
           ivar++;
           valevt=atoi(membertype->AsString(pointer));
       } 
   }

   if(training){
       *fileOut << "99999999" << sr << "9" << sr << "9999" << sr << "9" << sr;
       pvar=pvar+4;
   }
   else{
     *fileOut << valrun << sr << valsub << sr << valsnr << sr << valevt << sr; 
     pvar=pvar+4; 
   }
   if (outFormat=="SAS" || outFormat=="SPR"){ 
     if(!training) {*fileOut << fileNo << sr; pvar++;}
     else {*fileOut << 9 << sr; pvar++;}
     if(recalc==1){*fileOut << valweight << sr; pvar++; ivar++;}
     else{*fileOut << 1 << sr; pvar++; ivar++;}
   }

   *fileOut << classStr; pvar++;
   *fileOut << endl; 

   if(ivar!=nvar||pvar!=nvar) { 

       MSG("NuePrint",Msg::kError) << " printing less vars than in header" 
                                   <<endl;

       MSG("NuePrint",Msg::kError) << "Header: " << nvar 
                               << " Not excluded: " << ivar
                               << " Printed: " << pvar << endl;

   }

}
std::string NuePrint::GetClassLabel(const char &str){
    string cl;
    if(outFormat=="WEKA"){
        if(atoi(&str)==ClassType::NC) cl="nc";
        else if(atoi(&str)==ClassType::numu) cl="numu";
        else if(atoi(&str)==ClassType::nue) cl="nue";
        else if(atoi(&str)==ClassType::nutau) cl="nutau";
        else if(atoi(&str)==ClassType::bnue) cl="bnue";
        else {
            MSG("NuePrint",Msg::kWarning) << "Event with no class" << endl;
            cl="empty";
        } 
    }else if(outFormat=="SAS"){
        if(atoi(&str)==ClassType::NC) cl="ncu";
        else if(atoi(&str)==ClassType::numu) cl="num";
        else if(atoi(&str)==ClassType::nue) cl="nue";
        else if(atoi(&str)==ClassType::nutau) cl="nut";
        else if(atoi(&str)==ClassType::bnue) cl="bnu";
        else {
            MSG("NuePrint",Msg::kWarning) << "Event with no class" << endl;
            cl="non";
        }
    }else if(outFormat=="SPR"&&training){
        if(atoi(&str)==ClassType::NC) cl="0";
        else if(atoi(&str)==ClassType::numu) cl="0";
        else if(atoi(&str)==ClassType::nue) cl="1";
        else if(atoi(&str)==ClassType::nutau) cl="0";
        else if(atoi(&str)==ClassType::bnue) cl="0";
        else {
            MSG("NuePrint",Msg::kWarning) << "Event with no class" << endl;
            cl="0";
        }
    }else if(outFormat=="SPR"&&!training){
        if(atoi(&str)==ClassType::NC) cl="0";
        else if(atoi(&str)==ClassType::numu) cl="1";
        else if(atoi(&str)==ClassType::nue) cl="2";
        else if(atoi(&str)==ClassType::nutau) cl="3";
        else if(atoi(&str)==ClassType::bnue) cl="4";
        else {
            MSG("NuePrint",Msg::kWarning) << "Event with no class" << endl;
            cl="5";
        }
    }            
    
    return cl;
    
}
const Registry& NuePrint::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
    MSG("NuePrint",Msg::kDebug)<<"In NuePrint::DefaultConfig"<<endl;

    static Registry r; // Default configuration for module

    // Set name of config
    std::string name = this->GetName();
    name += ".config.default";
    r.SetName(name.c_str());

    // Set values in configuration
    r.UnLockValues();
    r.Set("excludeVars","templates.root");
    r.Set("outFilename","training.arff");
    r.Set("outFormat","WEKA");      
    r.Set("includeFile","");


    r.Set("HiPlaneTrackCut",25);   
    r.Set("printAll","0");      
    r.Set("training","0");      
    r.Set("recalcOsc",0);
    r.Set("defaultTrkVal",0);

    r.Set("LoPlaneEventCut",-1);
    r.Set("HiTrackLikeCut",-1);
    r.Set("LoPhNStripCut",-1);
    r.Set("LoPhNPlaneCut",-1);
    
    r.Set("MeuEnergyCut",0.);
    r.Set("PhProngCut",0.);
    
    r.Set("HiEnergyCut",-1);
    r.Set("LoEnergyCut",-1);
    r.Set("HiEnergyShowerCut",-1);
    r.Set("LoEnergyShowerCut",-1);

    r.Set("SigClass",2);
    r.Set("BgNCClass",0);
    r.Set("BgNumuClass",1);
    r.Set("BgNutauClass",3);
    r.Set("BgBnueClass",4);
    r.Set("SigResCode",1004);
    r.Set("EmFrac",-1);


    r.Set("DM2", 0.0025);   
    r.Set("Theta23", TMath::Pi()/4.0);
    r.Set("UE32", 0.01);

    r.LockValues();

    return r; 

}
void NuePrint::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
    MSG("NuePrint",Msg::kDebug)<<"In NuePrint::Config"<<endl;

    const char* tmps;


    if(r.Get("excludeVars",tmps)) { excludeVars = tmps;}
    if(r.Get("outFilename",tmps)) { outFilename = tmps;}
    if(r.Get("outFormat",tmps)) { outFormat = tmps;}
    if(r.Get("includeFile",tmps)) {incFile = tmps;}



    int imps;
    if(r.Get("HiPlaneTrackCut",imps)) { kHiPlaneTrackCut=imps;}
    if(r.Get("printAll",imps)) { outAll = imps;}
    if(r.Get("training",imps)) { training = imps;}
    if(r.Get("recalcOsc",imps)) { recalc = imps;}
    if(r.Get("defaultTrkVal",imps)) {deftrk = imps;}

    if(r.Get("LoPlaneEventCut",imps)) { kLoPlaneEventCut=imps;}
    if(r.Get("HiTrackLikeCut",imps)) { kHiTrackLikeCut=imps;}
    if(r.Get("LoPhNStripCut",imps)) { kLoPhNStripCut=imps;}
    if(r.Get("LoPhNPlaneCut",imps)) { kLoPhNPlaneCut=imps;}

    if(r.Get("SigClass",imps)){kSigClass=imps;}
    if(r.Get("BgNCClass",imps)){kBgNCClass=imps;}
    if(r.Get("BgNumuClass",imps)){kBgNumuClass=imps;}
    if(r.Get("BgNutauClass",imps)){kBgNutauClass=imps;}
    if(r.Get("BgBnueClass",imps)){kBgBnueClass=imps;}
    if(r.Get("SigResCode",imps)){kSigResCode=imps;}

    double fmps;
    if(r.Get("MeuEnergyCut",fmps)){ kMeuEnergyCut=fmps; }                    
    if(r.Get("PhProngCut",fmps)){ kPhProngCut=fmps; }                    

    if(r.Get("HiEnergyCut",fmps)) { kHiEnergyCut=fmps;}
    if(r.Get("LoEnergyCut",fmps)) { kLoEnergyCut=fmps;}
    if(r.Get("HiEnergyShowerCut",fmps)) { kHiEnergyShowerCut=fmps;}
    if(r.Get("LoEnergyShowerCut",fmps)) { kLoEnergyShowerCut=fmps;}

    if(r.Get("EmFrac",fmps)){kEmFrac=fmps;}

    if(r.Get("DM2",fmps)) { kDM2 = fmps;}
    if(r.Get("Theta23",fmps)) { kTheta23 = fmps;}
    if(r.Get("UE32",fmps)) { kUE32 = fmps;}

}
