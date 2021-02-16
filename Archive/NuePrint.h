/**
 *
 * $Id: NuePrint.h,v 1.1 2006/09/18 21:41:22 boehm Exp $
 *
 * \class NuePrint
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: M. Sanchez
 *
 * Created on: Thu Apr 21 19:49:07 2005
 *
 */

#ifndef NUEPRINT_H
#define NUEPRINT_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif


//C++ headers
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

class NueRecord;

class NuePrint : public JobCModule
{

public:

    NuePrint();
    ~NuePrint();

    // Analysis and Reconstruction methods
    JobCResult Ana(const MomNavigator* mom);
    void BeginJob();
    void EndJob();

    // Module configuration
    const Registry& DefaultConfig() const;
    void            Config(const Registry& r);  
    
    // Methods to print stuff
    TList *GetVarList(NueRecord *nr);
    void PrintHeader(TList* vlist);
    void PrintValues(NueRecord *nr, std::string &sr, std::string &defval);
    Bool_t IsVarExcluded(std::string &sr);
    Bool_t IsVarIncluded(std::string &sr);
    std::string GetClassLabel(const char &str);
    void FillIncludeVec();

private:
    ofstream *fileOut;
    int counter;
    int passcounter;  
    Int_t nvar;
    std::string outFilename;
    std::string outFormat;
    std::string excludeVars;
    Int_t kHiPlaneTrackCut;
    Float_t kPhProngCut;
    Float_t kMeuEnergyCut;
    Int_t kLoPlaneEventCut;
    Int_t kHiTrackLikeCut;
    Int_t kLoPhNStripCut;
    Int_t kLoPhNPlaneCut;
    Float_t kHiEnergyCut;
    Float_t kLoEnergyCut;
    Float_t kHiEnergyShowerCut;
    Float_t kLoEnergyShowerCut;
    Int_t kSigClass;
    Int_t kBgNCClass;
    Int_t kBgNumuClass;
    Int_t kBgNutauClass;
    Int_t kBgBnueClass;
    Int_t kSigResCode;
    Float_t kEmFrac;

    Int_t outAll;
    Int_t training;
    Int_t recalc;
    Int_t deftrk;

    Float_t kDM2;
    Float_t kTheta23;
    Float_t kUE32;

    std::string incFile;
    std::vector<std::string> excludeVec;
    std::vector<std::string> includeVec;

};                              // end of class NuePrint

#endif  // NUEPRINT_H
