/// $Id: MdaDiscrimAna.cxx,v 1.9 2007/03/01 16:38:49 rhatcher Exp $ 
///
/// class MdaDiscrimAna
///
/// NueAna package
///
/// Purpose: Calculate PID from SAS discriminant functions and perform 
///          posterior event classification. 
///
/// Author: Alex Sousa <a.sousa@physics.ox.ac.uk>
///
/// Created on:  Fri Feb 17 15:50:47 2006

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;

#include <string>
using std::string;

#include "MessageService/MsgService.h"
#include "NueAna/MdaDiscrimAna.h"
#include "NueAna/StringUtil.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

#include "TList.h"
#include "TClass.h"
#include "TDataType.h"
#include "TDataMember.h"
#include "TRealData.h"

ClassImp(MdaDiscrimAna)

CVSID("$Id: MdaDiscrimAna.cxx,v 1.9 2007/03/01 16:38:49 rhatcher Exp $");

Int_t MdaDiscrimAna::nClass;
std::deque <TMatrixD> MdaDiscrimAna::quadCoeff;
std::deque <TVectorD> MdaDiscrimAna::linCoeff;
std::deque <TVectorD> MdaDiscrimAna::constCoeff;
std::deque <TVectorD> MdaDiscrimAna::meanVec;   

std::deque < std::deque <std::string> > MdaDiscrimAna::typeClassVec;
std::deque < std::string> MdaDiscrimAna::typeCoeffVec;
std::deque < std::string> MdaDiscrimAna::varNameVec;
    
MdaDiscrimAna::MdaDiscrimAna(NueRecord &nr, MdaDiscrim &md):
    nueRec(nr),
    fMdaDiscrim(md)
{

}

MdaDiscrimAna::~MdaDiscrimAna()
{

}


void MdaDiscrimAna::Analyze()
{

    MSG("MdaDiscrimAna",Msg::kDebug)<<"In MdaDiscrimAna::Analyze"<<endl;

    //Fill SAS Calibration Arrays
    if (!isFillDone()) isFilled = FillCalibArrays();
    
    //Fill P-Vector values from NueRecord
    FillPVector();
    
    //Calculate nue PID from discriminator function
    PIDCalc();
    
    //Classify events according to a probability threshold cut  
    MdaClassify();
            
}

Bool_t MdaDiscrimAna::FillCalibArrays()
{
    
    //Reads the SAS coefficient input file 
    //(containing quadratic, linear and constant coeffs)
    //and fills the necessary matrices and vectors
    //Also fills a vector 

    static const Int_t kLineSize=2000;

    ifstream sasInput;
    sasInput.open(sasFileName.c_str());
    if (sasInput.fail()) {
        
        MSG("MdaDiscrimAna",Msg::kInfo)<< "SAS coefficient file not found! " 
                                        << "Not calculating MDA PID..." 
                                        << endl;
        return false;
    }

    Char_t oneLine[kLineSize];

    string sepaCol=" ";
    vector<string> oneVec;
         
    //Read the first line (mean values) and get the number 
    //of columns in the file
    sasInput.getline(oneLine,kLineSize);
    string topLine=oneLine;
    StringUtil::SplitString(topLine,sepaCol,oneVec);
    
    //Read the mean values into a vector
    TVectorD meanTemp(oneVec.size()-2);
    for(UInt_t col=2; col < oneVec.size(); col++){
        meanTemp(col-2)=(atof((oneVec.at(col)).c_str()));
    }
    
    for(Int_t m=0; m < meanTemp.GetNoElements(); m++){
        MSG("MdaDiscrimAna",Msg::kDebug)<< "Mean:" 
                                        << m 
                                        << "="
                                        << meanTemp(m) 
                                        << endl; 
    }

    meanVec.push_back(meanTemp);

    //Declare matrices and vectors for coefficient reading
    TMatrixD quadTemp(oneVec.size()-2, oneVec.size()-2);
    
    TVectorD linTemp(oneVec.size()-2);
    
    TVectorD constTemp(oneVec.size()-2);

    deque <string> tempClassVec;

    //Count the rows for each matrix
    Int_t rowCounter=0;
    
    //Reset the number of classes counter
    nClass=0;

    //Reset the ifstream pointer to the beginning of the file
    // sasInput.seekg(0,ios::beg);
    
    //Reset the line vector
    oneVec.clear();
    
    while (! sasInput.eof() ){
        sasInput.getline(oneLine,kLineSize);
        topLine=oneLine;
        StringUtil::SplitString(topLine,sepaCol,oneVec);
        
        string coeffType;            
        
        for(UInt_t col=0; col < oneVec.size(); col++){
            
            if (col==0) {
            
                tempClassVec.push_back(oneVec.at(col));
                
            }else if (col==1){ 
                
                coeffType=oneVec.at(col);
                
                if(coeffType!="_LINEAR_" && coeffType !="_CONST_"){
                    if (varNameVec.size()<(oneVec.size()-2)) 
                        varNameVec.push_back(oneVec.at(col));
                    typeCoeffVec.push_back("_QUAD_");
                }else typeCoeffVec.push_back(oneVec.at(col));
             
            }else{
             
                if (coeffType=="_LINEAR_"){
                    
                    linTemp(col-2)=(atof((oneVec.at(col)).c_str()));
                 
                }else if (coeffType=="_CONST_"){
                
                    constTemp(col-2)=(atof((oneVec.at(col)).c_str()));
                    
                }else{
                    quadTemp(rowCounter,(col-2))
                        =(atof((oneVec.at(col)).c_str()));
                }
            }
        }
        oneVec.clear();
        rowCounter++;
        
        if (coeffType=="_CONST_"){
            nClass++;
            rowCounter=0;
            typeClassVec.push_back(tempClassVec);
            quadCoeff.push_back(quadTemp);
            linCoeff.push_back(linTemp);
            constCoeff.push_back(constTemp);            
        
            tempClassVec.clear();

        }else continue;
    }
    
    if(sasInput) sasInput.close();

    return true;
}

void MdaDiscrimAna::FillPVector()
{
    //Now use varNameVec to fill the variable values p-vector
    //Shameless rip-off of CompareAll::FillFromList

    NueRecord *nueR=&nueRec;
    
    if(varNameVec.size() == 0) return;
    TString hname;
    UInt_t count = 0;
    
    TClass *cl;
    TRealData *rd;
    string vName;
    TDataMember *member;
    TDataType *membertype;
    Double_t value = 0.0;
    Int_t   valueI = 0;
    
    cl=nueR->IsA();
    TIter  next(cl->GetListOfRealData());                                                                                 
    TVectorD valueTemp(varNameVec.size());

    while ((rd =dynamic_cast<TRealData*>(next()))) {
        member = rd->GetDataMember();
        membertype = member->GetDataType();
        vName=rd->GetName();
        
        TString vNameOrig=vName;

        //Replace "." with "_" in the ntuple variable names
        if (vName.length() > vName.find_first_of("."))
            vName.replace(vName.find_first_of("."),1,"_");
        
        Int_t offset = rd->GetThisOffset();
        Char_t *pointer = reinterpret_cast<Char_t*>(nueR)  + offset;
        
        for(UInt_t i = 0; i < varNameVec.size();i++){
            MSG("MdaDiscrimAna",Msg::kDebug)<<"Found variable "
                                            << "NueRec: " << vName 
                                            << " SAS: " << varNameVec.at(i) 
                                            << endl;
            if(vName == varNameVec.at(i)){
	      value = ANtpDefVal::kDouble;
	      valueI = ANtpDefVal::kInt;
	      if(!NeedsSpecialAttention(vNameOrig, nueR, value))

	      if (!strcmp(membertype->GetTypeName(),"float") ||
		  !strcmp(membertype->GetTypeName(),"Float_t") || 
		  !strcmp(membertype->GetTypeName(),"double")  ||
		  !strcmp(membertype->GetTypeName(),"Double_t")){
		value=atof(membertype->AsString(pointer));
		valueI=1;
	      }
	      else if(!strcmp(membertype->GetTypeName(),"int") ||
		      !strcmp(membertype->GetTypeName(),"Int_t")){
		value=atoi(membertype->AsString(pointer));
		valueI=atoi(membertype->AsString(pointer));
	      }
	      
	      else MSG("MdaDiscrimAna",Msg::kWarning)<<"Found variable "
						   << "NueRec: " << vName
						   << " of unknown type "
						   << membertype->GetTypeName()
						   << endl ;
		
	      MSG("MdaDiscrimAna",Msg::kDebug)<<"Found variable "
					     <<vName
					     <<" with value "
					     <<value
					     <<endl;
		
                if(!ANtpDefVal::IsDefault(value) && 
		   !ANtpDefVal::IsDefault(valueI)){
                    
                    valueTemp(i)=value;

                }else {
                    valueTemp(i)=(meanVec.at(0))(i);
                }
                MSG("MdaDiscrimAna",Msg::kDebug)<<"Found variable "
                                                <<vName<<" with value "
                                                <<value
                                                <<endl;                
                count++;
                i = varNameVec.size();
            }
        }
        if(count == varNameVec.size()) break;
    }
    
    varPVector.push_back(valueTemp);

    return;     
    
}

bool MdaDiscrimAna::isFilled;

void MdaDiscrimAna::MdaClassify()
{

    //Use the calculated class probabilities and user supplied threshold cut
    //to perform posterior event classification
    
    //Determine the maximum probability and its index

    if(probClass.size() < static_cast<UInt_t>(nClass)) return;

    IterDeqDouble_t maxIter =
        max_element(probClass.begin(),probClass.end());

    Int_t maxPos= distance(probClass.begin(),maxIter);
    MSG("MdaDiscrimAna",Msg::kDebug)<< "Prob Max found for entry " 
                                    << maxPos << ": " 
                                    << probClass.at(maxPos) 
                                    << " with class " 
                                    << (typeClassVec.at(maxPos)).at(0)
                                    << endl;
    
    //Classify into class ()
    Double_t probMax=probClass.at(maxPos);
    string classMax=(typeClassVec.at(maxPos)).at(0);
    
    if (probMax < threshCut) fMdaDiscrim.fMdaClass=-1;
    else if (classMax=="nue") fMdaDiscrim.fMdaClass=ClassType::nue;
    else if (classMax=="ncu") fMdaDiscrim.fMdaClass=ClassType::NC;
    else if (classMax=="num") fMdaDiscrim.fMdaClass=ClassType::numu;
    else if (classMax=="nut") fMdaDiscrim.fMdaClass=ClassType::nutau;
    else {
        MSG("MdaDiscrimAna",Msg::kError)
             <<"Unknown class in SAS calibration set" << endl;
        fMdaDiscrim.fMdaClass=ANtpDefVal::kInt;
    }
    return;

}


Bool_t MdaDiscrimAna::NeedsSpecialAttention(TString name
                                            , NueRecord *nr
                                            , Double_t &value)
{
   
    //All the fHeaders and four of the MST vars require special effort
    if(name == "fHeader.fSnarl") {
        value = nr->GetHeader().GetSnarl();
    }if(name == "fHeader.fRun") {
        value = nr->GetHeader().GetRun();
    }if(name == "fHeader.fSubRun") {
        value = nr->GetHeader().GetSubRun();
    }if(name == "fHeader.fEvtNo") {
         value = nr->GetHeader().GetEventNo();
    }if(name == "fHeader.fEvents") {
        value = nr->GetHeader().GetEvents();
    }if(name == "fHeader.fTrackLength") {
        value = nr->GetHeader().GetTrackLength();
    }

    if(name == "mstvars.eallw1") {
        if(nr->mstvars.enn1 > 0) value = 0.0;
        for(int i=0;i<nr->mstvars.enn1;i++){
            value += nr->mstvars.eallw1[i];
        }
    }
    if(name == "mstvars.oallw1") {
        if(nr->mstvars.onn1 > 0) value = 0.0;
        for(int i=0;i<nr->mstvars.onn1;i++){
            value += nr->mstvars.oallw1[i];
        }
    } 
    if(name == "mstvars.eallm1") {
        if(nr->mstvars.enn1 > 0) value = 0.0;
        for(int i=0;i<nr->mstvars.enn1;i++){
            value += nr->mstvars.eallm1[i];
        }
    }
    if(name == "mstvars.oallm1") {
        if(nr->mstvars.onn1 > 0) value = 0.0;
        for(int i=0;i<nr->mstvars.onn1;i++){
            value += nr->mstvars.oallm1[i];
        }
    }    
    
    if(value > -9999) return true;
    return false;
}


void MdaDiscrimAna::PIDCalc()
{

    DeqDouble_t discrimValue;
    
    //Calculate discriminant function Q for each class
    for (Int_t n=0; n < nClass; n++){
        
        //Create an auxiliary copy of the Pvector
        TVectorD varPClone=(varPVector.at(0));

        //Note that varPClone is modified by the "*=" operator
        Double_t quadTerm=((varPClone)*=(quadCoeff.at(n)))*(varPVector.at(0));

        Double_t linTerm=(linCoeff.at(n))*(varPVector.at(0));

        Double_t constTerm=(constCoeff.at(n))(0);
        
        discrimValue.push_back(quadTerm+linTerm+constTerm);
        
        MSG("MdaDiscrimAna",Msg::kDebug)<< "Q("
                                        << (typeClassVec.at(n)).at(0) 
                                        <<") = " 
                                        << discrimValue.at(n) <<endl;
    }
     
    TVectorD tempDiff(nClass-1);

    //Calculate discriminator difference vectors
    //P(a) = e^a/(e^(a)+e^(b)+e^(c)+...) = 1/(1+e^(b-a)+e^(c-a)+...)
    //Avoids e^X where X is big => machine-dep floating exception.
  
    for (Int_t n=0; n < nClass; n++){

        Int_t diffCount=0;
        
        for (Int_t m=0; m < nClass; m++){
            if(m!=n){
                tempDiff(diffCount)=(discrimValue.at(m)-discrimValue.at(n));
                diffCount++;
            } else continue;   
        }

        if (TMath::Abs(tempDiff.Min())>600. || TMath::Abs(tempDiff.Max())>600){
            MSG("MdaDiscrimAna",Msg::kDebug)<<"Bypass floating exception "
                                            << endl;
            continue;
        }else discrimDiff.push_back(tempDiff); 

        tempDiff.Zero();
    }
  
    if(discrimDiff.size()!=static_cast<UInt_t>(nClass)) return;

    //Calculate probability value for each class
    for (Int_t n=0; n < nClass; n++){
        
        Double_t denomSum=0.0;
        
        for (Int_t m=0; m < (nClass-1); m++){            
            denomSum+=TMath::Exp((discrimDiff.at(n))(m));
        }

        probClass.push_back(1.0/(1.0+denomSum));

        MSG("MdaDiscrimAna",Msg::kDebug) << "P("
                                         << (typeClassVec.at(n)).at(0) 
                                         <<") = " 
                                         << probClass.at(n) <<endl;
        
        if((typeClassVec.at(n)).at(0)=="nue") 
            fMdaDiscrim.fMdaPIDnue=probClass.at(n);        
        if((typeClassVec.at(n)).at(0)=="ncu") 
            fMdaDiscrim.fMdaPIDnc=probClass.at(n);        
        if((typeClassVec.at(n)).at(0)=="num") 
            fMdaDiscrim.fMdaPIDnumu=probClass.at(n);        
        if((typeClassVec.at(n)).at(0)=="nut") 
            fMdaDiscrim.fMdaPIDnutau=probClass.at(n);        
    }
    
    return;
}


void MdaDiscrimAna::Reset()
{
 

}
