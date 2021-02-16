#include "NueSenseConfig.h"
#include "NueAna/NueAnaTools/NueConvention.h"
                                                                               
#include <iostream>            
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <cstdlib>  // atoi, atof

NSCErrorParam::NSCErrorParam(){
  bg_systematic = 0.0;
  sig_systematic = 0.0;
  for(int i = 0; i < 5; i++) scale[i] = 0.0;
}
                                                                                          
NueSenseConfig::NueSenseConfig(){
  Reset();
}

NueSenseConfig::NueSenseConfig(std::string input)
{
  Reset();
  ReadInput(input);
}

void NueSenseConfig::Reset()
{ 
   fDataMethod = 0;
   fPOT = 9.25;
   fOscUe3Square = -1.0;
   fOscDeltaMS23 = -1.0;
   fMustDeOsc = false;
   DeltaMS12 = 8.7;
   SinS2Th12 = 0.816;
   SinS2Th23 = 1.0;

   for(int i = 0; i < 5; i++){
     number[i] = 0.0;
     found[i] = true;
   }

   DeltaMS23.isfixed = true;
   DeltaMS23.start = DeltaMS23.end = 2.7;

   Delta.isfixed = false;
   Delta.start = 0;  Delta.end = 2.0;

   SinS2Th13.isfixed = false;
   SinS2Th13.start = 0.0;  SinS2Th13.end = 0.4;

   density = 2.65;
}

void NueSenseConfig::ReadInput(std::string input)
{
  std::ifstream stream;
  stream.open(input.c_str());
  assert(stream.is_open() && "Can't open input file");
  std::cout <<"Reading Nue Data from "<<input<<std::endl;

  std::string line, temp;
  char cline[200];
                                                                                
  while (stream.getline(cline,200)){
    line=cline;
    if (line.find("/",0)==0||line.size()==0) continue;  //ignore comments
   
 
    if(line.find("BEGINDATA") != std::string::npos){
       while (line.substr(0,7) != "ENDDATA"){
         stream.getline(cline,200);
         line=cline;
         if (line.find("/",0)==0||line.size()==0) continue;

         if (line.find("METHOD") != std::string::npos) {
           fDataMethod = atoi((line.substr(7,line.size()-7)).c_str());
         }
         if (line.find("POT") != std::string::npos) {
           fPOT = atof((line.substr(3,line.size()-3)).c_str());
           continue;
         }
         if (line.find("Numu")!=std::string::npos){
            SetDataInfo(ClassType::numu, line);  continue;
         }
         if (line.find("NC")!=std::string::npos){
            SetDataInfo(ClassType::NC, line); continue;
         }
         if (line.find("BNue")!=std::string::npos){
            SetDataInfo(ClassType::bnue, line); continue;
         }
         if (line.find("Nutau")!=std::string::npos) {
            SetDataInfo(ClassType::nutau, line); continue;
         }
         if (line.find("SigNue")!=std::string::npos){
            SetDataInfo(ClassType::nue, line);  continue;
         }
         if (line.find("NOOSC")!=std::string::npos){
           fOscUe3Square = fOscDeltaMS23 = -1;
           fMustDeOsc = false;  continue;
         }
         if (line.find("OSCPAR")!=std::string::npos){
           temp = line.substr(line.find("{")+1,
                     line.find_first_of(",")-line.find("{")-1);
           fOscUe3Square = atof(temp.c_str());
           temp = line.substr(line.find_last_of(",")+1,
                     line.find("}")-line.find_last_of(",")-1);
           fOscDeltaMS23 = atof(temp.c_str());
           fMustDeOsc = true;  continue;
         }
         if(line.find("FILE")!=std::string::npos){
           if(fDataMethod != 3 && fDataMethod != 4){
             std::cout<<"Error in config file, invalid data input"<<std::endl;
             continue;
           }
           temp = line.substr(line.find("{")+1,
                     line.find_first_of("}")-line.find("{")-1);
           datafiles.push_back(temp);
           continue;
         }
      }
    }  //End of Parsing the input data

    if(line.find("BEGINPAR") != std::string::npos){
       while (line.substr(0,6) != "ENDPAR"){
         stream.getline(cline,200);
         line=cline;
         if (line.find("/",0)==0||line.size()==0) continue;
         
         if (line.find("DELTAM2_23")!=std::string::npos) { 
            SetParamInfo(DeltaMS23,  line);  continue;
         }
         if (line.find("SIN2(2TH13)")!=std::string::npos){
            SetParamInfo(SinS2Th13,  line);  continue;
         }
         if (line.find("DELTACP")!=std::string::npos) {
            SetParamInfo(Delta,  line);      continue;
         }
         if (line.find("DELTAM2_12")!=std::string::npos){
            DeltaMS12 = atof((line.substr(10,line.size()-10)).c_str());
            continue;
         }
         if (line.find("SIN2(2TH23)")!=std::string::npos) {
            SinS2Th23 = atof((line.substr(11,line.size()-11)).c_str());
            continue;
         }
         if (line.find("SIN2(2TH12)")!=std::string::npos) {
            SinS2Th12 = atof((line.substr(11,line.size()-11)).c_str());
            continue;
         }
         if (line.find("DENSITY")!=std::string::npos) {
            density = atof((line.substr(8,line.size()-7)).c_str());
            continue;
         }

       } //End of while loop
    } //End of Par data
    
    if(line.find("BEGINSET") != std::string::npos){
       NSCErrorParam set;
       while (line.substr(0,6) != "ENDSET"){
         stream.getline(cline,200);
         line=cline;
         if (line.find("/",0)==0||line.size()==0) continue;  //ignore comments

         if (line.find("BG_SYSTEMATIC") != std::string::npos) {
           set.bg_systematic = atof((line.substr(13,line.size()-13)).c_str());
           continue;
         }
         if (line.find("SIG_SYSTEMATIC") != std::string::npos) {
           set.sig_systematic = atof((line.substr(14,line.size()-14)).c_str());
           continue;
         }
         if (line.find("NCSCALE") != std::string::npos) {
           set.scale[ClassType::NC] = atof((line.substr(7,line.size()-7)).c_str());
           continue;
         }
         if (line.find("NUMUSCALE") != std::string::npos) {
           set.scale[ClassType::numu] = atof((line.substr(9,line.size()-9)).c_str());
           continue;
         }
         if (line.find("BNUESCALE") != std::string::npos) {
           set.scale[ClassType::bnue] = atof((line.substr(9,line.size()-9)).c_str());
           continue;
         }
         if (line.find("NUTAUSCALE") != std::string::npos) {
           set.scale[ClassType::nutau] = atof((line.substr(10,line.size()-10)).c_str());
           continue;
         }
         if (line.find("SIGNUESCALE") != std::string::npos) {
           set.scale[ClassType::nue] = atof((line.substr(12,line.size()-12)).c_str());
           continue;
         }
       }
       //create a name for the set?
       // put the details in a separate tree in each subdir the way chris would do things?
       errorsets.push_back(set);
    }                           
  }  
  std::cout <<"Finished Reading Nue Data "<<std::endl;

  CheckConfig();                                                                     
} 

void NueSenseConfig::SetDataInfo(int nuClass, std::string line)
{
  if(fDataMethod != 1 && fDataMethod != 2){
    std::cout<<"Error in config file, invalid data input"<<std::endl;
    return;
  }

  std::string temp;
  found[nuClass] = true;

  while (line.find(" ",0)<line.size())
           line.replace(line.find(" ",0),1,"");
                                                                               
  if(line.find("Number") != std::string::npos){
     temp = line.substr(line.find("{")+1,line.find("}")-line.find("{")-1);
     number[nuClass] = atof(temp.c_str());
  } else if(line.find("Hist") != std::string::npos){
     if(fDataMethod != 2) {
       std::cout<<"Error in config file, invalid data input"<<std::endl;
       return;
     }
     signuehfile = line.substr(line.find("{")+1,
                     line.find_first_of(",")-line.find("{")-1);
     signuehname = line.substr(line.find_last_of(",")+1,
                     line.find("}")-line.find_last_of(",")-1);
   }
   else{
      found[nuClass] = false;
  }
}

void NueSenseConfig::SetParamInfo(NSCDataParam &par, std::string line)
{
  std::string temp;
  while (line.find(" ",0)<line.size())
              line.replace(line.find(" ",0),1,"");
                                                                                   
  temp = line.substr(line.find("{")+1,
              line.find_first_of(",")-line.find("{")-1);
  par.start = atof(temp.c_str());
  temp = line.substr(line.find_last_of(",")+1,
                     line.find("}")-line.find_last_of(",")-1);
  par.end = atof(temp.c_str());
                                                                                   
  par.isfixed = false;
  if(line.find("fix")!=std::string::npos) par.isfixed = true;
}

bool NueSenseConfig::CheckConfig(){

  if(fDataMethod < 1 || fDataMethod > 4) return false;
  if(fDataMethod == 1 || fDataMethod == 2){
    for(int i = 0; i < 5; i++) { if(!found[i]) return false; }
  }
  if(fDataMethod == 1 || fDataMethod == 2){
    if(fPOT < 0) return false;
    if(number[0] < 0) return false;
    if(number[1] < 0) return false;
    if(number[3] < 0) return false;
    if(number[4] < 0) return false;
    if(fMustDeOsc && (fOscUe3Square < 0 || fOscDeltaMS23 < 0)) 
       return false;
  }
  if(fDataMethod == 1){
    if(!fMustDeOsc) return false;    //Can't just have a number without some oscillation
    if(number[2] < 0) return false;
    DeltaMS23.start = 2.7;    DeltaMS23.isfixed = true;
    Delta.start = 0.0;        Delta.isfixed = true;
  }
    
  if(fDataMethod == 2 && (signuehfile.size() == 0 || signuehname.size() == 0))
      return false;
 
  if(fDataMethod == 3 && datafiles.size() == 0) return false;
  if(fDataMethod == 4 && datafiles.size() == 0) return false;

  if(DeltaMS12 < 0) return false;
  if(SinS2Th12 < 0) return false;
  if(SinS2Th23 < 0) return false;
  if(!DeltaMS23.isfixed && DeltaMS23.start > DeltaMS23.end) return false;
  if(!Delta.isfixed && Delta.start > Delta.end) return false;
  if(!SinS2Th13.isfixed && SinS2Th13.start > SinS2Th13.end) return false;

  if(DeltaMS23.isfixed) { DeltaMS23.end = DeltaMS23.start; }
  if(Delta.isfixed)     { Delta.end = Delta.start; }
  if(SinS2Th13.isfixed) { SinS2Th13.end = SinS2Th13.start; }

  if(errorsets.size() == 0){ 
    NSCErrorParam set;
    errorsets.push_back(set);
  }
  
  return true;
}
 
  
