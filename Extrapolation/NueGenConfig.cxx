#include "NueGenConfig.h"
#include "NueAna/NueAnaTools/NueConvention.h"
                                                                               
#include <iostream>            
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <cstdlib>  // atoi, atof

using namespace std;

NueGenConfig::NueGenConfig(){
  Reset();
}

NueGenConfig::NueGenConfig(std::string input)
{
  Reset();
  ReadInput(input);
}

void NueGenConfig::Reset()
{ 
   fDataMethod = NueGenConfig::kSetList;

   for(int i = 0; i < 6; i++){
     fErrors[i] = 0.0;
     found[i] = true;
   }

   fReadUntil = -1;
   fOffSet = 0;
}

void NueGenConfig::ReadInput(std::string input)
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
//           std::cout<<"Setting method "<<fDataMethod<<std::endl;
         }
         if(line.find("FILE")!=std::string::npos){
           if(fDataMethod != NueGenConfig::kFileList){
              std::cout<<"Error in config file, invalid data input"<<std::endl;
              continue;
           }
           temp = line.substr(line.find("{")+1,
                     line.find_first_of("}")-line.find("{")-1);
           files.push_back(temp);
           continue;
         }
//        std::cout<<"On line: "<<line<<std::endl;         
       }
    }
    if(line.find("BEGINERROR") != std::string::npos){
       while (line.substr(0,8) != "ENDERROR"){
         stream.getline(cline,200);
         line=cline;
         if (line.find("NuMu")!=std::string::npos){
           SetErrorVal(ClassType::numu, line); continue;
         }
         if (line.find("NC")!=std::string::npos){
           SetErrorVal(ClassType::NC, line); continue;
         }
         if (line.find("BNue")!=std::string::npos){
           SetErrorVal(ClassType::bnue, line); continue;
         }
         if (line.find("NuTau")!=std::string::npos) {
           SetErrorVal(ClassType::nutau, line); continue;
         }
         if (line.find("SigNue")!=std::string::npos){
           SetErrorVal(ClassType::nue, line);  continue;
         }
	 if (line.find("Total")!=std::string::npos){
           SetErrorVal(5, line);  continue;
         }
//         std::cout<<" err stuff "<<line<<std::endl;

       }
       std::cout<<"Err section"<<std::endl;
    }  //End of Parsing the input data
  
    if(line.find("BEGINSET") != std::string::npos){
       GCCompleteSet set;
       std::cout<<"beg set"<<std::endl;
       while (line.substr(0,6) != "ENDSET"){
         stream.getline(cline,200);
         line=cline;
         if (line.find("/",0)==0||line.size()==0) continue;
         if (line.find("NuMu")!=std::string::npos){
           SetDataRange(ClassType::numu, line, set); continue;
         }
         if (line.find("NC")!=std::string::npos){
           SetDataRange(ClassType::NC, line, set); continue;
         }
         if (line.find("BNue")!=std::string::npos){
           SetDataRange(ClassType::bnue, line, set); continue;
         }
         if (line.find("NuTau")!=std::string::npos) {
           SetDataRange(ClassType::nutau, line, set); continue;
         }
         if (line.find("SigNue")!=std::string::npos){
           SetDataRange(ClassType::nue, line,set);  continue;
         }
//         std::cout<<" set stuff "<<line<<std::endl;

       } //End of while loop
       EventGroups.push_back(set);
       std::cout<<" set stuff "<<line<<std::endl;
    } //End of Par data
       
  }  
  std::cout <<"Finished Reading Nue Data "<<std::endl;
 
  CheckConfig();                                                                     
} 

void NueGenConfig::SetDataRange(int nuClass, std::string line, GCCompleteSet &set)
{
  if(fDataMethod != NueGenConfig::kSetList){
//    std::cout<<"Error in config file, invalid data input"<<std::endl;
    return;
  }

  std::string temp;
  while (line.find(" ",0)<line.size())
              line.replace(line.find(" ",0),1,"");
                                                                                     
  temp = line.substr(line.find("{")+1,
              line.find_first_of(",")-line.find("{")-1);
  set.EventClass[nuClass].start = atof(temp.c_str());

  temp = line.substr(line.find_first_of(",")+1,
                     line.find_last_of(",")-line.find_first_of(",")-1);
  set.EventClass[nuClass].end = atof(temp.c_str());

  temp = line.substr(line.find_last_of(",")+1,
                     line.find("}")-line.find_last_of(",")-1);
  set.EventClass[nuClass].step = atof(temp.c_str());
                                                                                     
  set.EventClass[nuClass].isfixed = false;
  if(line.find("fix")!=std::string::npos) set.EventClass[nuClass].isfixed = true;
                                                                               
}

void NueGenConfig::SetErrorVal(int pos, std::string line)
{
   std::string temp = line.substr(line.find("{")+1,
                      line.find("}")-line.find("{")-1);
   fErrors[pos] = atof(temp.c_str());
   found[pos] = true;

}

bool NueGenConfig::CheckConfig(){

  if(fDataMethod < 0 || fDataMethod > 2) return false;
  if(fDataMethod == 1 || fDataMethod == 0){
    for(int i = 0; i < 6; i++) { if(!found[i]) return false; }
  }
    
  return true;
}
 
bool NueGenConfig::LoadNextNumberSet(double *num)
{
  static std::ifstream ins;
  static int count = 0;
  //  static place pos;
//  std::cout<<"Entering number loader"<<std::endl;

  if(fDataMethod == NueGenConfig::kFileList)
  {
    std::string dum;
    if(count >= fReadUntil && fReadUntil > 0) return false;

    if(!ins.is_open()){
        ins.open(files[0].c_str());
        int skipped = 0;

        if(fOffSet > 0){
           do{
             ins>>fOscPar[0]>>fOscPar[1]>>fOscPar[2]>>fOscPar[3];
             ins>>num[0]>>num[1]>>num[4]>>num[3]>>num[2]>>dum;
             skipped++;
           }while(skipped < fOffSet && !ins.eof());

        }
      std::cout<<"Reading: "<<fReadUntil<<" skipped: "<<fOffSet<<std::endl;
    }
    ins>>fOscPar[0]>>fOscPar[1]>>fOscPar[2]>>fOscPar[3];
    ins>>num[0]>>num[1]>>num[4]>>num[3]>>num[2]>>dum;
    count++;
    if(ins.eof()){  ins.close();  return false; }
  }
  if(fDataMethod == NueGenConfig::kSetList)
  {
     
  }




        

  return true;
}  
