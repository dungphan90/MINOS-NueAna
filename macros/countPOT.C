/////////////////////////////////////////////////////////////////////////////////////////
//
// Task:  
//       Return the total POT of the AnaNue files
//
// Inputs:
//        This macro can take 3 types of inputs:
//         1)A single AnaNue file that must end with ".root"
//           for example: "AnaNue.root"
//         2)A collection of AnaNue file name using a "*"
//           for example: "MyFiles/AnaNue*"
//         3)The name of an ASCII file containing a list of AnaNue files ('\n' delimited)
//           for example: "AnaNue.dat"
//
// Outputs:
//        The number of POT
//
// How to run it:
//      #Setup the nue analysis software
//        <YourPrompt> setup_nue            
//      #Compile the LOON macro
//        <YourPrompt> CompileLoonMacro countPOT.C 
//      #Run the LOON macro
//        <YourPrompt> loon -bq 'countPOT.C+("MyFiles/AnaNue*")'
//      
//
// Author: Gregory Pawloski 
// Created: August, 2010
////////////////////////////////////////////////////////////////////////


//------------------------------------------------------------------------------
//include minossoft related headers
#include "NueAna/NueUtilities.h"
//------------------------------------------------------------------------------

using namespace std;

void countPOT(string inputFileList) 
{
  NueUtilities::AnaNueProcessor anp(inputFileList);
  cout << "\n\n Processing " << anp.GetTotalPOT()/100000000.0 << " x 10^{20} POT\n" << endl;
}








