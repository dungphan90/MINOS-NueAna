//////////////////////////////////////////////////////////////////////////////
// TrimParah.h
//
// This helps with setting up a large number of trimmming runs at once.
//   This is meant for use with the TrimEngine.C 
//
// Created:  J Boehm -- March, 2007
//
// $Author: boehm $
//
// $Revision: 1.1 $
//
// $Name:  $
//
// $Id: TrimParams.h,v 1.1 2007/04/03 17:13:35 boehm Exp $
//
// $Log: TrimParams.h,v $
// Revision 1.1  2007/04/03 17:13:35  boehm
// Brand new trimmer which runs outside of job control.  This allows it to run about 5x faster than the old trimmer.  I've had to employ some tricks to make everything work together reasonably - so if you want to run this yourself you may want to contact me first - Josh
//
//
//
/////////////////////////////////////////////////////////////////////


#include <vector>
#include <string>
#include <iostream>

using namespace std;

void AddNearMCSet(string name);
void AddFirstYearNearData();


std::string base = "/afs/fnal.gov/files/data/minos/d139/Chimaera/Cedar/";
int count = 0;

std::vector<std::string> fileSets[100];

void AddNearMCSet(std::string name)
{
   fileSets[count].push_back(base+name+"0*.root");
   fileSets[count].push_back(base+name+"1*.root");
   fileSets[count].push_back(base+name+"2*.root");
   fileSets[count].push_back(base+name+"3*.root");
   fileSets[count].push_back(base+name+"4*.root");
   count++;

   fileSets[count].push_back(base+name+"5*.root");
   fileSets[count].push_back(base+name+"6*.root");
   fileSets[count].push_back(base+name+"7*.root");
   fileSets[count].push_back(base+name+"8*.root");
   fileSets[count].push_back(base+name+"9*.root");
   count++;
}

void AddFirstYearNearData(){
  fileSets[count].push_back(base+"NearData/2005-05/*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2005-06/*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2005-07/*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2005-08/*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2005-09/AnaNue-N000084*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2005-09/AnaNue-N000085*.root");
  fileSets[count].push_back(base+"NearData/2005-09/AnaNue-N000086*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2005-10/*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2005-11/AnaNue-N000090*.root");
  fileSets[count].push_back(base+"NearData/2005-11/AnaNue-N000091*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2005-11/AnaNue-N000092*.root");
  fileSets[count].push_back(base+"NearData/2005-11/AnaNue-N000093*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2005-12/*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N000095*.root");
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N0000960*.root");
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N0000961*.root");
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N0000962*.root");
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N0000963*.root");
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N0000964*.root");
  count++;
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N0000965*.root");
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N0000966*.root");
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N0000967*.root");
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N0000968*.root");
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N0000969*.root");
  fileSets[count].push_back(base+"NearData/2006-01/AnaNue-N000097*.root");
  count++; 
  fileSets[count].push_back(base+"NearData/2006-02/*.root");
  count++;
}

void AddFarMCSet(int i)
{
   string t = "0";
   if(i == 0) t = "0";
   if(i == 3) t = "3";
   if(i == 4) t = "4";

   fileSets[count].push_back(base+"Daikon/FarL010185/AnaNue-f21"+t+"110*.root");
   count++;
   fileSets[count].push_back(base+"Daikon/FarL010185/AnaNue-f21"+t+"111*.root");
   count++;
   fileSets[count].push_back(base+"Daikon/FarL010185/AnaNue-f21"+t+"112*.root");
   count++;
   fileSets[count].push_back(base+"Daikon/FarL010185/AnaNue-f21"+t+"113*.root");
   count++;
   fileSets[count].push_back(base+"Daikon/FarL010185/AnaNue-f21"+t+"114*.root");
   count++;
   fileSets[count].push_back(base+"Daikon/FarL010185/AnaNue-f21"+t+"115*.root");
   fileSets[count].push_back(base+"Daikon/FarL010185/AnaNue-f21"+t+"116*.root");
   count++;
}

