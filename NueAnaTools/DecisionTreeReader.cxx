#include <cassert>
#include <fstream>
#include <iostream>
#include "DecisionTreeReader.h"

using namespace std;

DecisionTreeReader::DecisionTreeReader()
{
   Tree = 0;
   nTree = 0;
   nNodes.clear();
}

DecisionTreeReader::~DecisionTreeReader()
{
   for(int i = 0; i < nTree; i++)
   {
        delete [] Tree[i];
   }
   delete [] Tree;
   nNodes.clear();

}

bool DecisionTreeReader::CreateFromFile(string file)
{
   ifstream ins(file.c_str());
   if(!ins.is_open()){
      cout<<"Unable to open file: "<<file<<endl;
      return false;
   }

   string dum1, dum2;
   int hold = 0;
   
   getline(ins, dum1);  //read and ignore first line
//   cout<<dum1<<endl;
   ins>>dum2>>hold;     //get Number of Trees
//   cout<<dum2<<hold<<endl;
   if(dum2 != "Classifiers:" || hold < 0){
     cout<<"Problem at: "<<dum2<<"  "<<hold<<endl;
     return false;
   }
   nTree = hold;
 
   Tree = new DTNode*[nTree]; 
   for(int i = 0; i < nTree; i++){  
       getline(ins, dum1);
       getline(ins, dum1);
       getline(ins, dum1);

       ins>>dum2>>hold>>dum2;  //get number of nodes
       if(dum2 != "nodes." || hold < 0){
          cout<<"Problem at (53): "<<dum2<<"  "<<hold<<endl;
          return false;
       }
       nNodes.push_back(hold);

       Tree[i] = new DTNode[hold];
	for(int j = 0; j < hold; j++){
            if(!ParseLine(ins, &Tree[i][j])) 
               return false;
            if(Tree[i][j].id != j) cout<<"Sheer problem"<<endl;
	}
   }

   return true;
//   then a few more lines giving the variables;

}

bool DecisionTreeReader::ParseLine(ifstream &ins, DTNode* set)
{
   string dum1, dum2, dum3, dum4;
   int id, left, right, var;
   float cut, score;

   ins>>dum1>>id>>dum2>>score>>dum3>>var>>dum3>>cut>>dum4>>left>>right;
   if(dum1 != "Id:" || dum2 != "Score:" || dum4 != "Daughters:")
      return false;
 
   set->score = score;
   set->id = id;
   set->var = var;
   set->cut = cut;
   set->left = left;
   set->right = right;

   return true;
}

double DecisionTreeReader::EvaluateTree(int index, vector<double> &input)
{
  DTNode* pos = &Tree[index][0];
  while( pos->var >= 0 ) {
    assert( pos->var < int(input.size()) );
    if( input[pos->var] < pos->cut )
      pos = &Tree[index][pos->left];
    else
      pos = &Tree[index][pos->right];
  }
  return pos->score;
}

double DecisionTreeReader::Evaluate(double *input, int size)
{
  vector<double> out;
  for(int i = 0; i < size; i++){
     out.push_back(input[i]);
  }
//  cout<<"here"<<endl;
  return Evaluate(out);
}

double DecisionTreeReader::Evaluate(vector<double> &input)
{
//  cout<<"Try again:  "<<nTree<<endl;
  double res = 0;

  for(int i = 0; i < nTree; i++){
//     cout<<res<<endl;
     res += EvaluateTree(i, input);     
  }
  res /= nTree;

  return res;
}
