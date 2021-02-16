#ifndef DTREADER_H
#define DTREADER_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

struct DTNode{
   int id;
   float score;
   int var;
   float cut;
   int left;
   int right;
};

class DecisionTreeReader{

   public: 
     DecisionTreeReader();
     ~DecisionTreeReader();

     bool CreateFromFile(string file);
     double Evaluate(double *input, int size);
     double Evaluate(vector<double> &input);
     double EvaluateTree(int index, vector<double> &input);

   private:
     bool ParseLine(ifstream &ins, DTNode* set);

     DTNode** Tree;
     int nTree;
     vector<int> nNodes;

};

#endif
