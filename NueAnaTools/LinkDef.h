#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all namespaces;

#pragma link C++ namespace NueConvention;
#pragma link C++ namespace ClassType;
#pragma link C++ namespace SntpHelpers;

#pragma link C++ namespace     Selection;
#pragma link C++ enum          Selection::ESelection;
#pragma link C++ nestedtypedef Selection::Selection_t;

#pragma link C++ class DCHit+;
#pragma link C++ class vector<DCHit>+;
#pragma link C++ class vector<DCHit>::iterator+;
#pragma link C++ class DCVertex<DCHit>+;
#pragma link C++ class DCEdge<DCHit>+;
#pragma link C++ class DCGraph<DCHit>+;
#pragma link C++ class DecisionTreeReader+;

//#pragma link C++ class OscWeight+;
//#pragma link C++ namespace OscPar+;
//#pragma link C++ class PMNS+;
//#pragma link C++ class OscCalc+;


//#pragma link C++ class Trimmer+;
//#pragma link C++ class NueAnalysisCuts+;

#endif
