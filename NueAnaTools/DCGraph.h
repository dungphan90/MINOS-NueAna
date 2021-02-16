#ifndef DCGRAPH_H
#define DCGRAPH_H
#include <iterator>
#include <vector>
#include <map>
#include <set>
#include "TCanvas.h"
#include "TObject.h"
#include "TH1F.h"
#include "NueAna/NueAnaTools/DCEdge.h"

using std::vector;
using std::cout;
using std::endl;
using std::multiset;

template<class T> class DCVertex;

template<class T> class DCGraph : public TObject
{
public:
   DCGraph(){first=0;nverts=0;nedges=0;
             summetric=0.;sumxmetric=0.;sumymetric=0.;maxmetric=0.;}

   virtual ~DCGraph() { if(first!=0) {delete first; first=0;}}
   
   DCVertex<T> *Find(T *);
   bool AddVertex(T*);
   bool AddEdge(T* ,T*);
   int GetNVerts() const {return nverts;}
   int GetNEdges() const {return nedges;}
   float GetSummedMetric() const {return summetric;}
   float GetXSummedMetric() const {return sumxmetric;}
   float GetYSummedMetric() const {return sumymetric;}
   float GetMaxMetric() const {return maxmetric;}
   float GetSummedZ() const;
   DCVertex<T> *GetFirst() {return first;}

   vector<float> GetAllWeights();
   vector<float> GetAllZ();
   void NewFirst(DCVertex<T> *nf) { first=nf; }
   virtual void Print(Option_t *t="") const;
   virtual void Draw(Option_t *t="");

   vector<DCGraph<T>* > FindAllMinSpan(float maxmetric=0., float minz=0., 
				       float maxlowzmetric=0.);
   DCGraph<T>* FindMinSpan(float maxmetric=0.,float minz=0., 
			   float maxlowzmetric=0.);
   //   bool AddClosestNode(DCGraph<T> *out, float maxmetric=0., float minz=0.,
   //		       float maxlowzmetric=0.);

   bool AddClosestNode(vector<T> &out, float maxmetric=0., float minz=0.,
   		       float maxlowzmetric=0.);

   vector<T> Complement(DCGraph<T>*);

   DCGraph<T> *GraphMaxHits(int nhits);

private:
   DCVertex<T> *first;  // a pointer to the first vertex of the graph
   //   DCVertex<T> *eventvertex; //a pointer to the "event vertex"
                             //(currently biggest hit in most upstream plane

   int nverts;
   int nedges;
   float summetric;
   float sumxmetric;
   float sumymetric;
   float maxmetric;

   ClassDef(DCGraph<T>,1)
};
ClassImpT(DCGraph,T)

// the method returns the pointer to the vertex with data theData 
// if such vertex occurs in the graph, otherwise returns NULL
template<class T> DCVertex<T> *DCGraph<T>::Find(T* theData) {
  DCVertex<T> *vPtr = first;

  while (vPtr != 0) {
     // check if the data of the vertex is theData
     if (*vPtr->GetData() == *theData)
	return vPtr;
     // go to the next vertex
     vPtr = vPtr->GetNextVertex();
  }

  // there has been no match in the loop, so there is no such vertex
  return 0;
}

// the method checks if there is a vertex with data theData in the graph
// if yes, the vertex cannot be added, the method returns 'false'
// otherwise the vertex is added, the method returns 'true'
template<class T> bool DCGraph<T>::AddVertex(T* theData) { 
  // checking if the vertex is already in the graph
  // if it is, it cannot be added again 
  //   cout<<"In DCGraph::AddVertex"<<endl;
  if (Find(theData) != 0 )
    return false;

  // allocate memory for new vertex with data theData,
  // make it point to the previous first vertex
  DCVertex<T> *newVertex = new DCVertex<T>(theData, first);
  // make the new vertex the first one in the list of vertexes
  first = newVertex;
  nverts++;
  return true;
}

// check if the vertexes with data Begin and End are present in the graph. 
// if yes, create an edge from Begin to End, return 'true'
// otherwise return 'false'
template<class T> bool DCGraph<T>::AddEdge(T* Begin, T* End) {
  //   cout<<"In DCGraph::AddEdge"<<endl;
  // find the pointer to the end vertex
   DCVertex<T> *vEnd = Find(End);

  // if the vertex is not in the graph, cannot add the edge
  if (vEnd == 0) return false;

  // find the pointer to the start vertex
  DCVertex<T> *vBegin = Find(Begin);

  // if the vertex is not in the graph, cannot add the edge
  if (vBegin == 0) return false;

  // connect the start vertex to the end one
  summetric+=vBegin->ConnectTo(vEnd);
  if(vBegin->GetFirstEdge()->GetMetric()>maxmetric){
     maxmetric=vBegin->GetFirstEdge()->GetMetric();
  }
  sumxmetric+=vBegin->GetFirstEdge()->GetXMetric();
  sumymetric+=vBegin->GetFirstEdge()->GetYMetric();

  nedges++;
  return true;
}
  
// print all vertexes and edges of the graph
template<class T> void DCGraph<T>::Print(Option_t */*t*/) const {
   cout<<"Tree has "<<nverts<<" Vertices and "<<nedges<<" Edges"<<endl;
   cout<<"Summed metric is "<<summetric
       <<" and SummedZ is "<<GetSummedZ()<<endl;
   cout<<"SummedXmetric is "<<sumxmetric
       <<" and SummedYmetric is "<<sumymetric<<endl;

  if (first == 0) 
    cout << "Graph has no vertexes " << endl;
  else
    first->PrintGraph();
}

template<class T> void DCGraph<T>::Draw(Option_t */*t*/)
{
//   TCanvas *pad = (TCanvas *)gROOT->GetSelectedPad();
   TCanvas *pad = (TCanvas *)(gPad);
   if(pad==0){
//      cout<<"no pad selected, making new"<<endl;
      pad = new TCanvas();
   }
//   cout<<"Drawing to canvas "<<pad->GetName()<<endl;
   pad->Clear();
   
   float x0=0.;
   float y0=0.;
   float z0=0.;
   first->GetData()->GetData(x0,y0,z0);
   DCVertex<T> *v=first;
   while(v!=0){
      float x=0;
      float y=0;
      float z=0;
      v->GetData()->GetData(x,y,z);
      if(x<x0){
	x0=x;
	y0=y;
      }
      v=v->GetNextVertex();
   }


   //   eventvertex->GetData()->GetData(x0,y0,z0);
   pad->DrawFrame(x0-140,y0-164,x0+700,y0+164);
   if (first == 0) 
      cout << "Graph has no vertexes " << endl;
   else{
      first->DrawGraph();
   }
}

template<class T> DCGraph<T>* DCGraph<T>::FindMinSpan(float maxmetric, 
						      float minz,
						      float maxlowzmetric)
{
//   cout<<"In DCGraph::FindMinSpan"<<endl;
   DCGraph<T> *msg = new DCGraph<T>();
   
   DCVertex<T> *v=first;
   msg->AddVertex(v->GetData());
   vector<T> out = msg->Complement(this);
    while(msg->AddClosestNode(out, maxmetric, minz, maxlowzmetric)){};
   return msg;
}

template<class T> bool DCGraph<T>::AddClosestNode(vector<T> &out, 
						  float maxmetric, 
						  float minz, 
						  float maxlowzmetric)
{
//   cout<<"In DCGraph::AddClosestNode"<<endl;
   float minmetric=-1.;
   //   DCVertex<T> *closestin=0;
   //   DCVertex<T> *closestout=0;
   T *closestin=0;
   T *closestout=0;
   int closestoutindex=0;
   DCVertex<T> *vin=first;
   while(vin!=0){      
     //      DCVertex<T> *vout=out->GetFirst();
     for(unsigned int c=0;c<out.size();c++){
       T *vout = &out[c];
       if(Find(vout)!=0){
	 continue;
       }
       float m = ((*vin->GetData())-(*vout)).abs();
       float x,y,z;
       vout->GetData(x,y,z);
       if(minmetric<0){
	 if(z<minz&&m>maxlowzmetric&&maxlowzmetric!=0){	    
	 }
	 else if(m>maxmetric&&maxmetric!=0){
	 }
	 else{   
	   minmetric = m;
	   closestout = vout;
	   closestin = vin->GetData();
	   closestoutindex=c;
	 }
       }
       if(m<minmetric){
	 if(z<minz&&m>maxlowzmetric&&maxlowzmetric!=0){	    
	 }
	 else if(m>maxmetric&&maxmetric!=0){
	 }
	 else{
	   minmetric = m;
	   closestout = vout;
	   closestin = vin->GetData();
	   closestoutindex=c;
	 }
       }
     }
     if(vin->GetNextVertex()){
       vin=vin->GetNextVertex();
     }
     else{
       vin=0;
     }
   }
   
   if(closestin==0||closestout==0){
     //      cout<<"couldnt find a closest node returning false"<<endl;
     return false;
   }

   else{
     //      cout<<"Adding ";
     //      closestout->Print();
     //      cout<<"with metric "<<minmetric<<endl;
     AddVertex(closestout);
     AddEdge(closestin,closestout);
     //     vector<T>::iterator vi(&(out[closestoutindex]));
     out.erase(out.begin()+closestoutindex);
     return true;
   }

   return false;
}

template<class T> float DCGraph<T>::GetSummedZ() const
{
//   cout<<"In DCGraph::GetSummedZ"<<endl;
   float sumz=0;
   if(this==0){
      return sumz;
   }
   DCVertex<T> *vert = first;
   while(vert){
      float x=0;
      float y=0;
      float z=0;
      vert->GetData()->GetData(x,y,z);
      sumz+=z;
      vert=vert->GetNextVertex();
   }
   return sumz;
}

template<class T> vector<T> DCGraph<T>::Complement(DCGraph<T> *g)
{
//   cout<<"In DCGraph::Complement"<<endl;
   vector<T> vout;
   DCVertex<T> *v = g->GetFirst();
   while(v!=0){
      if(Find(v->GetData())==0){
	 vout.push_back(*(v->GetData()));
      }
      v=v->GetNextVertex();
   }
   return vout;
}

template<class T> vector<DCGraph<T> *> DCGraph<T>::FindAllMinSpan(float maxmetric, float minz, float maxlowzmetric)
{
//   cout<<"In DCGraph::FindAllMinSpan"<<endl;
   vector<DCGraph<T> *> trees;
   DCGraph<T> *msg = FindMinSpan(maxmetric,minz,maxlowzmetric);
   trees.push_back(msg);
   vector<T> out=msg->Complement(this);
   while(out.size()>0){
//      cout<<"DCGraph::FindAllMinSpan() found "<<out.size()
//	  <<" vertices out"<<endl;
      DCGraph<T> ng;
      for(unsigned int i=0;i<out.size();i++){
//	 cout<<"adding vertices to ng"<<endl;
	 ng.AddVertex(&out[i]);
      }
      DCGraph<T> *msg2 = ng.FindMinSpan(maxmetric,minz,maxlowzmetric);
//      cout<<"Printing new msg2"<<endl;
//      msg2->Print();
//      cout<<"***"<<endl;
      //make newout vector
      vector<T> newout=out;
      out.clear();
//      cout<<"size of out now is "<<out.size()<<endl;
      for(unsigned int i=0;i<newout.size();i++){
	 if(msg2->Find(&newout[i])==0){
//	    cout<<"found a vertex in out that wasn't in the new tree"<<endl;
//	    newout[i].Print();
	    out.push_back(newout[i]);
	 }
      }

      trees.push_back(msg2);
//      cout<<"Found a new minspan, starting loop again."<<endl;
//      cout<<"Size of trees now: "<<trees.size()<<endl;
   }
//   cout<<"all done with AllMinSpan"<<endl;
   return trees;
}

template<class T> DCGraph<T> *DCGraph<T>::GraphMaxHits(int nhits)
{

   if(nhits==0){
      return 0;
   }

   DCVertex<T> *max[100]={0};
   float mz[100]={0.};
   int N=nhits;
   if(nverts<N){
      N=nverts;
   }
   if(N>100){
      N=100;
   }
//   cout<<"Going to look for "<<N<<" biggest hits"<<endl;

   DCVertex<T> *vin=first;
   while(vin!=0){
      float x=0.;
      float y=0.;
      float z=0.;
      vin->GetData()->GetData(x,y,z);
      for(int i=0;i<N;i++){
	 if(z>mz[i]){
	    for(int j=N-1;j>i;j--){
	       max[j]=max[j-1];
	       mz[j]=mz[j-1];
	    }
	    max[i]=vin;
	    mz[i]=z;
	    break;
	 }
      }
      vin=vin->GetNextVertex();
   }

   DCGraph<T> g;
   for(int i=0;i<N;i++){
      g.AddVertex(max[i]->GetData());
//      cout<<"Adding vertex to new graph"<<endl;
//      max[i]->GetData()->Print();
   }

   return g.FindMinSpan();
      
}
  
template<class T> vector<float> DCGraph<T>::GetAllWeights()
{
   multiset<float> w;
   DCVertex<T> *v = first;
   DCEdge<T> *e;
   while(v!=0){
//      cout<<"on vertex "<<endl;
//     v->Print();
      e = v->GetFirstEdge();
      while(e!=0){
	 w.insert(e->GetMetric());
//	 cout<<"current edge metric is "<<e->GetMetric()<<endl;
	 e = e->GetNext();
//	 if(e!=0){
//	    e->Print();
//	 }
//	 else{
//	    cout<<"edge zero, going on to next vertex"<<endl;
//	 }
      }
      v=v->GetNextVertex();
   }

   multiset<float>::iterator it(w.begin());
   vector<float> vw;
   while(it!=w.end()){
//      cout<<*it<<endl;
      vw.push_back(*it);
      it++;
   }

   return vw;

}


template<class T> vector<float> DCGraph<T>::GetAllZ()
{
   multiset<float> w;
   DCVertex<T> *v = first;
   while(v!=0){
      float x,y,z;
      v->GetData()->GetData(x,y,z);
      w.insert(z);
      v=v->GetNextVertex();
   }

   multiset<float>::iterator it(w.begin());
   vector<float> vw;
   while(it!=w.end()){
//      cout<<*it<<endl;
      vw.push_back(*it);
      it++;
   }

   return vw;

}



#endif //DCGRAPH_H

