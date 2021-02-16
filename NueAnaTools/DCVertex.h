#ifndef DCVERTEX_H
#define DCVERTEX_H

#include <iostream>
#include <vector>
#include "TObject.h"
#include "NueAna/NueAnaTools/DCEdge.h"

using std::vector;
using std::cout;
using std::endl;

template<class T> class DCGraph;

//for a class T to work in this template it must:
//have a < operater
//have a == operator
//have a - operator that returns a T
//have an abs() operator
//have a Print() method
//have a Draw() method
//have a SetData() method
//have a GetData() method

template <class T> class DCVertex : public TObject 
{
public:

   DCVertex(){data=0; edges=0; next=0;}
   DCVertex(T *idata, DCVertex<T> *nextVert) { data=new T(*idata); next=nextVert; edges=0;}
   
   virtual ~DCVertex(){ 
      if(data!=0){delete data; data=0;}
      if(next!=0){delete next; next=0;}
      if(edges!=0){delete edges; edges=0;}
   }

   T *GetData(){return data;}
   DCVertex<T> *GetNextVertex() {return next;}
   DCEdge<T> *GetFirstEdge() {return edges;}
   void ResetNextVertex(DCVertex<T> *newvert) { next=newvert;}
   void ResetFirstEdge(DCEdge<T> *newedge) { edges=newedge;}

   virtual void Print(Option_t *t="") const;
   void PrintEdges();
   void PrintGraph();
   float ConnectTo(DCVertex<T> *);
   virtual void Draw(Option_t* t="");
   void DrawEdges(Option_t *t="");
   void DrawGraph(Option_t *t="");
   void ClearConnections();

private:
   T *data;
   DCEdge<T> *edges;
   DCVertex<T> *next;
   ClassDef(DCVertex<T>,1)

};


ClassImpT(DCVertex,T)

template<class T> void DCVertex<T>::PrintEdges() 
{
   if (edges == 0){
      cout << "vertex ";
      GetData()->Print();
      cout<< " has no edges" << endl;
   }
  else {
    cout << "vertex ";
    GetData()->Print();
    cout << " has edges to: " << endl;
    edges->Print();
  }
}  

// prints the edges of the vertex, calls the printGraph method for 
// the next vertex
template<class T>void DCVertex<T>::PrintGraph()
{
   PrintEdges();
   if (next != 0){
      next->PrintGraph();
   }
}

template<class T>void DCVertex<T>::Draw(Option_t *t)
{
//   cout<<"in DCVertex::Draw on ";
//   data->Print();
//   cout<<endl;
   data->Draw(t);
}

template<class T>void DCVertex<T>::DrawEdges(Option_t *t)
{
//   cout<<"In DCVertex::DrawEdges "<<endl;
   data->Draw(t);
   if(edges!=0){
      edges->DrawThisEdge(this,t);
   }
}
  
template<class T> void DCVertex<T>::DrawGraph(Option_t *t)
{
//   cout<<"In DCVertex::DrawGraph"<<endl;
   DrawEdges(t);
   if (next != 0){
//      cout<<"Drawing next vertex"<<endl;
      next->DrawGraph(t);
   }  
}


template<class T> void DCVertex<T>::Print(Option_t *t) const
{
   data->Print(t);
}

// adds an edge to connect the vertex to the vertex pointed to by Vert
template<class T> float DCVertex<T>::ConnectTo(DCVertex<T> *Vert) 
{
//   cout<<"In DCVertex::ConnectTo"<<endl;
  // allocate memory for a new Edge, set its Vertex pointer to point 
  // to Vert, and its Edge pointer to point to the rest of edges 
   T minus = (*(Vert->GetData()))-(*data);
   float m=(*data-*(Vert->GetData())).abs();
   float xm=0.;
   float ym=0.; 
   float zm=0.;
   minus.GetData(xm,ym,zm);
   float wm = m*zm;
   float xwm=xm*zm;
   float ywm=ym*zm;
   DCEdge<T> *newEdge = new DCEdge<T>(Vert, edges, m, xm, ym, wm, xwm, ywm);
   
   // make the new edge the first edge of the vertex
   edges = newEdge; 
   return newEdge->GetMetric();
}

template<class T> void DCVertex<T>::ClearConnections()
{
   edges=0;
}


#endif //VERTEX_H
