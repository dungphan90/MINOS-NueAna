#ifndef DCEDGE_H
#define DCEDGE_H

#include <iostream>
#include "TLine.h"
#include "TText.h"
#include "TObject.h"

template<class T> class DCVertex;
template<class T> class DCGraph;

using std::cout;
using std::endl;

template<class T> class DCEdge : public TObject
{
public:
   DCEdge(){metric=0;end=0;next=0;}
   DCEdge(DCVertex<T> *Vert, DCEdge<T> *nextEdge, 
	  float m, float xm, float ym, float wm, float wxm, float wym)
      {metric=m; xmetric=xm, ymetric=ym; 
      wmetric=wm; xwmetric=wxm; ywmetric=wym;
      end=Vert; next=nextEdge;}

   virtual ~DCEdge() {delete next;}

   DCVertex<T> *GetEnd() {return end;}
   DCEdge<T> *GetNext() {return next;}
   float GetMetric() const {return metric;}
   float GetXMetric() const {return xmetric;}
   float GetYMetric() const {return ymetric;}

   float GetWeighMetric() const {return wmetric;}
   float GetXWeighMetric() const {return xwmetric;}
   float GetYWeighMetric() const {return ywmetric;}

   virtual void Print(Option_t* option="") const;
   void DrawThisEdge(DCVertex<T> *vert,Option_t* option="");

private:
   float metric;
   float xmetric;
   float ymetric;
   float wmetric;
   float xwmetric;
   float ywmetric;

   DCVertex<T> *end;
   DCEdge<T> *next;

   ClassDef(DCEdge<T>,1)
};

ClassImpT(DCEdge,T)

template<class T> void DCEdge<T>::Print(Option_t* option) const
{
   end->GetData()->Print(option);
   cout<<" Metric is "<<metric<<" X "<<xmetric<<" Y "<<ymetric<<endl;
   cout<<" Weighted Metric is "<<wmetric<<" WX "<<xwmetric<<" WY "<<ywmetric<<endl;
   if (next == 0)
      cout << endl;
   else {
      next->Print(option);
   }
}

template<class T> void DCEdge<T>::DrawThisEdge(DCVertex<T> *vert, Option_t* t)
{
//   cout<<"in DCEdge::DrawThisEdge"<<endl;
//   cout<<"Drawing Edge from ";
//   vert->Print();
//   cout<<" to ";
//   end->Print();
//   cout<<endl;
   end->GetData()->Draw();

   float x1=0.;
   float y1=0.;
   float z1=0.;
   vert->GetData()->GetData(x1,y1,z1);
   
   float x2=0.;
   float y2=0.;
   float z2=0.;
   end->GetData()->GetData(x2,y2,z2);
   TLine *l = new TLine(x1,y1,x2,y2);
   l->Draw();

   float xtext = (x1+x2)/2.;
   float ytext = (y1+y2)/2.;
   char dist[10];
   sprintf(dist,"d=%.1f",metric);
   TText *txt = new TText();
   txt->SetTextSize(.02);
   txt->DrawText(xtext,ytext,dist);

   if(next!=0){
//      cout<<"going on to next edge"<<endl;
      next->DrawThisEdge(vert,t); 
   }

}

#endif //DCEDGE_H
