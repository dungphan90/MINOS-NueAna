#include <iostream>
#include "TMarker.h"
#include "NueAna/NueAnaTools/DCVertex.h"
#include "NueAna/NueAnaTools/DCEdge.h"
#include "NueAna/NueAnaTools/DCGraph.h"
#include "NueAna/NueAnaTools/DCHit.h"
#include "CalDetDST/UberHit.h"

using namespace std;

ClassImp(DCHit)

template class vector<DCHit>;
template class DCEdge<DCHit>;
template class DCVertex<DCHit>;
template class DCGraph<DCHit>;

DCHit::DCHit():
   zpos(0),
   tpos(0),
   pulseheight(0),
   time(0),
   drawx(0.),
   drawy(0.),
   drawz(0.)
{}

DCHit::DCHit(const DCHit &d):
  TObject::TObject(),
  zpos(d.zpos),
  tpos(d.tpos),
  pulseheight(d.pulseheight),
  time(d.time),
  drawx(d.drawx),
  drawy(d.drawy),
  drawz(d.drawz)
{}

DCHit::DCHit(float zp, float tp, float s, float t):
   zpos(zp),
   tpos(tp),
   pulseheight(s),
   time(t),
   drawx(0.),
   drawy(0.),
   drawz(0.)
{}

DCHit::DCHit(UberHit *uh)
{
   zpos=uh->GetPlane()*5.94;
   tpos=uh->GetStrip()*4.1;
   pulseheight=uh->GetPosMIP()+uh->GetNegMIP();
   time=((uh->GetPosMIP()*uh->GetPosTime()+
	  uh->GetNegMIP()*uh->GetNegTime())/
	 (uh->GetPosMIP()+uh->GetNegMIP()));
   drawx=0.;
   drawy=0.;
   drawz=0.;
}

DCHit::~DCHit()
{}

bool operator< (const DCHit &left, const DCHit &right)
{
   if(left.GetZpos()<right.GetZpos()){
      return true;
   }
   else if(left.GetZpos()==right.GetZpos()){
      if(left.GetTpos()<right.GetTpos()){
	 return true;
      }
      else{
	 return false;
      }
   }
   else{
      return false;
   }
   return false;
}

void DCHit::Print(Option_t* /*option*/) const
{

   cout<<"Zpos "<<zpos<<" tpos "<<tpos
       <<" PulseHeight "<<pulseheight<<" Time "<<time<<endl;

}

   
void DCHit::Draw(Option_t* /*t*/)
{
   TMarker *m = new TMarker(1.*zpos,1*tpos,8);   
   m->Draw();
}

void DCHit::SetData(const char* cx, const char* cy, const char *cz)
{
   if((string)(cx)=="zpos"){
      drawx=1.*zpos;
   }
   else if((string)(cx)=="tpos"){
      drawx=1.*tpos;
   }
   else if((string)(cx)=="pulseheight"){
      drawx=pulseheight;
   }
   else if((string)(cx)=="time"){
      drawx=time;
   }
   else{
      cout<<"Error in set data, don't know option "<<cx<<endl;
      cout<<"Setting cx to zpos"<<endl;
      drawx=zpos;
   }

   if((string)(cy)=="zpos"){
      drawy=1.*zpos;
   }
   else if((string)(cy)=="tpos"){
      drawy=1.*tpos;
   }
   else if((string)(cy)=="pulseheight"){
      drawy=pulseheight;
   }
   else if((string)(cy)=="time"){
      drawy=time;
   }
   else{
      cout<<"Error in set data, don't know option "<<cy<<endl;
      cout<<"Setting cy to zpos"<<endl;
      drawy=zpos;
   }

   if((string)(cz)=="zpos"){
      drawz=1.*zpos;
   }
   else if((string)(cz)=="tpos"){
      drawz=1.*tpos;
   }
   else if((string)(cz)=="pulseheight"){
      drawz=pulseheight;
   }
   else if((string)(cz)=="time"){
      drawz=time;
   }
   else if((string)(cz)==""){
   }
   else{
      cout<<"Error in set data, don't know option "<<cz<<endl;
      cout<<"Setting cz to zpos"<<endl;
      drawz=zpos;
   }

}

