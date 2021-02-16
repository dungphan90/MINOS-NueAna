#ifndef DCHIT_H
#define DCHIT_H
#include <cmath>
#include "TObject.h"
class UberHit;

class DCHit : public TObject
{
public:
   DCHit();
   DCHit(const DCHit &d);
   DCHit(float zp, float tp, float s, float time);
   DCHit(UberHit *uh);

   ~DCHit();

   float GetZpos() const {return zpos;}
   float GetTpos() const {return tpos;}
   float GetPH() const {return pulseheight;}
   float GetTime() const {return time;}

   void Print(Option_t* option="") const;
   void Draw(Option_t* option="");
   void SetData(const char* cx, const char *cy, const char* cz="");
   void GetData(float &x, float &y, float &z) {x=drawx;y=drawy;z=drawz;}

   bool operator < (const DCHit &dc) const {
      if(zpos<dc.GetZpos()){
	 return true;
      }
      else if(zpos==dc.GetZpos()){
	 if(pulseheight>dc.GetPH()){
	    return true;
	 }
	 else if(pulseheight==dc.GetPH()){
	    if(tpos<dc.GetTpos()){
	       return true;
	    }	    
	 }
      }
      return false;
   }
	
   DCHit operator - (const DCHit &dc) const {
      DCHit newhit;
      newhit.zpos=zpos-dc.zpos;
      newhit.tpos=tpos-dc.tpos;
      if(pulseheight+dc.pulseheight!=0){
	 newhit.pulseheight=2*((pulseheight-dc.pulseheight)/
			       (pulseheight+dc.pulseheight));
      }
      else{
	 newhit.pulseheight=0.;
      }

      newhit.time=time-dc.time;
      newhit.drawx=drawx-dc.drawx;
      newhit.drawy=drawy-dc.drawy;
      if(drawz+dc.drawz!=0){
	 newhit.drawz=2*(drawz-dc.drawz)/(drawz+dc.drawz);
      }
      else{
	 newhit.drawz=0.;
      }
      return newhit;
   }

   DCHit operator + (const DCHit &dc) const {
      DCHit newhit;
      newhit.zpos=zpos+dc.zpos;
      newhit.tpos=tpos+dc.tpos;
      newhit.pulseheight=(pulseheight+dc.pulseheight);
      if(pulseheight+dc.pulseheight!=0){
	 newhit.time=((pulseheight*time+dc.pulseheight*dc.time)/
		      (pulseheight+dc.pulseheight));
      }
      else{
	 newhit.time=0;
      }
      //note, get data will work the same way for the new DCHit as it did
      //for the first dchit (not the added one)
      newhit.drawx=drawx+dc.drawx;
      newhit.drawy=drawy+dc.drawy;
      //this is going to get you into trouble.  Need a better way to 
      //set what data is going to be assigned to these variables.
      newhit.drawz=(drawz+dc.drawz);
      return newhit;
   }

   float abs() {return sqrt(GetZpos()*GetZpos()+GetTpos()*GetTpos());}

   bool operator== (const DCHit &dc) const {
     if((zpos==dc.zpos)&&
	(tpos==dc.tpos)&&
	(pulseheight==dc.pulseheight)&&
	(time==dc.time)) return true;
     return false;
   }


private:
   float zpos;
   float tpos;
   float pulseheight;
   float time;

   float drawx;
   float drawy;
   float drawz;

   ClassDef(DCHit,1)
};
#endif //DCHIT_H
