#include "Pixel.h"
#include "PixelSpot.h"
#include "PMT.h"

#include "MessageService/MsgService.h"

#include <cstdio>

int Pixel::printit=-1;

Pixel::Pixel()
{
        id=-1;
}

Pixel::Pixel(int id)
{
        if(printit<0)
        {
                if(MsgService::Instance()->IsActive("XTalkFilter",
                                           Msg::kDebug))printit=1;
                else printit=0;
        }

        this->id=id;
        parentPMT=0;
        stripE=0;
        e=0;
        pixelspots.clear();
        status=0;
        type=-1;
        hopstomixed=-1;
}

Pixel::~Pixel(){}

std::vector<int>  Pixel::GetPixelSpotVector()
{
        std::vector<int> ps; 
        for(unsigned int i=0;i<pixelspots.size();i++)  
                ps.push_back(pixelspots[i].GetId()); 
        return ps;
}


void Pixel::TakeE(double ee)
{

     //   if(printit)printf("pixel take %f\n",ee);
        if(e<=0)return;
        
// take the energy from the spots as well by proportion
// that they contribute
        
        for(unsigned int i=0;i<pixelspots.size();i++)
        {
                if(e<1e-5)break; //already took the energy
                double totake = ee * pixelspots[i].GetE()/e;
                pixelspots[i].TakeE(totake);            
       //         if(printit)printf("call pixelspot take %f from %d\n",
         //                         totake,i);
        }

//      if(ee<=e)e=0; 
//      else e-=ee; 

}

PixelSpot * Pixel::GetPixelSpot(int id, int create)
{
        for(unsigned int i=0;i<pixelspots.size();i++)
        {
                if(pixelspots[i].GetId()==id)
                {
// make sure the address is correct...somehow it can get shifted!
                        pixelspots[i].SetPixel(this);
                        return &(pixelspots[i]);
                }
        }
        
        if(!create)return 0;
        
// if its not there... create it
        PixelSpot newpixelspot(id);
        newpixelspot.SetPixel(this);
        
// if(printit)printf("creating pixelspot id %d with pixel addr %p\n",id,this); 
        pixelspots.push_back(newpixelspot);
        return & (pixelspots[pixelspots.size()-1]);
}


void Pixel::Dump()
{
        if(printit)printf("\t Pixel %d e %f\n",id,e);
        for(unsigned int i=0;i<pixelspots.size();i++)pixelspots[i].Dump();
}

int Pixel::GetNStrips()
{
        int n=0;
        for(unsigned int i=0;i<pixelspots.size();i++)
                n+=pixelspots[i].GetNStrips(); 
        return n;
}
