#include "PMT.h"

#include "MessageService/MsgService.h"

#include <cstdio>

int PMT::printit=-1;

void PMT::Init(int id)
{
        if(printit<0)
        {
                if(MsgService::Instance()->IsActive("XTalkFilter",
                                           Msg::kDebug))printit=1;
                else printit=0;
        }

        this->id=id;
        pixels.clear();
        nstrips=0;
        status=0;
}

PMT::~PMT(){}

std::vector<int>  PMT::GetPixelVector()
{
        std::vector<int>p; 
        for(unsigned int i=0;i<pixels.size();i++)
                p.push_back(pixels[i].GetId());
        return p;
}

Pixel * PMT::GetPixel(int id,int create)
{
        for(unsigned int i=0;i<pixels.size();i++)
        {
                pixels[i].SetPMT(this);
                if(pixels[i].GetId()==id)return &(pixels[i]);
        }
        
        if(!create)return 0;
//if its not there... create it
        Pixel newpixel(id);
        newpixel.SetPMT(this);
        pixels.push_back(newpixel);
        return & (pixels[pixels.size()-1]);
}

void PMT::Dump()
{
        if(nstrips<2)return;
        printf("PMT %d with %d strips\n",id,nstrips);
        for(unsigned int i=0;i<pixels.size();i++)pixels[i].Dump();
}
