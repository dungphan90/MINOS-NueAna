#include "MessageService/MsgService.h"
#include "PixelSpot.h"
#include "Pixel.h"
#include "PMT.h"

#include <cstdio>

int PixelSpot::printit=-1;

PixelSpot::~PixelSpot(){}

PixelSpot::PixelSpot()
{
        Init(-1);
}
 
PixelSpot::PixelSpot(int id)
{
        Init(id);
}

void PixelSpot::SetPixel(Pixel * p)
{
        parentPixel=p;
}
 
Pixel * PixelSpot::GetPixel()
{
        return parentPixel;
}

void PixelSpot::Init(int id)
{
        if(printit<0)
        {
                if(MsgService::Instance()->IsActive("XTalkFilter",
                                           Msg::kDebug))printit=1;
                else printit=0;
        }

        this->id=id;
        strips.clear();
        q.clear();
        sum_q=0;
        n=0;
        strip_e.clear();
        sum_strip_e=0;
        spot_e=0;
        status=0;

}


void PixelSpot::TakeE(double e)
{
     //   if(printit)printf("pixelspot take %f\n",e);
        if(sum_strip_e<=0)return;
                
       // if(printit)printf("sum_strip_e %f\n",sum_strip_e);      
        for(unsigned int i=0;i<strip_e.size();i++)
        {
                double totake = e * strip_e[i]/sum_strip_e;
        
                if(printit)printf("pixelspot adjusting strip %d e from %f to %f\n",
                                 strips[i],strip_e[i],strip_e[i]-totake);
                strip_e[i]-=totake;
                q[i]-= + e*q[i]/sum_q;
        }

	//recalculate sum_q
	sum_q=0;
        sum_strip_e=0;
	for(unsigned int i=0;i<q.size();i++)sum_q+=q[i];
	for(unsigned int i=0;i<strip_e.size();i++)sum_strip_e+=strip_e[i];

// mark the change
        status=1;
        GetPixel()->SetStatus(1);
        GetPixel()->GetPMT()->SetStatus(1);
        parentPixel->AddE(-e);
        parentPixel->AddStripE(-e);
}

void PixelSpot::Dump()
{
        if(printit)printf("\t\t PixelSpot %d  n %d  e %f strips ",
                                                      id,n,sum_q);
        for(unsigned int i=0;i<strips.size();i++)
                       if(printit)printf(" %d %f",strips[i],q[i]);
        if(printit)printf("\n");
}

void PixelSpot::AddStrip(int strip, double q, double strip_e)
{
// if(printit)printf("adding strip with idx %d\n",strip);

        this->strips.push_back(strip);
        this->strip_e.push_back(strip_e);
        this->q.push_back(q);
        sum_q+=q;
        sum_strip_e+=strip_e;
        n++;
        parentPixel->AddE(q);
        parentPixel->AddStripE(strip_e);
        parentPixel->GetPMT()->IncStripCount();
}

void PixelSpot::AddE(double e)
{
        spot_e+=e;
        status=1;
        GetPixel()->SetStatus(1);
        GetPixel()->GetPMT()->SetStatus(1);
        for(unsigned int i=0;i<strips.size();i++)
        {
                strip_e[i]=strip_e[i] + e*strip_e[i]/sum_strip_e;
                q[i]=q[i] + e*q[i]/sum_q;

        }
        sum_q+=e;
        sum_strip_e+=e;
        parentPixel->AddE(e);
        parentPixel->AddStripE(e);
        
// if(printit)printf("adding to pixel with addr %p\n",parentPixel);

}
