#include "NueAna/ParticlePID/Analysis/TH1DStat.h"
#include "math.h"

TH1DStat::TH1DStat(TH1D *t)
{
	th1d = t;
	
	td.Set(t->GetXaxis()->GetNbins()+2);
	td.Reset();

	pw.Set(t->GetXaxis()->GetNbins()+2);
	pw.Reset();
		
	wasErrorAdjusted=false;
}

//value, weight, and scaling for statistical error (POTweight of this event)
void TH1DStat::Fill(double v, double w, double s)
{
	int b = th1d->Fill(v,w);

	if(b<0)
	{
		b=0; //underflow
		if(v>th1d->GetXaxis()->GetXmax())
			b=pw.GetSize()-1; //overflow
	}
	

	td[b]++;
	pw[b]+=s;
	wasErrorAdjusted=false;
	
}

void TH1DStat::Dump()
{
	for(int i=0;i<td.GetSize();i++)
	{
		printf("%d %d %f\n",i,(int)td.At(i),pw.At(i));
	}

}

TH1D * TH1DStat::GetTH1D(int withstats)
{
	if(!withstats || wasErrorAdjusted)return th1d;
	
	//put in statistical errors!
	for(int i=0;i<td.GetSize();i++)
	{
		th1d->SetBinError(i,sqrt(td.At(i))*pw.At(i));
	}
	wasErrorAdjusted=true;
	return th1d;
}


