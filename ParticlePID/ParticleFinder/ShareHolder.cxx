#include "MessageService/MsgService.h"



#include "ShareHolder.h"



ShareHolder::ShareHolder()
{
	hits.clear();
}

ShareHolder::~ShareHolder()
{

}



void ShareHolder::Reset()
{
	hits.clear();
}


void ShareHolder::BulkInsert(std::vector<Particle3D *> &p3d)
{
	Reset();
	
	for(unsigned int i=0;i<p3d.size();i++)
	for(int j=0;j<p3d[i]->entries;j++)
			Insert(p3d[i]->chain[j],p3d[i]->chainhit[j],p3d[i]->e[j]);
		
	for(unsigned int i=0;i<p3d.size();i++)
	for(int j=0;j<p3d[i]->entries;j++)	
		if(GetNumShared(p3d[i]->chain[j],p3d[i]->chainhit[j])>1)
		p3d[i]->SetShared(p3d[i]->chain[j],p3d[i]->chainhit[j]);
			
	
}


void ShareHolder::Insert(int chain, int chainhit,  double e)
{
	int i=find(chain,chainhit);
	if(i>-1)
	{
		//hits[i].e+=e;
		hits[i].nshared++;
	}else{
		hit h;
		h.e=e;
		h.chain=chain;
		h.chainhit=chainhit;
		h.nshared=1;
		hits.push_back(h);
	}
}

int ShareHolder::GetTotRemaining()
{
	int cnt=0;
	for(unsigned int i=0;i<hits.size();i++)
		if(hits[i].nshared>1)cnt++;
	return cnt;
}

int ShareHolder::GetNumShared(int chain, int chainhit)
{
	int i=find(chain,chainhit);
	if(i>-1)
		return hits[i].nshared;
	else
		return 0;
}

double ShareHolder::GetEShared(int chain, int chainhit)
{
	int i=find(chain,chainhit);
	if(i>-1)
		return hits[i].e;
	else
		return 0;
}

void ShareHolder::Take(int chain, int chainhit, double e)
{

	//printf("taking hits from %d %d %f    ",chain,chainhit,e);

	int i=find(chain,chainhit);
	if(i>-1)
	{
	
		hits[i].e-=e;
		hits[i].nshared--;
		
	//	printf("   now sharing %d\n",hits[i].nshared);	
		
	}
}	


int ShareHolder::find(int chain, int chainhit)
{
	int idx=-1;
	for(int i=0;i<(int)hits.size();i++)
		if(hits[i].chain==chain && hits[i].chainhit==chainhit )
		{
			idx=i;
			break;
		}
		
	return idx;
}



void ShareHolder::dump()
{
return;

	printf("\ndumping shareholder\n");
	for(int i=0;i<(int)hits.size();i++)
		printf("%d %d  #%d  e %f\n",hits[i].chain, hits[i].chainhit, hits[i].nshared,hits[i].e);
	printf("\n");
}
