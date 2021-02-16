#include "NueAna/ParticlePID/ParticleFinder/Reporter/ParticleReportHelper.h"


ParticleReportHelper::ParticleReportHelper()
{

//	sorter_map.clear();
}

ParticleReportHelper::~ParticleReportHelper()
{
}


void ParticleReportHelper::addtruth(ParticleTruthObject *h)
{

	truth_map[h->plane][h->strip][h->pid]+=h->sumenergy;

}


void ParticleReportHelper::addfound(ParticleObject *h)
{

	for(int i=0;i<h->plane.size();i++)
{		found_map[h->plane[i]][h->strip[i]][h->particle_id]+=h->energy[i];


	printf("adding particle to record matcher p,s,pid %d %d %d\n",h->plane[i],h->strip[i],h->particle_id);
}
}




void ParticleReportHelper::Process(ParticleObjectHolder & p)
{

	//iterate over map, and for each entry, make a new entry to add in p


	std::map<int, std::map<int, std::map<int, double> > >::iterator it1;
	std::map<int, std::map<int, double> > ::iterator it2;
	std::map<int, double>::iterator it3;


        std::map<int, std::map<int, double> > found_maxe_map;
        std::map<int, std::map<int, double> > found_tote_map;
	std::map<int, std::map<int, int> > found_maxpid_map;

	for(it1=found_map.begin();it1!=found_map.end();it1++)
	for(it2=it1->second.begin();it2!=it1->second.end();it2++)
	{
		double found_tote=0;
		double found_maxe=0;
		int found_maxpid=0;

		for(it3=it2->second.begin();it3!=it2->second.end();it3++)
		{
			found_tote+=found_map[it1->first][it2->first][it3->first];

			if (found_map[it1->first][it2->first][it3->first]>found_maxe)
			{
				found_maxe = found_map[it1->first][it2->first][it3->first];
				found_maxpid = it3->first;
			}
		}	

		int plane = it1->first;
		int strip = it2->first;
	
		found_maxe_map[plane][strip]=found_maxe;
                found_tote_map[plane][strip]=found_tote;
                found_maxpid_map[plane][strip]=found_maxpid;
	}


        std::map<int, std::map<int, double> > true_maxe_map;
        std::map<int, std::map<int, double> > true_tote_map;
        std::map<int, std::map<int, int> > true_maxpid_map;

        for(it1=truth_map.begin();it1!=truth_map.end();it1++)
        for(it2=it1->second.begin();it2!=it1->second.end();it2++)
        {
                double true_tote=0;
                double true_maxe=0;
                int true_maxpid=0;

                for(it3=it2->second.begin();it3!=it2->second.end();it3++)
                {
                        true_tote+=truth_map[it1->first][it2->first][it3->first];

                        if (truth_map[it1->first][it2->first][it3->first]>true_maxe)
                        {
                                true_maxe = truth_map[it1->first][it2->first][it3->first];
                                true_maxpid = it3->first;
                        }
                }

                true_maxe_map[it1->first][it2->first]=true_maxe;
                true_tote_map[it1->first][it2->first]=true_tote;
                true_maxpid_map[it1->first][it2->first]=true_maxpid;
        }


	//now iterate over found/true maps to get stats


	int totstrip=0;
	int matched_strips=0;

	double tote=0;
	double matched_e=0;

	for(it2=found_maxe_map.begin();it2!=found_maxe_map.end();it2++)
        for(it3=it2->second.begin();it3!=it2->second.end();it3++)
	{
		totstrip++;
		tote+=found_tote_map[it2->first][it3->first];
		
		if(found_maxpid_map[it2->first][it3->first]==true_maxpid_map[it2->first][it3->first])
		{
			matched_strips++;
			matched_e+=found_maxe_map[it2->first][it3->first];
		}

	}

   /*     ParticleReportObject *d = (ParticleReportObject*)p.particlematch->New(p.particlematch->GetEntries());


	d->totstrips=totstrip;
	d->matched_strips=matched_strips;
	d->tote=tote;
	d->matched_e=matched_e;
*/

}


