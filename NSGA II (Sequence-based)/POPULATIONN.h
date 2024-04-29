# ifndef _POPULATIONN_H_
# define _POPULATIONN_H_


#include "INDIVIDUAL.h"


using namespace std;

class POPULATIONN
{
	public:
		INDIVIDUAL **ind;
		int pop;

		POPULATIONN(int size, int geneSize, int alt)
		{
			pop = size;

			this->geneSize = geneSize;

			ind = (INDIVIDUAL**)malloc(pop*sizeof(INDIVIDUAL*));

			
			if(alt)
			{
				for(int i = 0; i < pop; i++)
					ind[i] = new INDIVIDUAL(geneSize, size);
			}
		}
		

		void renew_pop()
		{
			int i, j, k;
			int count = 0;
			int *gen1, *gen2;
			INDIVIDUAL *toDel;

			for(i = 0; i < pop-1; i++)
			{
				gen1 = ind[i]->get_gene();
				for(j = i+1; j < pop; j++)
				{
					gen2 = ind[j]->get_gene();

					for(k = 1; k < geneSize+1; k++)
					{
						if(gen1[k]==gen2[k]){count++;}
						else if(gen1[k] == 0){break;}
					}

					if(count == k-1)
					{
						toDel = ind[j];
						ind[j] = new INDIVIDUAL(geneSize, pop);

						toDel->kill_indv_data();
						free(toDel);
					}

					count = 0;
				}
			}
		}
	private:
		int geneSize;

};

#endif