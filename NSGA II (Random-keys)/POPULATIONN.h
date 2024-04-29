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
			srand(time(NULL));
			pop = size;

			ind = (INDIVIDUAL**)malloc(pop*sizeof(INDIVIDUAL*));

			
			if(alt)
			{
				for(int i = 0; i < pop; i++)
					ind[i] = new INDIVIDUAL(geneSize, size);
			}
		}
		
	private:

};

#endif