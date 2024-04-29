# ifndef _INDIVIDUAL_H_
# define _INDIVIDUAL_H_

#include <iostream>
#include <string>
#include <cstring>
#include <string.h>
#include <ctype.h>
#include <vector>
#include <map>
#include <math.h> 
#include <ctime>
#include <cfloat>
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "POPULATIONN.h"

using namespace std;

class INDIVIDUAL
{
	public:
		//For a initial individual
		INDIVIDUAL(int genSize, int pop)
		{
			sizePop = pop;

			geneSize = genSize;
			dominates  = 0;
			dominators  = 0;

			obj[0] = 999999;
			obj[1] = 999999;
			obj[2] = 0;

			gen_rand_keys(pop);

			decode_gen();
		}

		//For a individual resulting from a crossover between two others
		INDIVIDUAL(int genSize, int pop, int theRandKeys[])
		{
			sizePop = pop; 

			geneSize = genSize;
			dominates  = 0;
			dominators  = 0;

			obj[0] = 999999;
			obj[1] = 999999;
			obj[2] = 0;

			set_randKeys(theRandKeys);

			decode_gen();
		}

		//For copy a oldest individual
		INDIVIDUAL(INDIVIDUAL* indv)
		{

			sizePop = indv->get_pop();
			geneSize = indv->get_gene_size();
			rank = indv->get_rank();
			crowdDist = indv->get_crowd_dist();

			dominates  = 0;
			dominators  = indv->get_dominators();

			obj[0] = indv->get_obj1();
			obj[1] = indv->get_obj2();
			obj[2] = indv->get_obj3();

			gene = (int*)malloc(sizeof(int)*(geneSize+1));
			randKeys = (int*)malloc(sizeof(int)*(geneSize+1));

			copy_keys(indv->get_rand_keys());
			copy_gen(indv->get_gene());
		}


		/**
		 * 	If returns:
		 * 		. 1 : The challenger dominates the challenged;
		 * 		.-1 : The challenger is dominateds by the challenged;
		 * 		. 0 : There is no dominance between them.
		 * */
		//int check_dominance(INDIVIDUAL *challenger)
		int check_dominance(float cObj1, float cObj2, float cObj3)
		{
			// FO Min cost Min slack days MAX visited cities


			/**/
			if(cObj1 < obj[0] || cObj3 > obj[2])
			{
				if(cObj1 <= obj[0] && cObj3 >= obj[2])
				{
					return 1;
				}
			}else if(obj[0] < cObj1 || obj[2] > cObj3)
			{
				if(obj[0] <= cObj1 && obj[2] >= cObj3)
				{
					return -1;
				}
			}

			return 0;
		}

		void increases_dominateds(){dominates++;}


		void increases_dominators(){dominators++;}

		/**
		 * 	If a mutation happens, then some new values will be assigned to the randKeys vector and 
		 * the previous gene sequence will be replaced by a new.  
		 * */
		void mutation(int pop, int interval, int chanceToChange, int lb)
		{
			int val; // drawn value
			int ver;
			int check;
			int i;
			int *auxRandKeys = (int*)malloc(sizeof(int)*geneSize);			

			ver = 0;

			for(i = 0; i < geneSize; i++)
			{
				check = rand()%interval;
				auxRandKeys[i] = randKeys[i];

				if(check > chanceToChange)
				{
					val = rand()%pop;

					check_att(i, lb, val);
				}

				if(randKeys[i] > 0)
					ver++;
			}


			if(ver < 2)
			{
				for(i = 0; i < geneSize; i++)
					randKeys[i] = auxRandKeys[i];
			}else
			{	
				free(gene);

				decode_gen();
			}

			obj[0] = 999999;
			obj[1] = 999999;
			obj[2] = 0;

			
			free(auxRandKeys);
		}

		void reset_dom()
		{
			dominates  = 0;
			dominators  = 0;
		}



		void set_crowd_dist(double crowd){crowdDist = crowd;}

		void set_obj1(float obj1){obj[0] = obj1;}

		void set_obj2(float obj2){obj[1] = obj2;}

		void set_obj3(float obj3){obj[2] = obj3;}

		void set_rank(int r){rank = r;}




		double get_crowd_dist(){return crowdDist;}
	
		int get_dominators(){return dominators;}

		int* get_gene(){return gene;}
	
		int get_gene_size(){return geneSize;}

		float get_obj1(){return obj[0];}

		float get_obj2(){return obj[1];}

		float get_obj3(){return obj[2];}

		int get_pop(){return sizePop;}

		int* get_rand_keys(){return randKeys;} 

		int get_rank(){return rank;}

		void kill_indv_data()
		{
			free(randKeys);
			//if(obj[2] > 0)
			free(gene);
		}

	private:
		//................................................................................................................
		//............................................... Variables  .....................................................
		//................................................................................................................
		/**
		 * How many individuals this one dominates.
		 * */
		int dominates; 
		
		/**
		 * How many individuals dominates this one.
		 * */
		int dominators; 

		/**
	  	* 1 -> feasible
	  	* 2 -> not feasible
	  	*/
		int feasible;  

		double crowdDist;


		/**
		 * Order in which cities should be visited.
		 * */
		int *gene;
		
		/**
		 * Gene size. Number of cities in the instance.	
		 * */
		int geneSize;
		
		/**
		 * 1st - Min trip cost;	
		 * 2nd - Min used slack time 	
		 * 3rd - Max number of visited cities
		 * */
		float obj[3]; 
		
		int sizePop;

		/**
		 * 	Randomly generated keys to be used as a reference for ordering the cities to be visited
		 * */
		int *randKeys;		


		int rank;

		//................................................................................................................
		//................................................................................................................
		//................................................................................................................



		void check_att(int pos, int lb, int val)
		{
			if(val <= lb){randKeys[pos] = -val;}
			else{randKeys[pos] = val;}
		}		

		void copy_gen(int newGen[])
		{
			int i;

			for(i = 0; i < geneSize; i++)
				gene[i] = newGen[i];
		}

		void copy_keys(int newKeys[])
		{
			int i;

			for(i = 0; i < geneSize; i++)
				randKeys[i] = newKeys[i];
		}

		/**
		 * According to randKeys it registers the sequence of cities to be visited (gene).
		 * */
		void decode_gen()
		{
			int aux;

			int *auxRandKeys = (int*)malloc(sizeof(int)*(geneSize));

			gene = (int*)malloc(sizeof(int)*(geneSize+1));




			for(int i = 0; i < geneSize; i++)
				auxRandKeys[i] = randKeys[i];
			

			for(int i=0; i < geneSize; i++)
			{
				if(auxRandKeys[i] < 0){gene[i] = -i;}
				else{gene[i]=i;}
			}

			gene[geneSize] = geneSize;


			sort(auxRandKeys);

			for(int i = geneSize-1; ; i--)
			{
				if(gene[i-1] < 0)
				{
					aux = gene[i];
					gene[i] = gene[i-1];
					gene[i-1] = aux;
				}else{break;}
			}
		}


		/**
		 * 	Generates the random keys to be used
		 * */
		void gen_rand_keys(int pop)
		{
			/**
			 *  From the population size (pop), 0.2% must be negative values and 0.8% nonnegative.  
			 * */
			int lb = (int)(0.6*pop); // lower bound
			int val; // drawn value
			int i, ver;

			randKeys = (int*)malloc(sizeof(int)*(geneSize));

			do
			{
				ver = 0;
				for(i = 0; i < geneSize; i++)
				{
					val = rand()%pop;

					check_att(i, lb, val);
				}

				for(i = 0; i <geneSize; i++)
				{
					if(randKeys[i] > 0)
						ver++;
				}
			}while(ver < 2);


		}

		void set_randKeys(int theRandKeys[])
		{
			randKeys = (int*)malloc(sizeof(int)*(geneSize));
			int i, ver = 0;

			for(int i = 0; i < geneSize; i++){
				randKeys[i] = theRandKeys[i];
				if(randKeys[i] > 0)
					ver++;
			}

			if(ver < 2)
				gen_rand_keys(sizePop);

			
		}		

		void sort(int auxRandKeys[]) /** Heap Sort*/
		{
			int i, holderG, holderS;

			for(i = (geneSize/2) ; i >=1; i--)
				sub_sort(auxRandKeys, i, geneSize-1);

			for(i = geneSize-1; i>=2; i--)
			{
				holderG = auxRandKeys[1];
				holderS = gene[1];

				auxRandKeys[1] = auxRandKeys[i];
				gene[1] = gene[i];

				auxRandKeys[i] = holderG;
				gene[i] = holderS;

				sub_sort(auxRandKeys, 1, i-1);
			}
		}

		void sub_sort(int auxRandKeys[], int root, int end)
		{
			int ready, son, holderG, holderS;

			ready = 0;
			int verify = root*2;


			while( ( verify <= end ) && (!ready) )
			{
				if(verify == end || (auxRandKeys[verify] < auxRandKeys[verify+1]) ){son = verify;}
				else {son = verify + 1;}

				if(auxRandKeys[root] > auxRandKeys[son])
				{
					holderG = auxRandKeys[root];
					holderS = gene[root];

					auxRandKeys[root] = auxRandKeys[son];
					gene[root] = gene[son];			

					auxRandKeys[son] = holderG;
					gene[son] = holderS;

					root = son;
					verify = root*2;
				}else{ready=1;}
			}
		}
};

#endif