# ifndef _POPULATIONNTREATMENT_H_
# define _POPULATIONNTREATMENT_H_

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

#include "INDIVIDUAL.h"
#include "POPULATIONN.h"

#define INF 999999
#define ETAC 0.00000000001 //eta_c

using namespace std;

class POPULATIONNTREATMENT
{
	public:
		POPULATIONNTREATMENT(int pop, float cross)
		{
			lb = (int)pop*0.1; // 50 -> 10
			crossChance = cross;
		}

		/**
		 * 	Tournament routine 
		 * */
		void selection(POPULATIONN *oldOne, POPULATIONN *newOne, int pop)
		{
			srand(time(NULL));

			int *firstAux = (int*)malloc(sizeof(int)*pop);
			int *secondAux = (int*)malloc(sizeof(int)*pop);
			int i, randValIndex;

			for(i = 0; i < pop; i++)
			{
				firstAux[i] = i;
				secondAux[i] = i;
			}

			for(i = 0; i < pop-1; i++)
			{
				randValIndex = interval_rand(i, pop-1);
				aux_swap(i, randValIndex, firstAux);

				randValIndex = interval_rand(i, pop-1);
				aux_swap(i, randValIndex, secondAux);
			}

			/**
			 * ..........................................................
			 * ..........................................................
			 * */

			INDIVIDUAL* father;
			INDIVIDUAL* mother;

			for(i = 0; i < pop; i+=4)
			{
				if(i+1<pop){father = tournament(oldOne->ind[ firstAux[i]   ], oldOne->ind[ firstAux[i+1] ]);}
				else{father = oldOne->ind[ firstAux[i] ];}
				
				if(i+3<pop){mother = tournament(oldOne->ind[ firstAux[i+2] ], oldOne->ind[ firstAux[i+3] ]);}
				else if(i+2<pop){mother = oldOne->ind[ firstAux[i+2] ];}
				else if(i+1 < pop){mother = oldOne->ind[ firstAux[i+1] ];}
				else{mother = oldOne->ind[ firstAux[i] ];}

			
				crossover(father, mother, newOne, i, pop);
				
				//--------------------------------------------------------------------------------
				if(i+1<pop){father = tournament(oldOne->ind[ secondAux[i]   ], oldOne->ind[ secondAux[i+1] ]);}
				else{father = oldOne->ind[ secondAux[i] ];}
				
			
				if(i+3<pop){mother = tournament(oldOne->ind[ secondAux[i+2] ], oldOne->ind[ secondAux[i+3] ]);}
				else if(i+2<pop){mother = oldOne->ind[ secondAux[i+2] ];}
				else if(i+1 < pop){mother = oldOne->ind[ secondAux[i+1] ];}
				else{mother = oldOne->ind[ secondAux[i] ];}
				
			
				crossover(father, mother, newOne, i+2, pop);
			}



			free(firstAux);
			free(secondAux);
		}

		/**
		 * 									"fast-non-dominated-sort"
		 * Updates the dominance value of each member of the POPULATIONN and their respective ranks.
		 * */
		void classifies_population(POPULATIONN *population, int pop)
		{
			int i, j;
			int check;

			for(i = 0; i < pop; i++)
				population->ind[i]->reset_dom();

			for(i = 0; i < pop-1; i++)
			{
				for(j = i+1; j < pop; j++)
				{
					check = population->ind[i]->check_dominance(population->ind[j]->get_obj1(), population->ind[j]->get_obj2(), population->ind[j]->get_obj3());
					
					if(check==1)
					{
						population->ind[i]->increases_dominators();
						population->ind[j]->increases_dominateds();						
					}else if(check==-1)
					{

						population->ind[i]->increases_dominateds();
						population->ind[j]->increases_dominators();
					}
				}
			} 	
			/** 
			 * Sort by dominance  -- rever/ refazer
			 * */
			sort_pop(population, pop);


			/**
			 * Set the ranks
			 * */
			int record = 0;
			int rank = 1; /**
						  *	The rank 1 is the > best < classification
						  **/

			vector<int> qtdInEachRank; // Quantity of members of each rank
			int qtdInRank = 0;


			for(i = 0; i < pop; i++)
			{
				if(population->ind[i]->get_dominators() == record)
				{
					population->ind[i]->set_rank(rank);
					qtdInRank++;
				}else
				{
					qtdInEachRank.push_back(qtdInRank); // register the quantity of members of previous rank 
					record = population->ind[i]->get_dominators();
					rank++;

					population->ind[i]->set_rank(rank);
					qtdInRank = 1; // reset the counter
				}
			}


			qtdInEachRank.push_back(0);

			

			INDIVIDUAL *aux;
			int popStop = pop;
			int  k;

			for(i = 0; i < popStop-1; i++)
			{
				for(j = i+1; j < popStop; j++)
				{
					if((population->ind[i]->get_obj1() == population->ind[j]->get_obj1()) && (population->ind[i]->get_obj3() == population->ind[j]->get_obj3()) )
					{
						aux = population->ind[j];
						aux->set_rank(-1);

						for( k = j; k < pop-1; k++)
						{
							population->ind[k] = population->ind[k+1];
						}

						if( (population->ind[i]->get_obj1() == population->ind[j]->get_obj1()) && (population->ind[i]->get_obj3() == population->ind[j]->get_obj3()) )
							j--;

						population->ind[k] = aux;
						popStop--;
					}
				}
			}

			int old = 1;
			int reg = 0;
			int ck = 0;

			for(i = 0; i < pop; i++)
			{
				if(population->ind[i]->get_rank() == old)
				{
					reg++;
				}else if(population->ind[i]->get_rank() == -1)
				{
					if(!ck){
						old++;
						reg = 0;
						ck=1;
					}

					population->ind[i]->set_rank(old);
					reg++;
				}else
				{
					qtdInEachRank[old-1] = reg;
					reg = 0;
					old++;
				}
			}

			qtdInEachRank[old-1] = reg;


			/**
			 *************************************************************************** 
			 *                     assing_crowding_distance_indices
			 ***************************************************************************
			 **/
			int popCounter = 0, auxCounter; //Members already verified from the population / the next ones to be verified
			int qtdRanks = qtdInEachRank.size(); // how many ranks the population have
			//aqui
			double *crowd = (double*)malloc(sizeof(double)*pop); 

			double *obj1 = (double*)malloc(sizeof(double)*pop); 
			double *obj2 = (double*)malloc(sizeof(double)*pop);
			double *obj3 = (double*)malloc(sizeof(double)*pop);

			int dif;

			for(i = 0; i < pop; i++)
			{
				obj1[i] = population->ind[i]->get_obj1();
				obj2[i] = population->ind[i]->get_obj2();
				obj3[i] = -population->ind[i]->get_obj3();
			}

			//This loop can be reduced
			for(i = 0; i < qtdRanks; i++)
			{
				auxCounter = popCounter + qtdInEachRank[i];
				assing_crowding_distance_indices(population, popCounter, auxCounter, pop, obj1, obj2, obj3);

				dif = auxCounter - popCounter;

				for(j = popCounter; j < auxCounter; j++)
					crowd[j] = -population->ind[j]->get_crowd_dist(); // Changing the signal to do a descending sort

				if(dif > 2)
					sort_vect(crowd, crowd, crowd, 1, 0, 0, population, popCounter, auxCounter, 0);

				popCounter = auxCounter;
			}

			free(obj1);
			free(obj2);
			free(obj3);
			free(crowd);
		}

		
		void join_pop(POPULATIONN* mainPop, POPULATIONN* childPop, POPULATIONN *joinPop, int pop)
		{
			int i, j = pop;
			
			for(i = 0; i < pop; i++)
			{
				joinPop->ind[i] = mainPop->ind[i];
			}
			
			for(i = 0; i < pop; i++)
			{
				joinPop->ind[j] = childPop->ind[i];
				j++;
			}
		}
		
		void pop_mutation(POPULATIONN* population, int pop, int popMutInterval, int geneMutInterval)
		{
			srand(time(NULL));

			if(popMutInterval > 0 && geneMutInterval > 0)
			{
				int valRand = (rand()%popMutInterval)+1; // if % is 1 so the  "+1" must be removed.
				int randIndex;
				
				for(int i = 0; i < valRand; i++)
				{
					randIndex = rand()%pop;
					population->ind[randIndex]->mutation(pop, pop, geneMutInterval, lb);
				}
			}
		}

		void update_pop(POPULATIONN* mainPop, POPULATIONN *joinPop, int pop, int doublePop, float elite)
		{
			int sizeElite = (int)pop*elite;
			int i = sizeElite, j;
			
			//for(j = 0; j < sizeElite; j++)
			//	mainPop->ind[j] = new INDIVIDUAL(mainPop->ind[j]);

			for(i = 0; i < pop; i++)
			{
				mainPop->ind[i] = new INDIVIDUAL(joinPop->ind[i]);
			}

			for (i = 0; i < doublePop; i++)
			{
				joinPop->ind[i]->kill_indv_data();
				

				delete joinPop->ind[i];
			}
		}
		
		
	private:
		int lb;
		float crossChance;

		typedef struct sortpop
		{
			sortpop *next;
			sortpop *prev;

			INDIVIDUAL* idv;
		}auxSortPop;


		/**
		 * Sets received value to the range [-0.5*pop,0[ and ]0.5*pop,pop].
		 **/
		int adjust_value(int val, int pop, int lb)
		{
			if(val > pop){return adjust_value(val-pop, pop, lb);}
			if(val <= lb)
			{
				if(val > 0 || (-val <= lb)){return -val;}
				else{return adjust_value(val+lb, pop, lb);}
			}else{return val;}

		}


		int arithmetic(int father, int mother, int pop, int lb)
		{
			return adjust_value( ((int)(father+mother)/2), pop, lb);
		}

		/**
		 * 							"crowding-distance-assignment"
		 * 	
		 * */
		void assing_crowding_distance_indices(POPULATIONN *population, int start, int end, int pop, double *obj1, double *obj2, double *obj3)
		{
			int i, j;

			for( i = start; i < end; i++)
				population->ind[i]->set_crowd_dist(0);

			//por enquanto, vou considerar a ordenação já feita na classifies_POPULATIONN(). Vou pensar como arrumar isto e deixar organizado pra ambos

			int dif = end-start;

			if(dif == 1)
			{
				population->ind[start]->set_crowd_dist(INF);
			}else if(dif == 2)
			{
				population->ind[start]->set_crowd_dist(INF);
				population->ind[end-1]->set_crowd_dist(INF);
			}else
			{
				
				double minObj, maxObj;
				//.................................................
				for(i = 0; i < 3; i++)
				{
					switch(i)
					{
						case 0:
							sort_vect(obj1, obj2, obj3, 1, 0, 0, population, start, end, 1);
							minObj = obj1[start];
							maxObj = obj1[end-1];
							for(j = start+1; j < end-1; j++)
								population->ind[j]->set_crowd_dist( population->ind[j]->get_crowd_dist() + abs( (population->ind[j+1]->get_obj1() - population->ind[j-1]->get_obj1())/(maxObj-minObj) ) );
							
							break;

						case 1:
							//sort_vect(obj1, obj2, obj3, 0, 1, 0, population, start, end, 1);
							//minObj = obj2[start];
							//maxObj = obj2[end-1];
							//for(j = start+1; j < end-1; j++)
						//		population->ind[j]->set_crowd_dist( population->ind[j]->get_crowd_dist() + abs( (population->ind[j+1]->get_obj2() - population->ind[j-1]->get_obj2())/(maxObj-minObj) ) );
							break;

						case 2:
							sort_vect(obj1, obj2, obj3, 0, 0, 1, population, start, end, 1);
							minObj = obj3[end-1];
							maxObj = obj3[start];
							for(j = start+1; j < end-1; j++)
								population->ind[j]->set_crowd_dist( population->ind[j]->get_crowd_dist() + abs( (population->ind[j+1]->get_obj3() - population->ind[j-1]->get_obj3())/(maxObj-minObj) ) );
							break;
					}
				}
				
			}

		}


		void aux_swap(int i, int randValIndex, int vet[])
		{
			int temp;

			temp = vet[randValIndex];
			vet[randValIndex] = vet[i];
			vet[i] = temp;
		}


		auxSortPop* create_aux_sort(INDIVIDUAL* idv)
		{
			auxSortPop *aux = (auxSortPop*)malloc(sizeof(auxSortPop));
			aux->next = NULL;
			aux->prev = NULL;
			aux->idv = idv;

			return aux;
		}

		/**
		 * Crossover for REAL variables. -- Using the: 
		 * 			(1)- geometric mean -> root(p1*p2); 
		 * 			(0)- arithmetic mean-> (p1+p2)/2.
		 * */
		void crossover(INDIVIDUAL *father, INDIVIDUAL *mother, POPULATIONN *childPop, int index, int pop)
		{
			int check = interval_rand(0, pop)/pop;
			int size = father->get_gene_size();
			int *fatherKeys = father->get_rand_keys();
			int *motherKeys = mother->get_rand_keys();
			
			if(check <= crossChance)
			{				
				cross_steps(fatherKeys, motherKeys, childPop, index, index+1, pop, lb, size); // Second child
			}else
			{	
				childPop->ind[index] = new INDIVIDUAL(father);
				childPop->ind[index+1] = new INDIVIDUAL(mother);
			}
		}

		//https://gist.github.com/decheng-zhang/25910f9abaf6fbc42509ba854348cc50
		//sbx x-over
		void cross_steps(int *father, int *mother, POPULATIONN *childPop, int firstChild, int secondChild, int pop, int lb, int size)
		{
			int i;
			double rand;

			int y1, y2; // y1 -> first
			double upper = pop, lower = 0; // the the upper limit and lower limit for the variable; 
			double alpha, beta, betaq; // eta_c is a distribution index, a really small value
			double check;


			int* keysFirst = (int*)malloc(sizeof(int)*size);
			int* keysSecond = (int*)malloc(sizeof(int)*size);


			for(i = 0; i < size; i++)
			{
				check = interval_rand(0, pop)/pop;

				if(check <= 0.5)
				{
					if(father[i] < mother[i])
					{
						y1 = father[i];
						y2 = mother[i];
					}else
					{
						y1 = mother[i];
						y2 = father[i];
					}

					if(father[i] != mother[i])
					{
						//First Child
						rand = interval_rand(0, pop)/pop;

						beta = 1.0 + (2.0*(y1-lower)/(y2-y1));
						alpha = 2.0 - pow(beta, -(ETAC+1.0));

						if(rand <= 1/alpha){betaq = pow(rand*alpha, (1.0/(ETAC+1)));}
						else{betaq = pow(  (1/(2.0- rand*alpha)), (1.0/(ETAC+1)));}

						keysFirst[i] = (int)0.5*((y1+y2) - betaq*(y2-y1));

						//Second Child
						beta = 1.0 + (2.0*(upper-y2)/(y2-y1) );
						alpha = 2.0 - pow(beta, -(ETAC+1.0));

						if(rand <= 1/alpha){betaq = pow(rand*alpha, (1.0/(ETAC+1)));}
						else{betaq = pow(  (1/(2.0- rand*alpha)), (1.0/(ETAC+1)));}

						keysSecond[i] = (int)0.5*((y1+y2) - betaq*(y2-y1));
					}else
					{
						keysFirst[i] = father[i];
						keysSecond[i] = mother[i];
					}

				}else
				{
					keysFirst[i] = father[i];
					keysSecond[i] = mother[i];
				}
			}




			for(i = 0; i < size; i++)
			{
				keysFirst[i] = adjust_value(keysFirst[i], pop, lb);
				keysSecond[i] = adjust_value(keysSecond[i], pop, lb);
			}



			childPop->ind[firstChild] = new INDIVIDUAL(size, pop, keysFirst);
			childPop->ind[secondChild] = new INDIVIDUAL(size, pop, keysSecond);

			free(keysFirst);
			free(keysSecond);
		}

		int interval_rand(int base, int pop)
		{
			int interval = pop-base;

			return base + rand()%interval;	
		}

		/**
		 *  An ascending order will be made in relation to how many dominate them.
		 **/
		void sort_pop(POPULATIONN *population, int pop)
		{
			int i;
			INDIVIDUAL *holder, *auxHolder;

			//sort by dominance
			for(i = (int)(pop)/2; i > 0; i--)
				sub_sort_pop(population, i-1, pop);


			for(i = pop-1; i>=1; i--)
			{
				holder = population->ind[0];
				population->ind[0] = population->ind[i];
				population->ind[i] = holder;

				sub_sort_pop(population, 0, i-1);
			}


			auxSortPop *aux = create_aux_sort(population->ind[0]);
			auxSortPop *holdAux, *swapAux;


			for(i = 1; i < pop; i++)
			{
				holdAux = aux;
				
				do
				{
					if(holdAux->idv->get_dominators() == population->ind[i]->get_dominators())
					{
						if(holdAux->idv->get_obj1() > population->ind[i]->get_obj1())
						{

							if(holdAux->prev == NULL)
							{
								aux = create_aux_sort(population->ind[i]);
								aux->next = holdAux;
								holdAux->prev = aux;

								holdAux = NULL;
							}else
							{
								swapAux = holdAux->prev;
								swapAux->next = create_aux_sort(population->ind[i]);
								swapAux->next->prev = swapAux;
								swapAux->next->next = holdAux;
								holdAux->prev = swapAux->next;

								holdAux = NULL;
							}
						}else if(holdAux->idv->get_obj1() == population->ind[i]->get_obj1())
						{
							if(holdAux->idv->get_obj2() > population->ind[i]->get_obj2())
							{

								if(holdAux->prev == NULL)
								{
									aux = create_aux_sort(population->ind[i]);
									aux->next = holdAux;
									holdAux->prev = aux;

									holdAux = NULL;
								}else
								{
									swapAux = holdAux->prev;
									swapAux->next = create_aux_sort(population->ind[i]);
									swapAux->next->prev = swapAux;
									swapAux->next->next = holdAux;
									holdAux->prev = swapAux->next;

									holdAux = NULL;
								}


							}else if(holdAux->idv->get_obj1() == population->ind[i]->get_obj1())
							{
								if(holdAux->idv->get_obj3() <= population->ind[i]->get_obj3())
								{
									if(holdAux->prev == NULL)
									{
										aux = create_aux_sort(population->ind[i]);
										aux->next = holdAux;
										holdAux->prev = aux;

										holdAux = NULL;
									}else
									{
										swapAux = holdAux->prev;
										swapAux->next = create_aux_sort(population->ind[i]);
										swapAux->next->prev = swapAux;
										swapAux->next->next = holdAux;
										holdAux->prev = swapAux->next;

										holdAux = NULL;
									}


								}else if(holdAux->next == NULL)
								{
									holdAux->next = create_aux_sort(population->ind[i]);
									holdAux->next->prev = holdAux;

									holdAux = NULL;
								}else{holdAux = holdAux->next;}


							}else if(holdAux->next == NULL)
							{
								holdAux->next = create_aux_sort(population->ind[i]);
								holdAux->next->prev = holdAux;

								holdAux = NULL;
							}else{holdAux = holdAux->next;}

						}else if(holdAux->next == NULL)
						{
							holdAux->next = create_aux_sort(population->ind[i]);
							holdAux->next->prev = holdAux;

							holdAux = NULL;
						}else{holdAux = holdAux->next;}


					}else if(holdAux->next == NULL)
					{
						holdAux->next = create_aux_sort(population->ind[i]);
						holdAux->next->prev = holdAux;

						holdAux = NULL;
					}else{holdAux = holdAux->next;}

				}while(holdAux!=NULL);
			}

			holdAux = aux;

			for(i = 0; i < pop; i++)
			{
				population->ind[i] = holdAux->idv;
				swapAux = holdAux->next;
				free(holdAux);
				holdAux = swapAux;
			}
		}

		void sub_sort_pop(POPULATIONN *population, int root, int end)
		{
			INDIVIDUAL *holder;
			int ready = 0, son;
			int verify = root*2;

			while(verify < end && !ready)
			{
				if((population->ind[verify]->get_dominators() > population->ind[verify+1]->get_dominators()) ){son = verify;}
				else {son = verify + 1;}


				if( population->ind[root]->get_dominators() < population->ind[son]->get_dominators() )
				{
					holder = population->ind[root];
					population->ind[root] = population->ind[son];
					population->ind[son] = holder;

					root = son;
					verify = root*2;
				}else{ready = 1;}
			}
		}






		/**
		 *  An ascending order
		 **/
		void sort_vect(double* firstVec, double* secondVec, double* thirdVec, int checkFirst, int checkSecond, int checkThird, POPULATIONN* population, int start, int end, int check)
		{
			int i;
			INDIVIDUAL *holder;
			double auxHolder;
			int half = start + (int)(end-start)/2;

			for(i = half; i > start; i--)
			{
				if(checkFirst){sub_sort_vect(firstVec, secondVec, thirdVec, population, i-1, end, check);}
				else if(checkSecond){sub_sort_vect(secondVec, firstVec, thirdVec, population, i-1, end, check);}
				else if(checkThird){sub_sort_vect(thirdVec, secondVec, firstVec, population, i-1, end, check);}
			}

			for(i = end-1; i > start; i--)
			{
				holder = population->ind[start];
				population->ind[start] = population->ind[i];
				population->ind[i] = holder;

				auxHolder = firstVec[start];
				firstVec[start] = firstVec[i];
				firstVec[i] = auxHolder;

				if(check)
				{
					//auxHolder = secondVec[start];
					//secondVec[start] = secondVec[i];
					//secondVec[i] = auxHolder;

					auxHolder = thirdVec[start];
					thirdVec[start] = thirdVec[i];
					thirdVec[i] = auxHolder;
				}

				if(checkFirst){sub_sort_vect(firstVec, secondVec, thirdVec, population, start, i-1, check);}
				else if(checkSecond){sub_sort_vect(secondVec, firstVec, thirdVec, population, start, i-1, check);}
				else if(checkThird){sub_sort_vect(thirdVec, secondVec, firstVec, population, start, i-1, check);}
			}
		}

		//sub_sort_vect(firstVec, secondVec, thirdVec, population, i-1, end, check);
		void sub_sort_vect(double* firstVec, double* secondVec, double* thirdVec, POPULATIONN* population, int root, int end, int check)
		{
			INDIVIDUAL *holder;
			double auxHolder;
			int ready = 0, son;
			int verify = root*2;

			if(check)
			{
				while(verify < end && !ready)
				{
					if(firstVec[verify] > firstVec[verify+1]){son = verify;}
					else{son = verify+1;}

					if(firstVec[root] < firstVec[son])
					{
						holder = population->ind[root];
						population->ind[root] = population->ind[son];
						population->ind[root] = holder;

						auxHolder = firstVec[root];
						firstVec[root] = firstVec[son];
						firstVec[son] = auxHolder;


						//auxHolder = secondVec[root];
						//secondVec[root] = secondVec[son];
						//secondVec[son] = auxHolder;

						auxHolder = thirdVec[root];
						thirdVec[root] = thirdVec[son];
						thirdVec[son] = auxHolder;	


						root = son;
						verify = root*2;
					}else{ready = 1;}
				}
			}else
			{
				while(verify < end && !ready)
				{
					if(firstVec[verify] > firstVec[verify+1]){son = verify;}
					else{son = verify+1;}

					if(firstVec[root] < firstVec[son])
					{
						holder = population->ind[root];
						population->ind[root] = population->ind[son];
						population->ind[root] = holder;

						auxHolder = firstVec[root];
						firstVec[root] = firstVec[son];
						firstVec[son] = auxHolder;

						root = son;
						verify = root*2;
					}else{ready = 1;}
				}
			}
		}


		INDIVIDUAL *tournament(INDIVIDUAL *one, INDIVIDUAL *two)
		{
			int flag = one->check_dominance(two->get_obj1(), two->get_obj2(), two->get_obj3()); 
			float lastResort = ((float)( (rand()%10) + 1)/10);

			if(flag == 1){return two;}
			else if(flag == -1){return one;}
			else if(one->get_crowd_dist() > two->get_crowd_dist()){return one;}
			else if(two->get_crowd_dist() > one->get_crowd_dist()){return two;}
			else if(lastResort <= 0.5){return one;}
			else{return two;}
		}

};


#endif