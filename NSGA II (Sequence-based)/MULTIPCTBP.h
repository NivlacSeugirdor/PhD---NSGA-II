#ifndef _MULTIPCTBP_
#define _MULTIPCTBP_

#include <sstream>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilomodel.h>
#include <iostream>
#include <string>
#include <cstring>
#include <string.h>
#include <ctype.h>
#include <vector>
#include <map>
#include <math.h> 
#include <random>
#include <ctime>
#include <cfloat>
#include <algorithm>

#include "Chrono.h"
#include "AIRPORT.h"
#include "DESTINATION.h"
#include "FLIGHT.h"
#include "FLIGHTTS.h"
#include "PERIOD.h"
#include "TIMELINE.h"
#include "VARIABLE_Xijt.h"

ILOSTLBEGIN
using namespace std;

/**
*	Model paper of Mariá, Leandro and Jacques (JORS - 04/12/2016)
*	Model paper Mariá and Horácio
*   Implemented on Java by Katia Nakamura and to C++ by Calvin Rodrigues
*	New features 10/03/2021 by Calvin Rodrigues 
*					- A new FO considering penalties for does not visit some of presented locations;
*					- How to choose the day to arrive on an location.
* 	New features 29/09/2022
* 					- Now there are 3 FOs:
* 						. 1st - Min of costs
* 						. 2nd - Max number of visited cities
* 						. 3rd - Min the use of the slack days
* 					
* 					- Change on cities index control. 
* 						- before - N cities + 1
* 						- Now - N cities
* 						- Reduction on more less 25% in execution time in some tests.
*/

class MULTIPCTBP
{	
	public:
	
		MULTIPCTBP()
		{
			obj1 = 999999;
			obj2 = 999999; 
			obj3  = 0; 
			best = 0;
		}


		//Just ignore gamma

		void call_model(vector<DESTINATION*> destinations, vector<FLIGHT*> flights, map<string, int> flightData, FLIGHTTS tripSearch, 
															TIMELINE tl, AIRPORT depot, float gamma, float norm, int* toSolve, float *oldOb1, float *oldOb2, float *oldOb3,
																																						float alphaCost, float alphaPrize)
		{ 
			//Class
			vector<VARIABLE_Xijt*> varVec; 
			vector<FLIGHT*> fromItoJ; 
			VARIABLE_Xijt* varXijt;
			FLIGHT* fly;

			best = 0;	

			normal = norm;

		    //cout.unsetf(ios::scientific);
			
			//Timer
			ChronoReal cl;
		   	double solvF = 0, readCons;

			//commom variables
			int slackPeriods, stayDays;
			int dSize, tlSize, travelTime, varVSize;
			char buffer[15];
			int pmin;
			int i, j, t, f;// i, j, t for x_ijt, f to flights
			int M = 100000, gijt;
			int s1, s2;
			string cpy, cpy2;
			int min_xfiNt;
			int ***xijtTq; // <<<< quantity of ticket in each ijt

			cl.start();
			//Getting data
			
			//Updating values for quantity of destinations, size of the time line and the measuring of how many time of slack period
			dSize = destinations.size();
			tlSize = tl.get_period_by_index().size();
			slackPeriods = tripSearch.get_period_value()*tripSearch.get_slack_days();
			
			fromItoJ.reserve(dSize);

			xijtTq = (int***)malloc(sizeof(int**)*(dSize+1));

			for(i = 0; i < dSize; i++)
			{
				xijtTq[i] = (int**)malloc(sizeof(int*)*(dSize+1));

				for(j = 0; j < dSize; j++) // + k
				{
					xijtTq[i][j] = (int*)malloc(sizeof(int)*(tlSize));

					for(t = 0; t < tlSize; t++)
						 xijtTq[i][j][t] = 0;
				}
			}

			
			for(i = 0; i < dSize; i++) // check
			{

				for(j = 0; j < dSize; j++) // check
				{

		            for(t = 0; t < tlSize; t++) // check
					{    

						fromItoJ.clear(); // check

		            	if(i!=j)
						{

							if(tl.get_period_by_index()[t+1] != NULL) // check
							{

								if(!tl.get_period_by_index().at(t+1)->get_flights_by_from().empty()) // check
			 					{

			 						cpy = destinations[i]->get_destination()->get_code(); 
			 						
			 						if(tl.get_period_by_index()[t+1]->get_flights_by_from().count(cpy) > 0) // check 
			 						{


			 							s1 = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy).size();
			 							
			 							for(f = 0; f < s1; f++) // check
			 							{
			 								fly = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy)[f];

			 								cpy2 = fly->get_to_where()->get_code();

			 								if( compare_string(cpy2, destinations[j]->get_destination()->get_code()) == 0) // check
				 							{
				 								fromItoJ.push_back(fly);
				 							}
			 							}
			 						}
								}
							}

							xijtTq[i][j][t] = fromItoJ.size();
						}
					}
				}
			}

			
			//cplex
			IloEnv env;
			IloModel mod(env);
			IloTimer timer(env);
			IloIntVar di[dSize+1];
			IloIntVar *st; // New -- refers to the slack time used in each stop
			IloIntVar ****xijtk; //****xijt
			IloIntVar *yi;

			IloExpr burn(env), _xijt_i(env), _xijt_j(env), sum_txijt(env);
			IloExpr sum_pyi(env), sum_yi(env); 


			

			

			//Initializing some variables of cplex
			try
			{
				xijtk = (IloIntVar****)malloc(sizeof(IloIntVar***)*(dSize+1));
				st = (IloIntVar*)malloc(sizeof(IloIntVar)*(dSize+1)); //new
				yi = (IloIntVar*)malloc(sizeof(IloIntVar)*(dSize+1));
				
				for(i = 0; i < dSize; i++)
				{
					xijtk[i] = (IloIntVar***)malloc(sizeof(IloIntVar**)*(dSize+1));
					
					sprintf(buffer, "y_%d", i);
					yi[i] = IloIntVar(env, 0, 1, buffer);
					sprintf(buffer, "st_%d", i);
					st[i] = IloIntVar(env, 0, tlSize, buffer);

					for(j = 0; j < dSize; j++) // + k
					{
						xijtk[i][j] = (IloIntVar**)malloc(sizeof(IloIntVar*)*(tlSize)); 

						for(t = 0; t < tlSize; t++)
						{
							xijtk[i][j][t] = (IloIntVar*)malloc(sizeof(IloIntVar)*(xijtTq[i][j][t]));
							

							for(int k = 0; k < xijtTq[i][j][t]; k++)
							{
								sprintf(buffer, "x_%d_%d_%d_%d", i, j, t, k);
								xijtk[i][j][t][k] = IloIntVar(env, 0, 1, buffer);  
								//cout << i << " " << j << " " << t << "\n";
								//scanf("%*c");
							}
							     
						}
					}
				}


				/*.............................................*/
				IloExpr ob(env), ob2(env), ob3(env);
				
				ob = 0*burn;
				ob2 = 0*st[0];
				ob3 = 0*yi[0];
				/*.............................................*/
				

				for(i = 0; i < dSize; i++)
				{
					sprintf(buffer, "t_%d", i+1);
					di[i] = IloIntVar(env, 0, tlSize);
					di[i].setName(buffer);
				}
				
		        // =====================================================================================================================
				/**
				*	Flow-based model
				*	FO construction
				*
				*	FO = Sum_{i\in N}Sum_{j\in N}Sum_{t\in T}Sum_{k\in P_{ijt}} c_{ij}^{tk}x_{ij}^{tk}
				*
				*	For this expression the three loops ahead will represent:
				*		- 1st loop (i):	the airports of departure; 
				*		- 2nd loop (j):	the airports of arrival;
				*		- 3rd loop (t): the period.
				*
				*   So, if:
				*		- i == j, there is nothing to do. T
				*		- i != j, add the cheapest from i to j on period t in expression;
				*		- j == |N|, check if there is a passage from i to the origin on period t.
				*				. If yes, add the cheapest passage in expression;
				*				. Otherwise, do nothing.
				*  ->rewrite
				*  An observation about this construction of FO:
				*	- The selection of only the cheapests passages reduce the number of possible 
				*	combinations to be verified. But this selection also allow the lost of possible better
				*	alternatives.
				*/

				//cout << "Obj: ";

			    for(i = 0; i < dSize; i++) // check
				{			
					//ob2 = ob2 + st[i]*(1/slackPeriods); // new  - normalized
					
					for(j = 0; j < dSize; j++) 
					{
		                
		                for(t = 0; t < tlSize; t++) 
						{    
							fromItoJ.clear(); 

							if(i!=j)
							{
								if(tl.get_period_by_index()[t+1] != NULL) 
								{

									if(!tl.get_period_by_index().at(t+1)->get_flights_by_from().empty()) 
				 					{

				 						cpy = destinations[i]->get_destination()->get_code(); 
				 						
				 						if(tl.get_period_by_index()[t+1]->get_flights_by_from().count(cpy) > 0) 
				 						{


				 							s1 = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy).size();
				 							
				 							for(f = 0; f < s1; f++) 
				 							{
				 								fly = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy)[f];

				 								cpy2 = fly->get_to_where()->get_code();

				 								if( compare_string(cpy2, destinations[j]->get_destination()->get_code()) == 0) 
					 							{
					 								fromItoJ.push_back(fly);
					 							}
				 							}
				 						}
									}
								}

								//Sort -> get the smallest value without the use of merge
								fromItoJ = merge_sort(fromItoJ, fromItoJ.size());
								
								// Introduce the available flights for i to j in period t
								for(int k = 0; k < xijtTq[i][j][t]; k++)
								{
									travelTime = travel_periods(tripSearch, fromItoJ, k);

									stayDays = destinations[j]->get_daysstay();

									gijt = travelTime + (stayDays*tripSearch.get_period_value()); // check
									
									ob = ob + (fromItoJ[k]->get_price())*xijtk[i][j][t][k];// + K[j]*lijtk[i][j][t][k];
									
									//cout << " + "<< fromItoJ[k]->get_price() << "*x[" << i << "][" << j << "][" << t << "][" << k << "]";
									varXijt = new VARIABLE_Xijt(i, j, t, k); // check

									varXijt->set_price(fromItoJ[k]->get_price()); // check
									varXijt->set_vars(fromItoJ[k]->get_from_where()->get_name(), fromItoJ[k]->get_to_where()->get_name(), fromItoJ[k]->get_departure_time(), fromItoJ[k]->get_arrive_time());
									varVec.push_back(varXijt); // check
								}
							}
						}
					}

					ob3 = ob3 + yi[i]*destinations[i]->get_prize();
				}
					
				//cout<< "\n";

				//ob3 = (ob3)/dSize; // new  - normalized
				
		        varVSize = varVec.size();


		        IloObjective fo;

//				fo = IloObjective(env, IloMinimize(env, ob + ob2 - ob3)); // Objective function about the costs

				fo = IloObjective(env, IloMinimize(env, (1)*ob - (1)*ob3)); // Objective function about the costs

		        
				mod.add(fo); // add it to the model
				
				/***  Restrictions ***/



				/*
					================================================================================================================================================
					================================================================================================================================================
				*/

				//1st change - its only for the purpose of reselve a already existing sequence
				for(i = 0; i <= dSize; i++) //ok
				{
					if(toSolve[i] >= 0){mod.add(yi[toSolve[i]] == 1);}
					else{mod.add(yi[-toSolve[i]] == 0);}
				}
			

				/*
					================================================================================================================================================
					================================================================================================================================================
				*/

				/**
				*	2nd: Sum_{t} Sum_{j} Sum_{k} x_ijt <= 1, for all i
				*
				*	These restrictions require a departure on each airport only once.
			    */ 
				//2nd change - mod construct
				for(i = 0; i < dSize && toSolve[i+1] >= 0; i++) //-> airports to depart
				{
					_xijt_i = IloExpr(env);	

					for(j = 0; j < varVSize; j++) // verified flights
					{	
						if( (varVec[j]->get_j() == toSolve[i+1]) && (varVec[j]->get_i() == toSolve[i]) )
						{
							_xijt_i = _xijt_i + xijtk[varVec[j]->get_i()][varVec[j]->get_j()][varVec[j]->get_t()][varVec[j]->get_k()]; 
						}
					}

					sprintf(buffer, "rest_1_i_%d", i);
					_xijt_i.setName(buffer);

			    	mod.add(_xijt_i == yi[toSolve[i]]);	
				}


				for(i = 0; i < dSize && toSolve[i+1] >= 0; i++) //-> airports to depart
				{
					_xijt_i = IloExpr(env);	

					for(j = 0; j < varVSize; j++) // verified flights
					{	
						if( (varVec[j]->get_j() != toSolve[i+1]) && (varVec[j]->get_i() == toSolve[i]) )
						{
							_xijt_i = _xijt_i + xijtk[varVec[j]->get_i()][varVec[j]->get_j()][varVec[j]->get_t()][varVec[j]->get_k()]; 
						}
					}

					sprintf(buffer, "rest_1_2_i_%d", i);
					_xijt_i.setName(buffer);

			    	mod.add(_xijt_i == 0);	
				}

				/*
					================================================================================================================================================
					================================================================================================================================================
				*/

				/**
			    *	3rd: Sum_{t} Sum_{i} Sum_{k} x_ijtk <= 1, for all j
				*	
				*	These restrictions require an arrival from each airport only once.
				*/
				//3rd change  - mod construct
				for(j = 1; j < dSize && toSolve[j] >= 0; j++) // -> airports to arrive
				{
					_xijt_j = IloExpr(env);

					for(i = 0; i < varVSize; i++) // -> verified flights
					{
						if(  (varVec[i]->get_j() == toSolve[j]) && varVec[i]->get_i() == toSolve[j-1])
						{
							_xijt_j = _xijt_j + xijtk[varVec[i]->get_i()][varVec[i]->get_j()][varVec[i]->get_t()][varVec[i]->get_k()];
						}
					}

					sprintf(buffer, "rest_2_j_%d", j);
					_xijt_j.setName(buffer);

			    	mod.add(_xijt_j == yi[toSolve[j]]);
				}


				for(j = 1; j < dSize && toSolve[j] >= 0; j++) // -> airports to arrive
				{
					_xijt_j = IloExpr(env);

					for(i = 0; i < varVSize; i++) // -> verified flights
					{
						if(  (varVec[i]->get_j() == toSolve[j]) && varVec[i]->get_i() != toSolve[j-1])
						{
							_xijt_j = _xijt_j + xijtk[varVec[i]->get_i()][varVec[i]->get_j()][varVec[i]->get_t()][varVec[i]->get_k()];
						}
					}

					sprintf(buffer, "rest_2_j_%d", j);
					_xijt_j.setName(buffer);

			    	mod.add(_xijt_j == 0);
				}
				/*
					================================================================================================================================================
					================================================================================================================================================
				*/
			


			    /**
			    *	4th: Sum_{i}Sum_{t} y_{it}p_{i} >= p_{min}
			    */

			    sum_pyi = IloExpr(env);

			    for(i = 0; i < dSize; i++)
			    {
			    	sum_pyi = sum_pyi + yi[i]*destinations[i]->get_prize();
			    }

			    mod.add(sum_pyi >= 10);//mod.add(sum_pyi >= pmin);


				/*
					================================================================================================================================================
					================================================================================================================================================
				*/


			    /**
				*	5th: Sum_{j} Sum_{t} (Sum_{k} t*x_{ijtk} <= d_{i}, i\in N
				*		Estabilish a relation between x[i][j][t] and the integer variable d_{i} of departure time by the sum of t periods.
			    */

			    for(i = 0; i < dSize; i++) // check
			    {
			    	sum_txijt = IloExpr(env);
			    	for(j = 0; j < varVSize; j++) // check 
			    	{
			    		if(varVec[j]->get_i() == i && varVec[j]->get_j() != i ) // check
			    		{
				    		sum_txijt = sum_txijt + (varVec[j]->get_t()+1)*xijtk[i][varVec[j]->get_j()][varVec[j]->get_t()][varVec[j]->get_k()]; // check
				    	}
			    	}
			    	sprintf(buffer, "rest_5_i_%d", i);
					sum_txijt.setName(buffer);

			    	mod.add(sum_txijt == di[i]);
			    }
			    

				/*
					================================================================================================================================================
					================================================================================================================================================
				*/


			    /**
				*	6th: d_i + M*x^{tk}_{ij} - d_j <= M - (g^{tk}_{ij}+1) --> [d_{jk} >= d_{i} - (1-x^{tk}_{ij})*M], i,j\in N, t\in T & k in P_{ijt}
				*		Ensure that the next arrival period must be ever greater or equal to the departure plus travel time of predecessor. 		
			    */


			    for(i = 0; i < dSize; i++) // check
			    {
			    	
			    	for(j = 1; j < dSize; j++) // check
			    	{
			    		
			    		for(t = 0; t < tlSize; t++) // check
			    		{
			    			fromItoJ.clear();
			    			
			    			if(i!=j) // check
				    		{
				    			if(tl.get_period_by_index()[t+1] != NULL) // check
				    			{
				    				if(!tl.get_period_by_index().at(t+1)->get_flights_by_from().empty()) // check
				    				{
				    					cpy = destinations[i]->get_destination()->get_code(); 
									
										if(tl.get_period_by_index()[t+1]->get_flights_by_from().count(cpy) > 0) // check
										{
											s1 = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy).size();

											for(f = 0; f < s1; f++) // check
											{
												fly = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy)[f];

												cpy2 = fly->get_to_where()->get_code();

												if( compare_string(cpy2, destinations[j]->get_destination()->get_code()) == 0) // check
												{
													fromItoJ.push_back(fly);
												}
											}
										}   
				    				} 		
								}

								//Sort -> get the smallest value without the use of merge
								fromItoJ = merge_sort(fromItoJ, fromItoJ.size());

								//ATRIBUTE
								for(int k = 0; k < xijtTq[i][j][t]; k++)
								{
									travelTime = travel_periods(tripSearch, fromItoJ, k);

									stayDays = destinations[j]->get_daysstay();

									gijt = travelTime + (stayDays*tripSearch.get_period_value()); // check
									
									mod.add(di[j] >= di[i] + gijt - (1-xijtk[i][j][t][k])*M);
								}
				    		}

			    		}
			    	}
			    }
				

				/*
					================================================================================================================================================
					================================================================================================================================================
				*/


			    /**
				*	7th: d[j] - d[i] + tlSize*x^{tk}_{ij} <= tlSize + g^{tk}_{ij} + slackPeriods
				*		
				*
			    */
			    

			
			    for(i = 0; i < dSize; i++) // check
			    {
			    	
			    	for(j = 0; j < dSize; j++) // check
			    	{

			    		for(t = 0; t < tlSize; t++) // check
			    		{
			    			fromItoJ.clear();

			    			if(i!=j)
					    	{
					    		if(tl.get_period_by_index()[t+1] != NULL) // check
			    				{
		    						if(!tl.get_period_by_index().at(t+1)->get_flights_by_from().empty()) // check
		    						{
										cpy = destinations[i]->get_destination()->get_code();

					    				if(tl.get_period_by_index()[t+1]->get_flights_by_from().count(cpy) > 0) // check
					    				{
					    					s1 = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy).size();

				 							for(f = 0; f < s1; f++) // check
				 							{
				 								fly = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy)[f];

					 							cpy2 = fly->get_to_where()->get_code();

					 							if( compare_string(cpy2, destinations[j]->get_destination()->get_code()) == 0) // check
					 							{
					 								fromItoJ.push_back(fly);
					 							}
				 							}
					    				}
		    						}
		    					}

					    		//Sort -> get the smallest value without the use of merge
			    				fromItoJ = merge_sort(fromItoJ, fromItoJ.size());

			    				for(int k = 0; k < xijtTq[i][j][t]; k++) // check
					    		{
					    			travelTime = travel_periods(tripSearch, fromItoJ, k);

					    			stayDays = destinations[j]->get_daysstay();

					    			gijt = travelTime + (stayDays*tripSearch.get_period_value());
					    			
					    			mod.add(di[j] <= di[i] + gijt + slackPeriods + (1-xijtk[i][j][t][k])*tlSize);
					    		}
					    	}
					    }

				    }

			    }
		 			

				/*
					================================================================================================================================================
					================================================================================================================================================
				*/

			    /**
				*	8th:  d_{i} \leq (T - \bar{g}_{i}), i\in \mathcal{N}
			    */

			    min_xfiNt = tlSize;

			    for(i = 1; i < dSize; i++) // check
			    {

			    	for(t = 0; t < tlSize; t++) // check
			    	{

			    		fromItoJ.clear();
			    		if(tl.get_period_by_index()[t+1] != NULL) // check
			    		{	
			    			if(!tl.get_period_by_index().at(t+1)->get_flights_by_from().empty()) // check
			    			{
					    	
				    			cpy = destinations[i]->get_destination()->get_code();

				    			if(tl.get_period_by_index()[t+1]->get_flights_by_from().count(cpy) > 0) // check
								{
									s1 = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy).size();

									for(f = 0; f < s1; f++) // check
									{
										fly = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy)[f];

										cpy2 = fly->get_to_where()->get_code();

										if( compare_string(cpy2, destinations[0]->get_destination()->get_code()) == 0) // check
										{
											fromItoJ.push_back(fly);
										}
									}
								}
			    			}
				    		
				    		//Sort -> get the smallest value without the use of merge
				    		fromItoJ = merge_sort(fromItoJ, fromItoJ.size());
			    			//atribute
				    		/*if(fromItoJ.size() > 0) // check
				    		{
				    			
				    			travelTime = travel_periods(tripSearch, fromItoJ, k);
				    			gijt = travelTime;
				    			
				    			if(gijt < min_xfiNt){min_xfiNt = gijt;} // check
				    		}*/
				    		for(int k = fromItoJ.size()-1; k >= 0 ; k--) // check
				    		{
				    			
				    			travelTime = travel_periods(tripSearch, fromItoJ, k);
				    			gijt = travelTime;
				    			
				    			if(gijt < min_xfiNt)
				    			{
				    				min_xfiNt = gijt;
				    			} // check
				    		}
			    		}


			    	}

			    	

			    	if(min_xfiNt < tlSize) // check
			    	{ 
			    		mod.add(di[i] <= (tlSize - min_xfiNt) ); // check
			    		min_xfiNt = tlSize;
			    	}else
			    	{			    			
			    		mod.add(di[i] <= tlSize);
			    	}
			    	

			    }

			    /*
					================================================================================================================================================
					================================================================================================================================================
				*/
			    //NEW
			    /**
				*	9th:  st_i_j^t <= slack*y[i], i,j\in \mathcal{N}, t\in \mathcal{T} 
			    */
			    /*
			    for(i = 0; i < dSize; i++) 
			    {
			    	mod.add(st[i] <= slackPeriods*yi[i]);//sum_txijt);
			    }
			    */
			    /*
					================================================================================================================================================
					================================================================================================================================================
				*/
			    //NEW
			    /**
				*	10th:  st_i_j^t >= d_j - d_i - g_i_j^t - (1-x_ijtk)*T, i\in \mathcal{N}
			    */
			    /**
			    for(i = 0; i < dSize; i++) // check
			    {
			    	
			    	for(j = 1; j < dSize; j++) // check j >=1 -> change
			    	{

			    		for(t = 0; t < tlSize; t++) // check
			    		{
			    			fromItoJ.clear();

			    			if(i!=j)
					    	{
					    		if(tl.get_period_by_index()[t+1] != NULL) // check
			    				{
		    						if(!tl.get_period_by_index().at(t+1)->get_flights_by_from().empty()) // check
		    						{
										cpy = destinations[i]->get_destination()->get_code();

					    				if(tl.get_period_by_index()[t+1]->get_flights_by_from().count(cpy) > 0) // check
					    				{
					    					s1 = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy).size();

				 							for(f = 0; f < s1; f++) // check
				 							{
				 								fly = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy)[f];

					 							cpy2 = fly->get_to_where()->get_code();

					 							if( compare_string(cpy2, destinations[j]->get_destination()->get_code()) == 0) // check
					 							{
					 								fromItoJ.push_back(fly);
					 							}
				 							}
					    				}
		    						}
		    					}

					    		//Sort -> get the smallest value without the use of merge
			    				fromItoJ = merge_sort(fromItoJ, fromItoJ.size());

			    				for(int k = 0; k < xijtTq[i][j][t]; k++) // check
					    		{
					    			travelTime = travel_periods(tripSearch, fromItoJ, k);

					    			stayDays = destinations[j]->get_daysstay();

					    			gijt = travelTime + (stayDays*tripSearch.get_period_value());
					    			
					    			mod.add(st[j] >= di[j] - di[i] - gijt - (1-xijtk[i][j][t][k])*tlSize); 
					    		}
					    	}
					    }

				    }

			    }
			    */
				/*
					================================================================================================================================================
					================================================================================================================================================
				*/


			    /**
				*	11th: d_{0} <= s
			    */
			    
			  	mod.add(di[0] <= slackPeriods); // check


				/*
					================================================================================================================================================
					================================================================================================================================================
				*/


			  	/**
			  	*	12th:
			  	*		(1-L[j])x[i][j][t]P[j][1] <= d[i] + f[i][t][j];
			  	*/
			    

			  	for(i = 0; i < dSize; i++) // check
			    {
			    	
			    	for(j = 0; j < dSize; j++) // check
			    	{
			    		
			    		for(t = 0; t < tlSize; t++) // check
			    		{
			    			fromItoJ.clear();
			    			if(i!=j) // check
				    		{
				    			if(tl.get_period_by_index()[t+1] != NULL) // check
				    			{
				    				if(!tl.get_period_by_index().at(t+1)->get_flights_by_from().empty()) // check
				    				{
				    					cpy = destinations[i]->get_destination()->get_code(); 
									
										if(tl.get_period_by_index()[t+1]->get_flights_by_from().count(cpy) > 0) // check
										{
											s1 = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy).size();

											for(f = 0; f < s1; f++) // check
											{
												fly = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy)[f];

												cpy2 = fly->get_to_where()->get_code();

												if( compare_string(cpy2, destinations[j]->get_destination()->get_code()) == 0) // check
												{
													fromItoJ.push_back(fly);
												}
											}
										}   
				    				} 		
								}

								//Sort -> get the smallest value without the use of merge
								fromItoJ = merge_sort(fromItoJ, fromItoJ.size());

								//ATRIBUTE
								for(int k = 0; k < xijtTq[i][j][t]; k++)
								{
								    travelTime = travel_periods(tripSearch, fromItoJ, k);

									gijt = travelTime; // check
									
									mod.add( xijtk[i][j][t][k]*destinations[j]->get_min_arr() <= di[i] + gijt);
								}
				    		}

			    		}
			    	}
			    }


				/*
					================================================================================================================================================
					================================================================================================================================================
				*/


			  	/**
			  	*	13th:
			  	*		d[i] + f[i][t][j] <= x[i][j][t]P[j][2] + L[j]M;
			  	*/
			  	for(i = 0; i < dSize; i++) // check
			    {
			    	
			    	for(j = 1; j < dSize; j++) // check
			    	{
			    		
			    		for(t = 0; t < tlSize; t++) // check
			    		{
			    			fromItoJ.clear();
			    			
			    			if(i!=j) // check
				    		{
				    		    if(tl.get_period_by_index()[t+1] != NULL) // check
				    			{
				    				if(!tl.get_period_by_index().at(t+1)->get_flights_by_from().empty()) // check
				    				{
				    					cpy = destinations[i]->get_destination()->get_code(); 
									
										if(tl.get_period_by_index()[t+1]->get_flights_by_from().count(cpy) > 0) // check
										{
											s1 = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy).size();

											for(f = 0; f < s1; f++) // check
											{
												fly = tl.get_period_by_index()[t+1]->get_flights_by_from().at(cpy)[f];

												cpy2 = fly->get_to_where()->get_code();

												if( compare_string(cpy2, destinations[j]->get_destination()->get_code()) == 0) // check
												{
													fromItoJ.push_back(fly);
												}
											}
										}   
				    				} 		
								}

								//Sort -> get the smallest value without the use of merge
								fromItoJ = merge_sort(fromItoJ, fromItoJ.size());

								//ATRIBUTE
								for(int k = 0; k < xijtTq[i][j][t]; k++)
								{
								    mod.add( di[j] <= xijtk[i][j][t][k]*destinations[j]->get_max_dep() + (1-xijtk[i][j][t][k])*tlSize);
								}
				    		}

			    		}
			    	}
			    }

			    

			    clock_t start_it, end_it, start, end;
			    IloCplex cplex(mod);
			    

			    //Use only one Thread
			    cplex.setParam(IloCplex::Threads, 1);
			    float gap = random_value();
			    cplex.setParam(IloCplex::EpGap, gap);
			    cplex.setOut(env.getNullStream()); 

			    cl.stop();
			    readCons = cl.getTime();

			    cl.reset();
			    //start = clock();
			    


			 	for(i = 0; i < 1; i++)
			 	{
			 		//start_it = clock();	
			 		cl.start();

					cplex.setParam(IloCplex::TiLim, 3600.0); 
					j = cplex.solve();
					//end_it = clock();
					cl.stop();

					//print_iteractions(argv[1], j, 0, readCons+cl.getTime(), cplex, i);
					//solvF += (readCons + cl.getTime());
					cl.reset();
			 	}

			 	//end = clock();
			 	if(cplex.isPrimalFeasible())
			 	{
			 		feasibleCheck = 1;
			 		
			 		update_values(destinations, flights, cplex, oldOb1, oldOb2, oldOb3, toSolve, ob, ob3);
			 		//char bf[] = "cd";
			 		//print_best_value(bf, destinations, cplex, xijtk, varVec, di, yi, ob, ob3);
			 	}


			 	free(st);
			 	free(yi);


			 	for(i = 0; i < dSize; i++)
			 	{
			 		for(j = 0; j < dSize; j++)
			 		{
			 			free(xijtTq[i][j]);

			 			for(t = 0; t < tlSize; t++)
			 				free(xijtk[i][j][t]);

			 			free(xijtk[i][j]);
			 		}
			 		free(xijtTq[i]);
			 		free(xijtk[i]);
			 	}

	 			free(xijtTq);
		 		free(xijtk);


		 		for(i = 0; i < varVec.size(); i++)
		 			delete varVec[i];

			 	varVec.clear();
				fromItoJ.clear();	

				env.end();
				//mod.end();
				//cplex.end();	

			}catch(IloException& ex) {
		     		cerr << "\nError: " << ex << endl;
		   	}catch (string er){
		      		cerr << "\nError \n" << er << "\n" << endl;
		   	}

		}


		int check(){return feasibleCheck;}

		float get_obj1(){return obj1;}
		float get_obj2(){return obj2;}
		float get_obj3(){return obj3;} 		

		void get_gene(int newBest[], int qtCities)
		{
			for(int i = 0; i <= qtCities; i++)
			{
				newBest[i] = newGene[i];
			}
			
			free(newGene);
		}

		/**
		 * if its true = 1
		 * else = 0 	
		 * */
		int im_the_best(){return best;}

	private:
		int feasibleCheck;
		int best;
		int *newGene;
		float obj1, obj2, obj3;
		float normal;

		int compare_string(string s1, string s2)
		{
			int check = -1;

			if(s1.size() == s2.size())
			{
				for(int i = 0; i < s2.size(); i++)
				{
					if(tolower(s1[i]) == tolower(s2[i])){check = 0;}
					else
					{
						check = -1;
						break;
					}
				}
			}

			return check;
		}

		static void exportLP(char arq[], IloCplex cplex) 
		{
			char fn[80];
			strcpy(fn, arq);
			strcat(fn, "_model_.lp");

			cplex.exportModel(fn);
		} 

		vector<FLIGHT*> merge(vector<FLIGHT*> f, int l, int m, int r)
		{
			int i, j, k;
			int n1, n2;
			n1 = m - l + 1;
			n2 = r-m;


			vector<FLIGHT*> L, R;



			for( i = 0; i < n1; i++)
				L.push_back(f[l+i]);


			for(i = 0; i < n2; i++)
				R.push_back(f[m+1+i]);

			i = 0;
			j = 0;
			k = l;

			while(i < n1 && j < n2)
			{
				if(L[i]->get_price() <= R[j]->get_price())
				{
					f[k] = L[i];
					i++;
				}else
				{
					f[k] = R[j];
					j++;
				}
				k++;
			}

			while(i < n1)
			{
				f[k] = L[i];
				i++;
				k++;
			}

			while(j < n2)
			{
				f[k] = R[j];
				j++;
				k++;
			}

			return f;
		}

		vector<FLIGHT*>merge_sort(vector<FLIGHT*> f, int n)
		{
			int currSize, leftStart;
			int mid, rightEnd;

			

			for(currSize = 1; currSize < n; currSize = 2*currSize)
			{
				for(leftStart = 0; leftStart < n-1; leftStart+= 2*currSize)
				{
					mid = min(leftStart+currSize-1, n-1);
					rightEnd = min(leftStart + 2*currSize -1, n-1);

					f = merge(f, leftStart, mid, rightEnd);
				}
			}

			return f;
		}

		//Must to adapt it to save the best result. 
		static void print_best_value(char arc[], vector<DESTINATION*> destinations, IloCplex cplex, IloIntVar ****xijtk, vector<VARIABLE_Xijt*> varVec, IloIntVar di[], IloIntVar yi[], IloExpr ob, IloExpr ob2) 
		{
			FILE *fp;
			FILE *fp2;

			char fn[50], fn2[50];//name of the archive to be read (fn = file's name)
			strcpy(fn, arc);
			strcat(fn, "_outputs.txt");

			strcpy(fn2, arc);
			strcat(fn2, ".ods");

			int size = varVec.size();
			int i;
			float sum=0;

			fp = fopen(fn, "a+");
			fp2 = fopen(fn2, "a+");

			fprintf(fp, "Time on seconds: ---\n");

			for(i = 0; i < size; i++)
			{
				if(cplex.getValue( xijtk[ varVec[i]->get_i() ][ varVec[i]->get_j() ][ varVec[i]->get_t() ][varVec[i]->get_k()]) > 0)
				{
					fprintf(fp, "x[%d][%d][%d][%d] - %f - Preço: %f - Dep %lf - Arr %lf\n", varVec[i]->get_i(), varVec[i]->get_j(), varVec[i]->get_t(), varVec[i]->get_k(), cplex.getValue(xijtk[varVec[i]->get_i()][varVec[i]->get_j()][varVec[i]->get_t()][varVec[i]->get_k()]), varVec[i]->get_price(), varVec[i]->get_dep(), varVec[i]->get_arr()); 	
				}
			}


			for(i = 0; i < destinations.size(); i++)
			{
				fprintf(fp, "di[%d]  - %d\n", i, (int)cplex.getValue(di[i]));
				//if((int)cplex.getValue(yi[i])==0){sum++;}
			}

			fprintf(fp, "FO total = %lf, Custo = %lf, Premio: %lf", (double)cplex.getObjValue(), cplex.getValue(ob), cplex.getValue(ob2));


			fclose(fp);
			fclose(fp2);
		} 

		// https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
		float random_value() {
			random_device rd;
			mt19937 gen(rd());
			uniform_real_distribution<float> dis(0, 100);
			  
			return ( ((int)dis(gen))*1.0)/100;
		}

		static int travel_periods(FLIGHTTS tripSearch, vector<FLIGHT*> fromItoJ, int k)
		{
			int arriveTime = (int)ceil((fromItoJ[k]->get_arrive_time() - tripSearch.get_departure_date())/14400.00); // 1 day equals to 6 periods of 14400 sec each
			int departureTime = (int)ceil((fromItoJ[k]->get_departure_time()  - tripSearch.get_departure_date())/14400.00);
			return (arriveTime - departureTime);
		}

		void update_values(vector<DESTINATION*> destinations, vector<FLIGHT*> flights, IloCplex cplex, float *oldOb1, float *oldOb2, float *oldOb3, int gene[], IloExpr ob, IloExpr ob3) 
		{
			int size_C = destinations.size();
			int i;
			float sum=0;


			obj1 = 0;
			obj2 = 0; 
			obj3  = 0; 

			obj1 = cplex.getValue(ob);;
			obj3  = cplex.getValue(ob3); 



			if(obj1 < *oldOb1 || obj3 > *oldOb3)
			{
				//if((obj1 <= *oldOb1) && (obj2 <= *oldOb2) && (obj3 >= *oldOb3))
				if((obj1 <= *oldOb1) && (obj3 >= *oldOb3))
				{
					best = 1;
					set_gene(size_C, gene);
				}
			}

			/*if(obj1 <= *oldOb1 && obj2 <= *oldOb2 && obj3 >= *oldOb3)
			{
				if( (obj1 < *oldOb1) || (obj2 < *oldOb2) || (obj3 > *oldOb3)  )
				{
					best = 1;
					set_gene(size_C, gene);
				}
			}*/
			
		} 

		void set_gene(int geneSize, int gene[])
		{
			newGene = (int*)malloc(sizeof(int)*geneSize);

			for(int i = 0; i < geneSize; i++)
				newGene[i] = gene[i];

		}

};

#endif