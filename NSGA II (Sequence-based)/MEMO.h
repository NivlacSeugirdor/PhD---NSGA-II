# ifndef _MEMO_H_
# define _MEMO_H_

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


#include "DESTINATION.h"
#include "AIRPORT.h"
#include "FLIGHT.h"
#include "FLIGHTTS.h"
#include "PERIOD.h"
#include "TIMELINE.h"
#include "VARIABLE_Xijt.h"
#include "MULTIPCTBP.h"

using namespace std;


class MEMO
{
	public:
		MEMO(char argv[])
		{
			read_archive(argv);
			float auxPrize;
			prizeMin = 999999;
			prizeTotal = 0;

			for(int i = 0; i < destinations.size(); i++)
			{
				auxPrize = destinations[i]->get_prize();

				if(auxPrize < prizeMin)
					prizeMin = auxPrize;

				prizeTotal = auxPrize;
			}

			//gammaVal = gammaVal*atof(argv[2]);
		}

		void evaluates_sol(int toSolve[], float *ob1, float *ob2, float *ob3, float alphaCost, float alphaPrize)
		{
			mod = new MULTIPCTBP();			

			mod->call_model(destinations, flights, flightData, tripSearch, tl, depot, gammaVal, norm, toSolve, ob1, ob2, ob3, alphaCost, alphaPrize);
		}

		int check_best(){return mod->im_the_best();}

		int check_sol(){return mod->check();}



		//get tickets
		void get_best_genes(int newBest[]){mod->get_gene(newBest, destinations.size());}

		//get val1 min obj1
		float get_firstObj(){return mod->get_obj1();}

		//get val2 min obj2
		float get_secondObj(){return mod->get_obj2();}
		
		//get val3 max obj3
		float get_thirdObj(){return mod->get_obj3();}

		float get_gammaVal(){return gammaVal;}
	
		int get_gene_size(){return destinations.size();}

		float get_min_prize(){return prizeMin;}
	
		float get_norm(){return norm;}

		int get_pmin(){return pmin;}

		void model_terminator(){delete mod;}

		float get_total_prize(){return prizeTotal;}

	private:
		AIRPORT depot;
		vector<DESTINATION*> destinations;

		map<string, int> flightData;
		vector<FLIGHT*> flights;
		
		float gammaVal;
		MULTIPCTBP *mod;

		float prizeMin;
		float prizeTotal;

		float norm;
		int pmin, qtdC;

		TIMELINE tl;
		FLIGHTTS tripSearch;		


		
		int check_city_number(char air[])
		{
			char aux[10];

			for(int i = 0; i < destinations.size(); i++)
			{
				strcpy(aux, destinations[i]->get_destination()->get_code().c_str() );
				if(strcmp(air, aux) == 0){return i;}
			}

			return 0;
		}


		int get_period(int val , int t0)
		{
			return (int)ceil( (val-t0)/14400.00 );
		}

		int get_period2(int val , int t0)
		{
			return (int)floor( (val-t0)/14400.00 );
		}

		void read_archive(char arc[])
		{
			FILE *fp;

			//--------------------------------------------------------------------
			//--------------------------------------------------------------------
			
			/**/	char fn[50];//name of the archive to be read (fn = file's name)
			/**/	DESTINATION *dest;
			/**/
			/**/	AIRPORT *air;
			/**/	FLIGHT *f;
			/**/
			/**/	PERIOD period;
			/**/	map<int, PERIOD*> example;
			/**/	
			/**/	int qtd, index;
			/**/	char dept[4], arr[4];
			/**/
			/**/	int slack, timeStay, destQt;
			/**/	int s1, s2;
			/**/
			/**/	int prize, pm;
			/**/	int p1, p2;
			/**/
			/**/	double price;
			/**/	int t0;
			/**/
			/**/	int arrival, departure;
			/**/	float ga, cmin = 99999, cmax = 0; 

			//--------------------------------------------------------------------
			//--------------------------------------------------------------------



			strcpy(fn, arc);
			strcat(fn, "conv2.txt");


			fp = fopen(fn, "r");
			
			if(fp == NULL)
				exit(1);

			fscanf(fp, "%d %d", &slack, &t0);// time of tolerance, initial time
			fscanf(fp, "%d", &pm); //Minimal prize
			fscanf(fp, "%s", dept); // name of initial airport
			fscanf(fp, "%d", &destQt); // quantity of airports to visit
			
			depot.set_code(dept); 
		    
			pmin = pm;

			//Get the names and the time to stay in each airport
			for(int i = 0; i < destQt; i++)
			{
				fscanf(fp, "%s", dept);
				fscanf(fp, "%d", &timeStay);
				fscanf(fp, "%d", &prize);

				air = new AIRPORT();
				

				fscanf(fp, "%d", &p1);
				fscanf(fp, "%d", &p2);

				p1 = get_period(p1, t0);
				p2 = get_period(p2, t0);

				

				air->set_code(dept);

				dest = new DESTINATION(air, timeStay, prize, p1, p2);
				flightData.insert(pair<string, int>(dept, timeStay)); 
				destinations.push_back(dest); 
			}

			//--------------------------------------------------------------------
			//--------------------------------------------------------------------

			dest = new DESTINATION(&depot, 0, 0, 0, 0);

			destinations.resize(destinations.size()); 

			destinations.push_back(new DESTINATION()); 
			
			qtdC = destinations.size();

			for(int i = (destinations.size()-1); i > 0; i--)
				destinations[i] = destinations[i-1];  

			destinations[0] = dest; 


			for(int i = 0; i < qtdC; i++)
				destinations[i]->set_number(i);


			//--------------------------------------------------------------------
			//--------------------------------------------------------------------
			tripSearch.set_slack_days(slack); 
			tripSearch.set_departure_date(t0); 
			tripSearch.set_origin(&depot); 
			tripSearch.set_dest(flightData); 
			tripSearch.set_destinations(destinations); 
			tripSearch.set_period_value(6); 

			fscanf(fp, "%d", &qtd);
			
			for(int i = 0; i < qtd; i++)
			{
				f = new FLIGHT();

				fscanf(fp, "%d", &index);
				fscanf(fp, "%s", dept);
				fscanf(fp, "%s", arr);
				fscanf(fp, "%d %d %lf", &departure, &arrival, &price);


				if(cmin > price)
					cmin = price;

				if(cmax < price)
					cmax = price;


				air = new AIRPORT();
				air->set_code(dept);
				air->set_number_code( check_city_number(dept) );

				f->set_from_where(air);
				
				air = new AIRPORT();
				air->set_code(arr);

				air->set_number_code( check_city_number(arr) );

				f->set_ticket(index);

				f->set_to_where(air);

				f->set_departure_time(departure);
				f->set_converted_departure_time( get_period2(departure, t0) );
				f->set_arrive_time(arrival);
				f->set_converted_arrive_time( get_period2(arrival, t0) );


				f->set_price(price);

				flights.push_back(f); 
			}


			flights.pop_back(); 

			
			norm = cmax;	

			tl.create_time_line(tripSearch.get_budget()+1); 

			example = tl.get_period_by_index(); 

			period.set_flights_by_period_1(flights, t0, example); 

			destinations[0]->set_max_dep(example.size()); 

			gammaVal = (cmin+cmax)/norm;


			fclose(fp);
		}
};



#endif



