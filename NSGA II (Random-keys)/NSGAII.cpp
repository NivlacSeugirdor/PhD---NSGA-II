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


#include "Chrono.h"
#include "MEMO.h"
#include "POPULATIONN.h"
#include "POPULATIONNTREATMENT.h"


using namespace std;

void final_registration(float bestOb1, float bestOb2, float bestOb3, int* bestGenes, int geneSize, char arc[], double timeExec, float **dists, int pop, float df[], float dl[], int stop[])
{
	char fn[50];

	float average[3];

	average[0] = 0;
	average[1] = 0;
	average[2] = 0;

	for(int  i = 1; i < pop-1; i++)
	{
		if(stop[0]>i)
			average[0] += dists[0][i];
		if(stop[1]>i)
			average[1] += dists[1][i];
		if(stop[2]>i)
			average[2] += dists[2][i]; 
	}

	average[0] = average[0]/(stop[0]);
	average[1] = average[1]/(stop[1]);
	average[2] = average[2]/(stop[2]);


	float sum[3];

	for(int  i = 1; i < pop-1; i++)
	{
		if(stop[0]>i)
			sum[0] += abs(dists[0][i] - average[0]);
		
		if(stop[1]>i)
			sum[1] += abs(dists[1][i] - average[1]);
		
		if(stop[2]>i)
			sum[2] += abs(dists[2][i] - average[2]); 
	}


	sprintf(fn, "final_records.txt");
	FILE *fp = fopen(fn, "a+");

	fprintf(fp, "%s & ", arc);

	//for(int i = 0; i < geneSize+1; i++)
		//fprintf(fp, "%d ", bestGenes[i]);

	float a = ((df[0] + dl[0] + sum[0])/(df[0] + dl[0] + (stop[0])*average[0]));
	float b = ((df[1] + dl[1] + sum[1])/(df[1] + dl[1] + (stop[1])*average[1]));
	float c = ((df[2] + dl[2] + sum[2])/(df[2] + dl[2] + (stop[2])*average[2]));

	fprintf(fp, "%lf ", timeExec);//, a, b,  c);

	
	//for(int i = 0; i < 3; i++)
	//	free(dists[i]);

	//free(dists);

	fprintf(fp, "\n");

	fclose(fp);
}

void frontier_registration(POPULATIONN *mainPop, int pop, int it, char arc[], float **dists, int step, float df[], float dl[], int geneSize, int stop[], float minPrize, float maxPrize)
{
	char fn[50];

	sprintf(fn, "%s_%d_frontier.txt", arc, (it+1) );

	FILE *fp = fopen(fn, "a+");

	fprintf(fp, "Frontier - Dom - Crowd Dist - Min Obj1 - Max Obj3 - Sequence\n");

	float obj1F = 999999, obj3F = maxPrize;
	float obj1L = 0, obj3L = minPrize;
	float check, closest = 999999;
	int *gen;

	for(int reg = 0; reg < pop; reg++)
	{
		fprintf(fp, "%d - %d - %.4f - %.4f - %.4f - ", mainPop->ind[reg]->get_rank(), mainPop->ind[reg]->get_dominators(), mainPop->ind[reg]->get_crowd_dist(), mainPop->ind[reg]->get_obj1(), mainPop->ind[reg]->get_obj3());
		
		gen = mainPop->ind[reg]->get_gene();

		for(int i = 0; i < geneSize; i++)
			fprintf(fp, "%d ", gen[i]);

		//fprintf(fp, " -  %f - %f\n", mainPop->ind[reg]->get_alpha_cost(), mainPop->ind[reg]->get_alpha_prize());
		fprintf(fp, "\n");

		if(obj1F > mainPop->ind[reg]->get_obj1()){obj1F = mainPop->ind[reg]->get_obj1();}
		if(obj1L < mainPop->ind[reg]->get_obj1()){obj1L = mainPop->ind[reg]->get_obj1();}


		if(reg > 0 && reg < pop-1 && mainPop->ind[reg]->get_rank() == 0)
		{
			for(int j = 0; mainPop->ind[reg]->get_rank() == 0; j++)
			{
				if(mainPop->ind[reg]->get_obj1()!=mainPop->ind[j]->get_obj1() && mainPop->ind[reg]->get_obj3() != mainPop->ind[j]->get_obj3())
				{
					check = (float)sqrt(pow(mainPop->ind[reg]->get_obj1() - mainPop->ind[j]->get_obj1(), 2) + pow(mainPop->ind[reg]->get_obj3() - mainPop->ind[j]->get_obj3(), 2));
				
					if(closest > check){closest = check;}
				}
			}


			dists[step][reg] = closest;
			//dists[step][reg] = mainPop->ind[reg]->get_crowd_dist(); 
			//cout << "dists[" << step << "][" << reg << "] = " << dists[step][reg] << "\n";
			closest = 999999;
		}else if(mainPop->ind[reg]->get_rank()==1){stop[step]=reg;}
	}

	df[step] = (float)sqrt(pow(mainPop->ind[0]->get_obj1() -  obj1F, 2) + pow(mainPop->ind[0]->get_obj3() - obj3F, 2) ); //extreme is "min cost" & "max cities"
	dl[step] = (float)sqrt(pow(mainPop->ind[stop[step]]->get_obj1() - obj1L, 2) + pow(mainPop->ind[stop[step]]->get_obj3() - obj3L, 2) ); //extreme is "max cost" & "min cities"

	//cout << df[step] << " - " << dl[step] << "\n===================================\n";

	fclose(fp);
}


void solve_pop(MEMO* core, int pop, POPULATIONN* population, float *bestOb1, float *bestOb2, float *bestOb3, int *bestGenes)
{
	int i;
	float auxOb1, auxOb2, auxOb3;

	for(i = 0; i < pop; i++)
	{
		if(population->ind[i]->get_obj1() == 999999 &&
					population->ind[i]->get_obj2() == 999999 &&
						population->ind[i]->get_obj3() == 0)
		{
			core->evaluates_sol(population->ind[i]->get_gene(), bestOb1, bestOb2, bestOb3);
		
			if(core->check_sol())
			{	
				auxOb1 = core->get_firstObj();
				auxOb2 = core->get_secondObj();
				auxOb3 = core->get_thirdObj();
				
				population->ind[i]->set_obj1(auxOb1);
				population->ind[i]->set_obj2(auxOb2);
				population->ind[i]->set_obj3(auxOb3);
				
				if(core->check_best())
				{
					*bestOb1 = auxOb1;
					*bestOb2 = auxOb2;				
					*bestOb3 = auxOb3;
					core->get_best_genes(bestGenes);
				}
			}else
			{
				population->ind[i]->set_obj1(999999);
				population->ind[i]->set_obj2(999999);
				population->ind[i]->set_obj3(0);
			}
			
			core->model_terminator();
		}	
	}
}

int main(int argc, char **argv)
{
	POPULATIONN *childPop;// new spring
	POPULATIONN *joinPop; // auxiliar population
	POPULATIONN *mainPop; // old spring
	
	POPULATIONNTREATMENT* popT;


	ChronoReal cl;
	double timeReading = 0, execTime = 0;

	cl.start();

	MEMO* core = new MEMO(argv[1]); // To use the model/read the archive



	float crossChance = atof(argv[2]); // Chance of crossover
	float elite = atof(argv[5]);
	
	int doublePop, pop; // size of auxiliar population, size of population
	int geneSize, genQtd; // quantity of cities, quantity of iterations
	
	float bestOb1 = 999999, bestOb2 = 999999, bestOb3 = 0;
	
	

	int* bestGenes;
	int mutationChanceGene, mutationChancePop;
	
	float **dists;
	float df[3], dl[3];
	int stop[3];


	geneSize = core->get_gene_size();
	
	pop = atoi(argv[3]);
	doublePop = 2*pop;



	mutationChanceGene = (int)geneSize*(atof(argv[5]));
	mutationChancePop = (int)pop*(atof(argv[5]));
	
	genQtd = atoi(argv[4]);
	
	bestGenes = (int*)malloc(sizeof(int)*geneSize);
	
	dists = (float**)malloc(sizeof(float*)*3);
	
	for(int i = 0; i < 3; i++){
		dists[i] = (float*)malloc(sizeof(float)*pop);
		for(int j = 0; j < pop; j++)
			dists[i][j] = 0;
	}
		//Create the populations
	mainPop = new POPULATIONN(pop, geneSize, 1);
	childPop = new POPULATIONN(pop, geneSize, 0);
	joinPop = new POPULATIONN(doublePop, geneSize, 0);	
	popT = new POPULATIONNTREATMENT(pop, crossChance);
	
	
	// solve
	solve_pop(core, pop, mainPop, &bestOb1, &bestOb2, &bestOb3, bestGenes);
	

	// classifies
	popT->classifies_population(mainPop, pop);
	
	cl.stop();

	timeReading = cl.getTime();

	cl.reset();

	int step = 0;

	for(int i = 0; i < genQtd; i++)
	{
		cl.start();
	
		//tournment
		popT->selection(mainPop, childPop, pop);
		
		

		popT->pop_mutation(childPop, pop, mutationChancePop, mutationChanceGene);
		
		// solve
		
		solve_pop(core, pop, childPop, &bestOb1, &bestOb2, &bestOb3, bestGenes);
	
		//merge
		
		popT->join_pop(mainPop, childPop, joinPop, pop);
		
		//classifies
		
		popT->classifies_population(joinPop, doublePop);
		
		//get the new  classified pop
		
		popT->update_pop(mainPop, joinPop, pop, doublePop, elite);
		

		//classifies
		
		//popT->classifies_population(mainPop, pop);

		cl.stop();

		execTime += cl.getTime();

		cl.reset();
		

		if( (i == 0) || (i == genQtd/2))
		{
			frontier_registration(mainPop, pop, i, argv[1], dists, step, df, dl, geneSize, stop, core->get_min_prize(), core->get_total_prize());
			step++;
		}
	}

	execTime += timeReading;

	frontier_registration(mainPop, pop, genQtd, argv[1], dists, step, df, dl, geneSize, stop, core->get_min_prize(), core->get_total_prize());
	final_registration(bestOb1, bestOb2, bestOb3, bestGenes, geneSize, argv[1], execTime, dists, pop, df, dl, stop);
	
	delete core;

	delete popT;

	//free(bestGenes);

	//free(mainPop);
	//free(childPop);
	//free(joinPop);

	return 1;
}

