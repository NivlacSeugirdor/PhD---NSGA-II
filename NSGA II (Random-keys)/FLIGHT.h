#ifndef _FLIGHT_H_
#define _FLIGHT_H_

#include <iostream>
#include <string>
#include "AIRPORT.h"
using namespace std;

// Check
class FLIGHT
{
	public:
		AIRPORT *fromWhere;
		AIRPORT *toWhere;
	 
	 	int departureTime;
	 	int convertDepartureTime;
	 	int arriveTime;
	 	int convertArriveTime;

	 	double price;

	 	int ticket;

	 	AIRPORT *get_from_where(){return fromWhere;}

		void set_from_where(AIRPORT* fromWhere){this->fromWhere = fromWhere;}
		
		AIRPORT *get_to_where(){return toWhere;}

		void set_to_where(AIRPORT* toWhere){this->toWhere = toWhere;}


		int get_converted_departure_time(){return convertDepartureTime;}

		void set_converted_departure_time(int convDepartureTime){this->convertDepartureTime = convDepartureTime;}

		int get_converted_arrive_time(){return convertArriveTime;}

		void set_converted_arrive_time(int convArriveTime){this->convertArriveTime = convArriveTime;}



		int get_departure_time(){return departureTime;}

		void set_departure_time(int departureTime){this->departureTime = departureTime;}

		int get_arrive_time(){return arriveTime;}

		void set_arrive_time(int arriveTime){this->arriveTime = arriveTime;}

		double get_price(){return price;}

		void set_price(int price){this->price = price;}


		int get_ticket(){return ticket;}

		void set_ticket(int ticket){this->price = ticket;}


		int compare_to(FLIGHT another) {
			if(this->price < another.get_price())
				return -1;
			else if(this->price > another.get_price())
				return 1;

			return 0;
		}
};

#endif