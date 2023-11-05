
#include "pch.h"
#include "utilities.h"
namespace ses {


	class InitialGuessBuilder
	{
	private:
		LocalType calc_mean_on_axis(int numRows, LocalType* pores_in_axis , int* boundaries , int bnd_type) {
			LocalType mean = 0;
			for (int n = 0; n < numRows; n++)
			{
				mean += pores_in_axis[n];
			}
			return mean / numRows;
		}
		LocalType calc_variance_on_axis(int numRows, LocalType* pores_in_axis, int* boundaries , int bnd_type) {
			LocalType var = 0;
			LocalType mean = this->calc_mean_on_axis(numRows, pores_in_axis, boundaries, bnd_type);
			for (int n = 0; n < numRows; n++)
			{
				if(boundaries[n] == bnd_type) // means inlet
				var += (pores_in_axis[n] - mean) * (pores_in_axis[n] - mean);
			}
			return var / numRows;
		}
		void coordinates_fill_vector_linear(int numRows, LocalType* x, LocalType* locActAxis, LocalType from, LocalType to) {
			LocalType quantum = to - from / (LocalType)numRows;
			LocalType accumulation = 0;
			for (int i = 0; i < numRows; i++) {
				x[i] = accumulation;
				accumulation += quantum;
			}
		}
	public:
		//defining options - options are public. feel free to change them
		LocalType COORDINATES_LINEAR_FROM = 0; // pressure starts from 0
		LocalType COORDINATES_LINEAR_TO = 1; // pressure goes to 1
		LocalType COORDINATES_INLET_VARIANCE_TRESHOLD = 0.5; // maximum variance - if inlet nodes variance on any axis break this limit, that axis won't be selected 
		LocalType COORDINATES_OUTLET_VARIANVE_TRESHOLD = 0.5; // maximum variance - if outlet nodes variance on any axis break this limit, that axis won't be selected 
		LocalType COORDINATES_MIN_NODES = 1000; // minimum nodes - if nodes are less than this value, initial guess will be canceled
		// define other initialization algorithms here - give your all declerations a prifix
		void build_with_coordinates(int numRows,int numRowsAct, LocalType* locX, LocalType* locY, LocalType* locZ, LocalType* locActX, LocalType* locActY, LocalType* locActZ, int* bnd, LocalType* x) {
			if (numRowsAct < this->COORDINATES_MIN_NODES) {
				return;
			}
			// calculating variance on x and y and z axis
			LocalType x_inlets_variance = this->calc_variance_on_axis(numRows, locX, bnd, 1); // 1 is for inlets
			LocalType x_outlets_variance = this->calc_variance_on_axis(numRows, locX, bnd, 2); // 2 is for outlets
			if (x_inlets_variance < this->COORDINATES_INLET_VARIANCE_TRESHOLD && x_outlets_variance < this->COORDINATES_OUTLET_VARIANVE_TRESHOLD) {
				//means x is selected as wanted axis
				this->coordinates_fill_vector_linear(numRowsAct, x, locActX, this->COORDINATES_LINEAR_FROM, this->COORDINATES_LINEAR_TO);
				return;

			}
			LocalType y_inlets_variance = this->calc_variance_on_axis(numRows, locY, bnd, 1); // 1 is for inlets
			LocalType y_outlets_variance = this->calc_variance_on_axis(numRows, locY, bnd, 2); // 2 is for outlets
			if (y_inlets_variance < this->COORDINATES_INLET_VARIANCE_TRESHOLD && y_outlets_variance < this->COORDINATES_OUTLET_VARIANVE_TRESHOLD) {
				//means y is selected as wanted axis
				this->coordinates_fill_vector_linear(numRowsAct, x, locActY, this->COORDINATES_LINEAR_FROM, this->COORDINATES_LINEAR_TO);
				return;
			}
			LocalType z_inlets_variance = this->calc_variance_on_axis(numRows, locZ, bnd, 1); // 1 is for inlets
			LocalType z_outlets_variance = this->calc_variance_on_axis(numRows, locZ, bnd, 2); // 2 is for outlets
			if (z_inlets_variance < this->COORDINATES_INLET_VARIANCE_TRESHOLD && z_outlets_variance < this->COORDINATES_OUTLET_VARIANVE_TRESHOLD) {
				//means z is selected as wanted axis
				this->coordinates_fill_vector_linear(numRowsAct, x, locActZ, this->COORDINATES_LINEAR_FROM, this->COORDINATES_LINEAR_TO);
				return;
			}
			
		}

	};
}
