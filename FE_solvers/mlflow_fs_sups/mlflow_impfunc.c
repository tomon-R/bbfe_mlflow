
#include "mlflow_impfunc.h"

#include <math.h>

/**********************************************************
 * Renew values based on levelset function
 **********************************************************/
void BBFE_fluid_renew_vals_by_levelset(
		double* levelset,
		double* val_vec,
		double val_l,
		double val_g,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++){
		val_vec[i] = 0.5 * (val_l + val_g) + levelset[i] * (val_l - val_g);
	}
}

void BBFE_fluid_convert_levelset2heaviside(
		double* levelset,
		const double mesh_size,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		double h1 = levelset[i]/mesh_size + 1/M_PI*sin(M_PI*levelset[i]/mesh_size);
		if(h1 > 1.0) h1 = 1.0;
		if(h1 < -1.0) h1 = -1.0;
		levelset[i] = 0.5 * h1;
	}
}

void BBFE_fluid_renew_levelset(
		double* v,
		double* ans_vec,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		v[i] = ans_vec[i];
	}
}
