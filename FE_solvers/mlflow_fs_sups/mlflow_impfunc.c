
#include "mlflow_impfunc.h"

#include <math.h>
#include <stdio.h>

/**********************************************************
 * Renew values based on levelset function
 **********************************************************/
void BBFE_mlflow_renew_vals_by_levelset(
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

void BBFE_mlflow_renew_vals_by_CLSM(
		double* levelset,
		double* val_vec,
		double val_l,
		double val_g,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++){
		val_vec[i] = val_l + (val_g - val_l) * levelset[i];
	}
}

void BBFE_mlflow_convert_levelset2heaviside(
		double* heaviside,
		double* levelset,
		const double size_interface,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		double h1 = levelset[i]/size_interface + 1.0/M_PI*sin(M_PI*levelset[i]/size_interface);
		if(h1 > 1.0) h1 = 1.0;
		if(h1 < -1.0) h1 = -1.0;
		heaviside[i] = 0.5 * h1;
	}
}

void BBFE_mlflow_convert_levelset2CLSM(
		double* levelset,
		const double size_interface,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		levelset[i] = 1.0 / (1.0 + exp(levelset[i]/size_interface));
		//levelset[i] = 0.5 * (tanh(levelset[i]/(2*size_interface)) + 1);
	}
}

void BBFE_mlflow_renew_levelset(
		double* v,
		double* ans_vec,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		v[i] = ans_vec[i];
	}
}