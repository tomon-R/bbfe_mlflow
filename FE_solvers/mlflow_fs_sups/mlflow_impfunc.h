#pragma once

#include "mlflow_elemmat.h"

void BBFE_fluid_renew_vals_by_levelset(
		double* levelset,
		double* val_vec,
		double val_l,
		double val_g,
		const int total_num_nodes);

void BBFE_fluid_convert_levelset2heaviside(
		double* levelset,
		const double mesh_size,
		const int total_num_nodes);

void BBFE_fluid_renew_levelset(
		double*  v,
		double*  ans_vec,
		const int total_num_nodes);
