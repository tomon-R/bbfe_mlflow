#pragma once

#include "mlflow_elemmat.h"

void BBFE_mlflow_renew_vals_by_levelset(
		double* levelset,
		double* val_vec,
		double val_l,
		double val_g,
		const int total_num_nodes);

void BBFE_mlflow_renew_vals_by_CLSM(
		double* levelset,
		double* val_vec,
		double val_l,
		double val_g,
		const int total_num_nodes);

void BBFE_mlflow_convert_levelset2heaviside(
		double* heaviside,
		double* levelset,
		const double mesh_size,
		const int total_num_nodes);

void BBFE_mlflow_convert_levelset2CLSM(
		double* levelset,
		const double size_interface,
		const int total_num_nodes);

void BBFE_mlflow_renew_levelset(
		double*  v,
		double*  ans_vec,
		const int total_num_nodes);

void BBFE_mlflow_renew_acceleration(
	double* accel, 
	double* accel_amp,
	double* accel_angle_vel,
	double t);