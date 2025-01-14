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

void BBFE_mlflow_renew_mesh_velocity(
		double** v_mesh,
		double* accel_inertia,
		const int total_num_nodes,
		const double dt);

void BBFE_mlflow_renew_mesh_position(
		double** x,
		double** v_mesh,
		const int total_num_nodes,
		const double dt);

void BBFE_mlflow_clear_surface_tension(
		double** surf_tension,
		const int total_num_nodes);