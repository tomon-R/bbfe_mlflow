#pragma once

double elemmat_supg_coef_ml(
		const double v[3],
		const double h_e,
		const double dt);

double elemmat_lsic_coef(
		const double h_e,
		const double v[3]);

double BBFE_elemmat_mat_levelset(
		const double N_i,
		const double N_j,
		const double grad_N_i[3],
		const double v[3],
		const double tau,
		const double v_mesh[3]);

double BBFE_elemmat_vec_levelset(
		double         vec[3],
		const double   N_i,
		const double   grad_N_i[3],
		const double   v[3],
		double**       grad_v,
		const double   phi,
		const double   grad_phi[3],
		const double   density,
		const double   viscosity,
		const double   tau_supg_ml,
		const double   tau_lsic,
		const double   dt,
		const double   v_mesh[3]);

void BBFE_elemmat_vec_surface_tension(
		const double   N_i,
		const double   grad_N_i[3],
		const double   phi,
		const double   grad_phi[3],
		const double   sigma,
		double*        surf_tension_vec,
		const double   size_interface);

double BBFE_elemmat_vec_levelset_reinitialize(
		const double N_i,
		const double phi,
		const double phi_zero,
		const double grad_phi[3],
		const double dt,
		const double epsilon);

double BBFE_elemmat_mat_CLSM_reinitialize(
		const double N_i,
		const double N_j,
		const double grad_N_i[3],
		const double grad_N_j[3],
		const double phi,
		const double grad_phi[3],
		const double dt,
		const double epsilon);

double BBFE_elemmat_vec_CLSM_reinitialize(
		const double N_i,
		const double grad_N_i[3],
		const double phi,
		const double normal_vec[3],
		const double grad_phi[3],
		const double dt,
		const double epsilon);

double BBFE_elemmat_vec_grad_phi_L2_projection(
		double vec[3],
		const double N_i,
		const double grad_phi[3]);