#pragma once

/**********************************************************
 * Stabilized FEM related
 **********************************************************/
double elemmat_supg_coef(
		const double density,
		const double viscosity,
		const double v[3],
		const double h_e,
		const double dt);

/**********************************************************
 * supg + pspg method
 **********************************************************/

double BBFE_elemmat_fluid_sups_coef(
		const double density,
		const double viscosity,
		const double v[3],
		const double h_e,
		const double dt);

void BBFE_elemmat_fluid_sups_mat(
		double         mat[4][4],
		const double   N_i,
		const double   N_j,
		const double   grad_N_i[3],
		const double   grad_N_j[3],
		const double   v[3],
		const double   density,
		const double   viscosity,
		const double   tau,
		const double   dt);

void BBFE_elemmat_fluid_sups_vec(
		double         vec[4],
		const double   N_i,
		const double   grad_N_i[3],
		const double   v[3],
		const double   density,
		const double   tau,
		const double   dt,
		const double*  gravity,
		const double*  surf_tension_vec,
		const double*  accel_inertia);

void BBFE_elemmat_fluid_sups_mat_crank_nicolson(
		double         mat[4][4],
		const double   N_i,
		const double   N_j,
		const double   grad_N_i[3],
		const double   grad_N_j[3],
		const double   v[3],
		const double   density,
		const double   viscosity,
		const double   tau,
		const double   dt);

void BBFE_elemmat_fluid_sups_vec_crank_nicolson(
		double         vec[4],
		const double   N_i,
		const double   grad_N_i[3],
		const double   v[3],
		double**       grad_v,
		const double   density,
		const double   viscosity,
		const double   tau,
		const double   dt,
		const double*  gravity,
		const double*  surf_tension,
		const double*  accel_inertia);