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

double elemmat_supg_coef_ml(
		const double v[3],
		const double h_e,
		const double dt);

double elemmat_lsic_coef(
		const double h_e,
		const double v[3]);

/**********************************************************
 * fractional step method
 **********************************************************/

double elemmat_mat_pred_expl(
		const double N_i,
		const double N_j,
		const double grad_N_i[3],
		const double v[3],
		const double tau);

void elemmat_vec_pred_expl(
		double       vec[3],
		const double N_i,
		const double grad_N_i[3],
		const double v[3],
		double**     grad_v,
		const double density,
		const double viscosity,
		const double tau,
		const double dt,
		const double* gravity);

double elemmat_vec_ppe(
		const double N_i,
		const double div_v,
		const double density,
		const double dt);

void elemmat_vec_corr(
		double       vec[3],
		const double N_i,
		const double grad_p[3],
		const double v[3],
		const double density,
		const double dt);

/**********************************************************
 * levelset
 **********************************************************/
double elemmat_vec_levelset(
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
		const double   dt);

/**********************************************************
 * 2 step fractional step method
 **********************************************************/
void elemmat_mat_pred_expl_1st_step(
	    double vec[3],
		const double N_i,
		const double N_j,
		const double grad_N_i[3],
		const double v[3],
		const double tau,
		double** grad_v,
		const double density,
		const double viscosity,
		const double dt);

void elemmat_vec_pred_expl_1st_step(
		double         vec[3],
		const double   N_i,
		const double   grad_N_i[3],
		const double   v[3],
		double**       grad_v,
		const double   density,
		const double   viscosity,
		const double   tau,
		const double   dt, 
		const double*  gravity,
		const double   p,
		const double   grad_p[3]);

void elemmat_mat_pred_expl_2nd_step(
		double vec[3],
		const double N_i,
		const double N_j);

void elemmat_vec_pred_expl_2nd_step(
		double         vec[3],
		const double   density,
		const double   N_i,
		const double   grad_N_i[3],
		const double   v[3],
		const double   dt,
		const double   p,
		const double   grad_p[3]);