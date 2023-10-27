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
		const double dt);

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

