#pragma once

double elemmat_supg_coef_ml(
		const double v[3],
		const double h_e,
		const double dt);

double elemmat_lsic_coef(
		const double h_e,
		const double v[3]);

double elemmat_mat_levelset(
		const double N_i,
		const double N_j,
		const double grad_N_i[3],
		const double v[3],
		const double tau);

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