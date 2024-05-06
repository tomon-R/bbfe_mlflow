#include "mlflow_elemmat.h"
#include "BB/calc.h"

#include <math.h>
#include <stdio.h>

const double ZERO_CRITERION  = 1.0e-10;

/**********************************************************
 * Stabilized FEM for levelset convection
 **********************************************************/
double elemmat_supg_coef_ml(
		const double v[3],
		const double h_e,
		const double dt)
{
	double l_v = BB_calc_vec3d_length(v);
	if( l_v < ZERO_CRITERION) { return 0.0; }

	double denom = (2.0/dt)*(2.0/dt) + (2.0*l_v/h_e)*(2.0*l_v/h_e);

	double val = sqrt(1.0/denom);

	return (val);
}

double elemmat_lsic_coef(
		const double h_e,
		const double v[3])
{
	double l_v = BB_calc_vec3d_length(v);
	double val = h_e * l_v / 2;
	return (val);
}

/**********************************************************
 * Element Matrix for levelset convection 
 * This function is same as elemmat_mat_pred_expl(). Only the name is different. 
 **********************************************************/
double elemmat_mat_levelset(
		const double N_i,
		const double N_j,
		const double grad_N_i[3],
		const double v[3],
		const double tau)
{
	double val = 
		N_i * N_j + 
		tau * BB_calc_vec3d_dot(v, grad_N_i) * N_j;

	return val;
}

/**********************************************************
 * Vector for levelset convection
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
		const double   dt)
{
	double val = 0.0;

	val += - N_i *( 
			 v[0] * grad_phi[0] + 
			 v[1] * grad_phi[1] + 
			 v[2] * grad_phi[2] 
			 );

	val += - tau_supg_ml * BB_calc_vec3d_dot(v, grad_N_i) * BB_calc_vec3d_dot(v, grad_phi);

   	val += - tau_lsic * BB_calc_vec3d_dot(grad_N_i, grad_phi);

	val *= dt;

	val += N_i * phi;

	val += tau_supg_ml * BB_calc_vec3d_dot(v, grad_N_i) * phi;

	return val;
}