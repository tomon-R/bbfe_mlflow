#include "mlflow_elemmat.h"
#include "BB/calc.h"
#include "BB/std.h"

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
double BBFE_elemmat_mat_levelset(
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

/**********************************************************
 * Surface tension
 **********************************************************/
void BBFE_elemmat_vec_surface_tension(
		const double   N_i,
		const double   grad_N_i[3],
		const double   phi,
		const double   grad_phi[3],
		const double   sigma,
		double*        surf_tension_vec,
		const double   size_interface)
{
	double l_n = BB_calc_vec3d_length(grad_phi);

	double kappa; // curvature
	double delta;
	double alpha = size_interface;

	if( l_n < ZERO_CRITERION){
		kappa = 0;
	}else{
		kappa = BB_calc_vec3d_dot(grad_N_i, grad_phi) / l_n;
	}
	
	//* Conventional LSM
	if(abs(phi) <= alpha){
		delta = (1 + cos(M_PI*phi/alpha))/(2*alpha);
	}else{
		delta = 0;
	}
	for(int d=0; d<3; d++) {
	   	surf_tension_vec[d] = sigma * kappa * grad_phi[d] / l_n * delta;
	}
	//*/

	/* CLSM
	double dist = log(1.0/phi-1) * size_interface;
	if(abs(dist) <= alpha){
		delta = (1 + cos(M_PI*dist/alpha))/(2*alpha);
	}else{
		delta = 0;
	}
	for(int d=0; d<3; d++) {
	   	surf_tension_vec[d] = sigma * kappa * grad_phi[d] / l_n * delta;
	}
	//*/
}

/**********************************************************
 * Vector for Reinitialization for Levelset Method 
 **********************************************************/
double BBFE_elemmat_vec_levelset_reinitialize(
		const double N_i,
		const double phi,
		const double phi_tmp,
		const double grad_phi[3],
		const double dt,
		const double epsilon)
{
	double sign = phi_tmp / sqrt(phi_tmp * phi_tmp + epsilon * epsilon);
	double l_n = BB_calc_vec3d_length(grad_phi);
	double w_vec_grad_phi;
	if( l_n < ZERO_CRITERION){
		w_vec_grad_phi = 0;
	}else{
		w_vec_grad_phi = sign * BB_calc_vec3d_dot(grad_phi, grad_phi) / l_n;
	}

	double val = 0.0;

	val += N_i * phi;
	//val += - N_i * w_vec_grad_phi * dt;
	//val += N_i * sign * dt;
	val += N_i * (sign * (1 - l_n)) * dt;

	return val;
}


/**********************************************************
 * Element Matrix for Reinitialization for Conservative Levelset Method 
 **********************************************************/
double BBFE_elemmat_mat_CLSM_reinitialize(
		const double N_i,
		const double N_j,
		const double grad_N_i[3],
		const double grad_N_j[3],
		const double phi,
		const double grad_phi[3],
		const double dt,
		const double epsilon)
{
	double l_n = BB_calc_vec3d_length(grad_phi);
	double* n_vec;
	n_vec = BB_std_calloc_1d_double(n_vec, 3);
	for(int d=0; d<3; d++){
		n_vec[d] = grad_phi[d]/l_n;
	}

	double val = 0.0;

	val += N_i * N_j;
	val += BB_calc_vec3d_dot(grad_N_i, n_vec) * N_j * (-0.5 + phi) * dt;
	val += epsilon * BB_calc_vec3d_dot(grad_N_j, n_vec) * 0.5 * BB_calc_vec3d_dot(grad_N_i, n_vec) * dt;

	BB_std_free_1d_double(n_vec, 3);
	//if(abs(val)>ZERO_CRITERION)printf("CLSM mat val: %f\n", val);
	return val;
}

/**********************************************************
 * Vector for Reinitialization for Conservative Levelset Method 
 **********************************************************/
double BBFE_elemmat_vec_CLSM_reinitialize(
		const double N_i,
		const double grad_N_i[3],
		const double phi,
		const double normal_vec[3],
		const double grad_phi[3],
		const double dt,
		const double epsilon)
{
	double l_n = BB_calc_vec3d_length(normal_vec);
	double* n_vec;
	n_vec = BB_std_calloc_1d_double(n_vec, 3);
	for(int d=0; d<3; d++){
		n_vec[d] = normal_vec[d]/l_n;
	}

	double val = 0.0;

	val += N_i * phi;
	val += BB_calc_vec3d_dot(grad_N_i, n_vec) * 0.5 * phi * dt;
	val -= epsilon * BB_calc_vec3d_dot(grad_phi, n_vec) * 0.5 * BB_calc_vec3d_dot(grad_N_i, n_vec) * dt;

	BB_std_free_1d_double(n_vec, 3);
	//if(abs(val)>ZERO_CRITERION)printf("CLSM vec val: %f\n", val);
	return val;
}

/**********************************************************
 * Vector for L2 Projection of gradient of levelset
 **********************************************************/
double BBFE_elemmat_vec_grad_phi_L2_projection(
		double vec[3],
		const double N_i,
		const double grad_phi[3])
{
	for(int d=0; d<3; d++) {
		vec[d] = N_i * grad_phi[d];
	}
}