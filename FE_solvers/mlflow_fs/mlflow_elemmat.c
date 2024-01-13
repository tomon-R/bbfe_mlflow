
#include "mlflow_elemmat.h"
#include "BB/calc.h"

#include <math.h>
#include <stdio.h>

const double ZERO_CRITERION  = 1.0e-10;

/**********************************************************
 * Stabilized FEM related
 **********************************************************/
double elemmat_supg_coef(
		const double density,
		const double viscosity,
		const double v[3],
		const double h_e,
		const double dt)
{

	double l_v = BB_calc_vec3d_length(v);
	if( l_v < ZERO_CRITERION) { return 0.0; }

	double nu = viscosity/density;
	double denom = (2.0/dt)*(2.0/dt) + (2.0*l_v/h_e)*(2.0*l_v/h_e) 
		+ (4.0*nu/(h_e*h_e))*(4.0*nu/(h_e*h_e));

	double val = sqrt(1.0/denom);

	return (val);
}

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
 * fractional step method
 **********************************************************/

double elemmat_mat_pred_expl(
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


void elemmat_vec_pred_expl(
		double         vec[3],
		const double   N_i,
		const double   grad_N_i[3],
		const double   v[3],
		double**       grad_v,
		const double   density,
		const double   viscosity,
		const double   tau,
		const double   dt, 
		const double*  gravity)
{
	double dyn_vis = viscosity/density;

	for(int d=0; d<3; d++) {
		double val = 0.0;
		
		val += -N_i * ( 
				 v[0]*grad_v[d][0] + 
				 v[1]*grad_v[d][1] + 
				 v[2]*grad_v[d][2] 
				 );
	   /*
	   if(d==0){
			val += -dyn_vis * ( 
				 2*grad_N_i[0]*grad_v[d][0] +
				 grad_N_i[1]*grad_v[d][1] +
				 grad_N_i[2]*grad_v[d][2] +

				 grad_N_i[0]*grad_v[d][1] +
				 grad_N_i[2]*grad_v[d][1] +
				 grad_N_i[0]*grad_v[d][2] +
				 grad_N_i[1]*grad_v[d][2]
				 );
		}else if(d==1){
			val += -dyn_vis * ( 
				 grad_N_i[0]*grad_v[d][0] +
				 2*grad_N_i[1]*grad_v[d][1] +
				 grad_N_i[2]*grad_v[d][2] +

				 grad_N_i[1]*grad_v[d][0] +
				 grad_N_i[2]*grad_v[d][0] +
				 grad_N_i[0]*grad_v[d][2] +
				 grad_N_i[1]*grad_v[d][2]
				 );
		}else if(d==2){
			val += -dyn_vis * ( 
				 grad_N_i[0]*grad_v[d][0] +
				 grad_N_i[1]*grad_v[d][1] +
				 2*grad_N_i[2]*grad_v[d][2] +

				 grad_N_i[1]*grad_v[d][0] +
				 grad_N_i[2]*grad_v[d][0] +
				 grad_N_i[0]*grad_v[d][1] +
				 grad_N_i[2]*grad_v[d][1]
				 );
		}
		//*/
		//*
		val += -dyn_vis * ( 
				 grad_N_i[0]*grad_v[d][0] +
				 grad_N_i[1]*grad_v[d][1] +
				 grad_N_i[2]*grad_v[d][2] 
				 );
		//*/
	   	val += -tau * BB_calc_vec3d_dot(v, grad_N_i) * BB_calc_vec3d_dot(v, grad_v[d]);

	   	val += N_i * gravity[d];

		val *= dt;

		val += N_i * v[d];

		val += tau * BB_calc_vec3d_dot(v, grad_N_i) * v[d];

		vec[d] = val;
	}
}


double elemmat_vec_ppe(
		const double N_i,
		const double div_v,
		const double density,
		const double dt)
{
	double val = density/dt * div_v * N_i;
	
	return val;
}


void elemmat_vec_corr(
		double       vec[3],
		const double N_i,
		const double grad_p[3],
		const double v[3],
		const double density,
		const double dt)
{
	for(int d=0; d<3; d++) {
		vec[d] = -1.0/density * N_i * grad_p[d] * dt + v[d] * N_i;
	}
}

/**********************************************************
 * levelset function convection
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
		const double dt)
{
	double dyn_vis = viscosity/density;

	for(int d=0; d<3; d++){
		double val = 0.0;
		double mass = N_i * N_j;
		double mass_stab = tau * BB_calc_vec3d_dot(v, grad_N_i) * N_j;
		double adv = N_i * ( 
					 v[0]*grad_v[d][0] + 
					 v[1]*grad_v[d][1] + 
					 v[2]*grad_v[d][2] 
					 );
		double adv_stab = tau * BB_calc_vec3d_dot(v, grad_N_i) * BB_calc_vec3d_dot(v, grad_v[d]);
		double diff = dyn_vis * ( 
					 grad_N_i[0]*grad_v[d][0] +
					 grad_N_i[1]*grad_v[d][1] +
					 grad_N_i[2]*grad_v[d][2] 
					 );

		val = mass + mass_stab + 0.5*(adv + adv_stab + diff) * dt;

		vec[d] = val;
	}
}

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
		const double   grad_p[3])
{
	double dyn_vis = viscosity/density;

	for(int d=0; d<3; d++) {
		double val = 0.0;
		
		double adv = -N_i * ( 
				 v[0]*grad_v[d][0] + 
				 v[1]*grad_v[d][1] + 
				 v[2]*grad_v[d][2] 
				 );
	   
		double diff = -dyn_vis * ( 
				 grad_N_i[0]*grad_v[d][0] +
				 grad_N_i[1]*grad_v[d][1] +
				 grad_N_i[2]*grad_v[d][2] 
				 );

	   	double adv_stab = -tau * BB_calc_vec3d_dot(v, grad_N_i) * BB_calc_vec3d_dot(v, grad_v[d]);

	   	double grav = N_i * gravity[d];

		double mass = N_i * v[d];

		double mass_stab = tau * BB_calc_vec3d_dot(v, grad_N_i) * v[d];

		double pres = -1.0/density * N_i * grad_p[d];
		
		double pres_stab = -1.0/density * tau * BB_calc_vec3d_dot(v, grad_N_i) * grad_p[d];

		val = (0.5*(adv + adv_stab + diff) + pres + pres_stab) * dt + (mass + mass_stab) + grav * dt;

		vec[d] = val;
	}
}

void elemmat_mat_pred_expl_2nd_step(
		double vec[3],
		const double N_i,
		const double N_j)
{
	for(int d=0; d<3; d++) {
		double val = 0.0;

		double mass = N_i * N_j;
		val = mass;
		
		vec[d] = val;
	}
}

void elemmat_vec_pred_expl_2nd_step(
		double         vec[3],
		const double   density,
		const double   N_i,
		const double   grad_N_i[3],
		const double   v[3],
		const double   dt,
		const double   p,
		const double   grad_p[3])
{

	for(int d=0; d<3; d++) {
		double val = 0.0;

		double mass = N_i * v[d];
		double pres = 1.0/density * N_i * grad_p[d];

		val = mass + pres* dt;

		vec[d] = val;
	}
}