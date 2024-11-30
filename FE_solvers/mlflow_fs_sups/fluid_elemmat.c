
#include "fluid_elemmat.h"
#include "BB/calc.h"
#include "BB/std.h"

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
	double nu = viscosity/density;

	double denom = (2.0/dt)*(2.0/dt) + (2.0*l_v/h_e)*(2.0*l_v/h_e) 
	+ (4.0*nu/(h_e*h_e))*(4.0*nu/(h_e*h_e));

	if(fabs(denom) < ZERO_CRITERION) { return 0.0; }

	double val = sqrt(1.0/denom);
	return (val);
}

/**********************************************************
 * supg + pspg method
 **********************************************************/
//sups = supg + pspg
double BBFE_elemmat_fluid_sups_coef(
		const double density,
		const double viscosity,
		const double v[3],
		const double h_e,
		const double dt)
{
	double l_v = BB_calc_vec3d_length(v);
	double nu = viscosity/density;

	double denom = (2.0/dt)*(2.0/dt) + (2.0*l_v/h_e)*(2.0*l_v/h_e) 
	+ (4.0*nu/(h_e*h_e))*(4.0*nu/(h_e*h_e));

	if(fabs(denom) < ZERO_CRITERION) { return 0.0; }

	double val = sqrt(1.0/denom);
	return (val);
}

/**********************************************************
 * Shock capturing
 **********************************************************/
double BBFE_elemmat_mlflow_shock_capturing_coef(
		const double density,
		const double viscosity,
		const double v[3],
		const double h_e)
{
	double l_v = BB_calc_vec3d_length(v);
	double nu = viscosity / density;
	double re = l_v * h_e / (2 * nu);
	double xi = re/3;
	if(xi>1){
		xi = 1;
	}
	
	double val = (h_e/2) * l_v * xi;
	return (val);
}

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
		const double   tau_c,
		const double   dt,
		const double   v_mesh[3])
{
	double*  v_ale = NULL;
	v_ale = BB_std_calloc_1d_double(v_ale, 3);
	for(int d=0; d<3; d++){
		v_ale[d] = v[d] - v_mesh[d];
	}

	double M = density * N_i * N_j / dt;
	double A = density * N_i * BB_calc_vec3d_dot(v_ale, grad_N_j);

	double G1 = - grad_N_i[0] * N_j;
	double G2 = - grad_N_i[1] * N_j;
	double G3 = - grad_N_i[2] * N_j;

	double D_11 = viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[0]*grad_N_j[0] );
	double D_12 = viscosity * grad_N_i[1] * grad_N_j[0];
	double D_13 = viscosity * grad_N_i[2] * grad_N_j[0];
	double D_21 = viscosity * grad_N_i[0] * grad_N_j[1];
	double D_22 = viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[1]*grad_N_j[1] );
	double D_23 = viscosity * grad_N_i[2] * grad_N_j[1];
	double D_31 = viscosity * grad_N_i[0] * grad_N_j[2];
	double D_32 = viscosity * grad_N_i[1] * grad_N_j[2];
	double D_33 = viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[2]*grad_N_j[2] );

	double C1 = N_i * grad_N_j[0];
	double C2 = N_i * grad_N_j[1];
	double C3 = N_i * grad_N_j[2];

	//SUPG 項
	double M_s = density * tau * BB_calc_vec3d_dot(v_ale, grad_N_i) * N_j / dt;
	double A_s = density * tau * BB_calc_vec3d_dot(v_ale, grad_N_i) * BB_calc_vec3d_dot(v_ale, grad_N_j);

	double G_s1 = tau  * BB_calc_vec3d_dot(v, grad_N_i) * grad_N_j[0];
	double G_s2 = tau  * BB_calc_vec3d_dot(v, grad_N_i) * grad_N_j[1];
	double G_s3 = tau  * BB_calc_vec3d_dot(v, grad_N_i) * grad_N_j[2];

	//PSPG 項
	double M_p1 = tau * grad_N_i[0] * N_j / dt;
	double M_p2 = tau * grad_N_i[1] * N_j / dt;
	double M_p3 = tau * grad_N_i[2] * N_j / dt;

	double A_p1 = tau * grad_N_i[0] * BB_calc_vec3d_dot(v_ale, grad_N_j);
	double A_p2 = tau * grad_N_i[1] * BB_calc_vec3d_dot(v_ale, grad_N_j);
	double A_p3 = tau * grad_N_i[2] * BB_calc_vec3d_dot(v_ale, grad_N_j);

	double G_p = tau * (
			grad_N_i[0] * grad_N_j[0] +
			grad_N_i[1] * grad_N_j[1] +
			grad_N_i[2] * grad_N_j[2]) / density;

	//Shock Capturing項
	double C_s1 = density * tau_c * grad_N_i[0] * grad_N_j[0];
	double C_s2 = density * tau_c * grad_N_i[1] * grad_N_j[1];
	double C_s3 = density * tau_c * grad_N_i[2] * grad_N_j[2];

	mat[0][0] = M + M_s + A + A_s + D_11 + C_s1;
	mat[0][1] = D_12;
	mat[0][2] = D_13;
	mat[0][3] = G1 + G_s1;
	mat[1][0] = D_21;
	mat[1][1] = M + M_s + A + A_s + D_22 + C_s2;
	mat[1][2] = D_23;
	mat[1][3] = G2 + G_s2;
	mat[2][0] = D_31;
	mat[2][1] = D_32;
	mat[2][2] = M + M_s + A + A_s + D_33 + C_s3;
	mat[2][3] = G3 + G_s3;
	mat[3][0] = C1 + M_p1 + A_p1;
	mat[3][1] = C2 + M_p2 + A_p2;
	mat[3][2] = C3 + M_p3 + A_p3;
	mat[3][3] = G_p;

	BB_std_free_1d_double(v_ale, 3);
}


void BBFE_elemmat_fluid_sups_vec(
		double         vec[4],
		const double   N_i,
		const double   grad_N_i[3],
		const double   v[3],
		const double   density,
		const double   tau,
		const double   dt,
		const double*  gravity,
		const double*  surf_tension,
		const double*  accel_inertia,
		const double   v_mesh[3],
		const int      ale_option
		)
{
	double*  v_ale = NULL;
	v_ale = BB_std_calloc_1d_double(v_ale, 3);
	for(int d=0; d<3; d++){
		v_ale[d] = v[d] - v_mesh[d];
	}

	for(int d=0; d<3; d++) {
		double val = 0.0;

		val += density * N_i * v[d];

		val += density * tau * BB_calc_vec3d_dot(v_ale, grad_N_i) * v[d];

		/* external force */
		if(ale_option == 1){
			val += density * N_i * gravity[d] * dt;
			val += (density * N_i * gravity[d] * dt + surf_tension[d] * dt) * tau * BB_calc_vec3d_dot(v_ale, grad_N_i);
		}else{
			val += density * N_i * (gravity[d] + accel_inertia[d]) * dt;
			val += (density * N_i * (gravity[d] + accel_inertia[d]) * dt + surf_tension[d] * dt) * tau * BB_calc_vec3d_dot(v_ale, grad_N_i);
		}
		val += surf_tension[d] * dt;

		vec[d] = val;
	}

	vec[3] = tau * BB_calc_vec3d_dot(grad_N_i, v);

	for(int i=0; i<4; i++){
		vec[i] /= dt;
	}
	BB_std_free_1d_double(v_ale, 3);
}

/**********************************************************
 * A formulation where the Crank-Nicolson method is applied for time discretization 
 *   and the second-order Adams-Bashforth method is applied for linearization with respect to the advection velocity
 **********************************************************/
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
		const double   tau_c,
		const double   dt,
		const double   v_mesh[3])
{
	double*  v_ale = NULL;
	v_ale = BB_std_calloc_1d_double(v_ale, 3);
	for(int d=0; d<3; d++){
		v_ale[d] = v[d] - v_mesh[d];
	}

	double M = density * N_i * N_j;
	double A = dt * density * N_i * BB_calc_vec3d_dot(v_ale, grad_N_j);

	double G1 = - dt * grad_N_i[0] * N_j;
	double G2 = - dt * grad_N_i[1] * N_j;
	double G3 = - dt * grad_N_i[2] * N_j;

	double D_11 = dt * viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[0]*grad_N_j[0] );
	double D_12 = dt * viscosity * grad_N_i[1] * grad_N_j[0];
	double D_13 = dt * viscosity * grad_N_i[2] * grad_N_j[0];
	double D_21 = dt * viscosity * grad_N_i[0] * grad_N_j[1];
	double D_22 = dt * viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[1]*grad_N_j[1] );
	double D_23 = dt * viscosity * grad_N_i[2] * grad_N_j[1];
	double D_31 = dt * viscosity * grad_N_i[0] * grad_N_j[2];
	double D_32 = dt * viscosity * grad_N_i[1] * grad_N_j[2];
	double D_33 = dt * viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_N_j) + grad_N_i[2]*grad_N_j[2] );

	double C1 = dt * N_i * grad_N_j[0];
	double C2 = dt * N_i * grad_N_j[1];
	double C3 = dt * N_i * grad_N_j[2];

	//SUPG 項
	double M_s = density * tau * BB_calc_vec3d_dot(v_ale, grad_N_i) * N_j;
	double A_s = dt * density * tau * BB_calc_vec3d_dot(v_ale, grad_N_i) * BB_calc_vec3d_dot(v_ale, grad_N_j);

	double G_s1 = dt * tau  * BB_calc_vec3d_dot(v, grad_N_i) * grad_N_j[0];
	double G_s2 = dt * tau  * BB_calc_vec3d_dot(v, grad_N_i) * grad_N_j[1];
	double G_s3 = dt * tau  * BB_calc_vec3d_dot(v, grad_N_i) * grad_N_j[2];

	//PSPG 項
	double M_p1 = tau * grad_N_i[0] * N_j;
	double M_p2 = tau * grad_N_i[1] * N_j;
	double M_p3 = tau * grad_N_i[2] * N_j;

	double A_p1 = dt * tau * grad_N_i[0] * BB_calc_vec3d_dot(v_ale, grad_N_j);
	double A_p2 = dt * tau * grad_N_i[1] * BB_calc_vec3d_dot(v_ale, grad_N_j);
	double A_p3 = dt * tau * grad_N_i[2] * BB_calc_vec3d_dot(v_ale, grad_N_j);

	double G_p = dt * tau * (
			grad_N_i[0] * grad_N_j[0] +
			grad_N_i[1] * grad_N_j[1] +
			grad_N_i[2] * grad_N_j[2]) / density;

	//Shock Capturing項
	double C_s1 = dt * density * tau_c * grad_N_i[0] * grad_N_j[0];
	double C_s2 = dt * density * tau_c * grad_N_i[1] * grad_N_j[1];
	double C_s3 = dt * density * tau_c * grad_N_i[2] * grad_N_j[2];

	mat[0][0] = M + M_s + (A + A_s + D_11) * 0.5 + C_s1;
	mat[0][1] = D_12 * 0.5;
	mat[0][2] = D_13 * 0.5;
	mat[0][3] = G1 + G_s1;
	mat[1][0] = D_21 * 0.5;
	mat[1][1] = M + M_s + (A + A_s + D_22) * 0.5 + C_s2;
	mat[1][2] = D_23 * 0.5;
	mat[1][3] = G2 + G_s2;
	mat[2][0] = D_31 * 0.5;
	mat[2][1] = D_32 * 0.5;
	mat[2][2] = M + M_s + (A + A_s + D_33) * 0.5 + C_s3;
	mat[2][3] = G3 + G_s3;
	mat[3][0] = C1 + M_p1 + A_p1 * 0.5;
	mat[3][1] = C2 + M_p2 + A_p2 * 0.5;
	mat[3][2] = C3 + M_p3 + A_p3 * 0.5;
	mat[3][3] = G_p;

	BB_std_free_1d_double(v_ale, 3);
}

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
		const double*  accel_inertia,
		const double   v_mesh[3],
		const int      ale_option)
{
	double*  v_ale = NULL;
	v_ale = BB_std_calloc_1d_double(v_ale, 3);
	for(int d=0; d<3; d++){
		v_ale[d] = v[d] - v_mesh[d];
	}

	for(int d=0; d<3; d++) {
		double val = 0.0;

		val += density * N_i * v[d]; // M

		val += density * tau * BB_calc_vec3d_dot(v_ale, grad_N_i) * v[d]; // M_s

		val += density * dt * 0.5 * (-N_i) * BB_calc_vec3d_dot(v_ale, grad_v[d]); // A/2

		val += density * dt * 0.5 * (-tau) * BB_calc_vec3d_dot(v_ale, grad_N_i) * BB_calc_vec3d_dot(v_ale, grad_v[d]); // A_s/2

		// D
		double D_11 = viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_v[d]) + grad_N_i[0]*grad_v[d][0] );
		double D_12 = viscosity * grad_N_i[1] * grad_v[d][0];
		double D_13 = viscosity * grad_N_i[2] * grad_v[d][0];
		double D_21 = viscosity * grad_N_i[0] * grad_v[d][1];
		double D_22 = viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_v[d]) + grad_N_i[1]*grad_v[d][1] );
		double D_23 = viscosity * grad_N_i[2] * grad_v[d][1];
		double D_31 = viscosity * grad_N_i[0] * grad_v[d][2];
		double D_32 = viscosity * grad_N_i[1] * grad_v[d][2];
		double D_33 = viscosity * ( BB_calc_vec3d_dot(grad_N_i, grad_v[d]) + grad_N_i[2]*grad_v[d][2] );

		if(d==0){
			val += - density * dt * 0.5 * (D_11 + D_12 + D_13);
		}else if(d==1){
			val += - density * dt * 0.5 * (D_21 + D_22 + D_23);
		}else if(d==2){
			val += - density * dt * 0.5 * (D_31 + D_32 + D_33);
		}
		/* external force */
		if(ale_option == 1){
			val += density * N_i * gravity[d] * dt;
			val += (density * N_i * gravity[d] * dt + surf_tension[d] * dt) * tau * BB_calc_vec3d_dot(v_ale, grad_N_i);
		}else{
			val += density * N_i * (gravity[d] + accel_inertia[d]) * dt;
			val += (density * N_i * (gravity[d] + accel_inertia[d]) * dt + surf_tension[d] * dt) * tau * BB_calc_vec3d_dot(v_ale, grad_N_i);
		}
		val += surf_tension[d] * dt;

		vec[d] = val;
	}
	
	vec[3] = 0;

	vec[3] += tau * BB_calc_vec3d_dot(grad_N_i, v); // Mp1

	double A_p1 = tau * grad_N_i[0] * BB_calc_vec3d_dot(v_ale, grad_v[0]);
	double A_p2 = tau * grad_N_i[1] * BB_calc_vec3d_dot(v_ale, grad_v[1]);
	double A_p3 = tau * grad_N_i[2] * BB_calc_vec3d_dot(v_ale, grad_v[2]);

	vec[3] += - dt * 0.5 * (A_p1 + A_p2 + A_p3) ; //Ap

	BB_std_free_1d_double(v_ale, 3);
}