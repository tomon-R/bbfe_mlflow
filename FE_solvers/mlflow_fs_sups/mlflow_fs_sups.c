
#include "mlflow_impfunc.h"
#include "fluid_core.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 2;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-8;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 10000;
const char*              ID_DT  = "#time_spacing";
const double           DVAL_DT  = 0.01;
const char*     ID_FINISH_TIME  = "#finish_time";
const double  DVAL_FINISH_TIME  = 1.0;
const char* ID_OUTPUT_INTERVAL  = "#output_interval";
const int DVAL_OUTPUT_INTERVAL  = 1;
const char*       ID_DENSITY_L  = "#density_l";
const double    DVAL_DENSITY_L  = 1000.0;
const char*       ID_DENSITY_G  = "#density_g";
const double    DVAL_DENSITY_G  = 1000.0;
const char*     ID_VISCOSITY_L  = "#viscosity_l";
const double  DVAL_VISCOSITY_L  = 1.0;
const char*     ID_VISCOSITY_G  = "#viscosity_g";
const double  DVAL_VISCOSITY_G  = 1.0;
const char*         ID_GRAVITY  = "#gravity";
const double      DVAL_GRAVITY  = 0.0;
const char*    ID_SURF_TENSION_COEF = "#surf_tension_coef";
const double DVAL_SURF_TENSION_COEF = 0.0;
const char*    ID_SIZE_INTERFACE = "#size_interface";
const double DVAL_SIZE_INTERFACE = 0.1;
const char*         ID_DT_REINIT = "#dt_reinit";
const double      DVAL_DT_REINIT = 0.2;
const char*    ID_EPSILON_REINIT = "#epsilon_reinit";
const double DVAL_EPSILON_REINIT = 0.5;
const char*      ID_DELTA_REINIT = "#delta_reinit";
const double   DVAL_DELTA_REINIT = 1e-1;

const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc_p.dat";
static const char* INPUT_FILENAME_LEVELSET= "levelset.dat";

static const char* OUTPUT_FILENAME_VTK    = "result_%d_%06d.vtk";
static const char* OUTPUT_FILENAME_DAMBREAK = "result_dambreak.csv";
static const char* OUTPUT_FILENAME_BUBBLE = "result_bubble.csv";


typedef struct
{
	int    num_ip_each_axis;
	double mat_epsilon;
	int    mat_max_iter;

	double dt;
	double finish_time;
	int    output_interval;

	double density_l, density_g; 
	double viscosity_l, viscosity_g;
	double* gravity;

	double** v;
	double*  p;
	double*  density;
	double*  viscosity;

	double*  levelset;
	double*  levelset_tmp;
	double*  heaviside;

	double surf_tension_coef;
	double** surf_tension;
	double** grad_phi;
	double size_interface;

	double dt_reinit;
	double epsilon_reinit;
	double delta_reinit;

} VALUES;


typedef struct
{
	const char* directory;

} CONDITIONS;


typedef struct
{
	BBFE_BASIS   basis;
	BBFE_DATA    fe;

	CONDITIONS   cond;
	VALUES       vals;

	BBFE_BC      bc;

	MONOLIS      monolis;
	MONOLIS      mono_levelset;
	MONOLIS      mono_reinit;
	MONOLIS      mono_L2;

	MONOLIS_COM  mono_com;

} FE_SYSTEM;


void memory_allocation_nodal_values(
		VALUES*         vals,
		const int       total_num_nodes)
{
	vals->v = BB_std_calloc_2d_double(vals->v, total_num_nodes, 3);
	vals->p = BB_std_calloc_1d_double(vals->p, total_num_nodes);
	vals->density   = BB_std_calloc_1d_double(vals->density, total_num_nodes);
	vals->viscosity = BB_std_calloc_1d_double(vals->viscosity, total_num_nodes);
	vals->levelset  = BB_std_calloc_1d_double(vals->levelset, total_num_nodes);
	vals->levelset_tmp  = BB_std_calloc_1d_double(vals->levelset_tmp, total_num_nodes);
	vals->heaviside  = BB_std_calloc_1d_double(vals->heaviside, total_num_nodes);
	vals->surf_tension  = BB_std_calloc_2d_double(vals->surf_tension, total_num_nodes, 3);
	vals->grad_phi  = BB_std_calloc_2d_double(vals->grad_phi, total_num_nodes, 3);
}


void assign_default_values(
		VALUES*     vals)
{
	vals->num_ip_each_axis = DVAL_NUM_IP_EACH_AXIS;
	vals->mat_epsilon      = DVAL_MAT_EPSILON;
	vals->mat_max_iter     = DVAL_MAT_MAX_ITER;

	vals->dt               = DVAL_DT;
	vals->finish_time      = DVAL_FINISH_TIME;
	vals->output_interval  = DVAL_OUTPUT_INTERVAL;

	vals->density_l        = DVAL_DENSITY_L;
	vals->density_g        = DVAL_DENSITY_G;
	vals->viscosity_l      = DVAL_VISCOSITY_L;
	vals->viscosity_g      = DVAL_VISCOSITY_G;
	
	vals->gravity = BB_std_calloc_1d_double(vals->gravity, 3);
	vals->gravity[0] = DVAL_GRAVITY;
	vals->gravity[1] = DVAL_GRAVITY;
	vals->gravity[2] = DVAL_GRAVITY;

	vals->surf_tension_coef = DVAL_SURF_TENSION_COEF;

	vals->size_interface = DVAL_SIZE_INTERFACE;

	vals->dt_reinit = DVAL_DT_REINIT;
	vals->epsilon_reinit = DVAL_EPSILON_REINIT;
	vals->delta_reinit = DVAL_DELTA_REINIT;
}


void print_all_values(
		VALUES*  vals)
{
	printf("\n%s ---------- Calculation condition ----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, ID_NUM_IP_EACH_AXIS, vals->num_ip_each_axis);
	printf("%s %s: %e\n", CODENAME, ID_MAT_EPSILON,      vals->mat_epsilon);
	printf("%s %s: %d\n", CODENAME, ID_MAT_MAX_ITER,     vals->mat_max_iter);

	printf("%s %s: %e\n", CODENAME, ID_DT,               vals->dt);
	printf("%s %s: %e\n", CODENAME, ID_FINISH_TIME,      vals->finish_time);
	printf("%s %s: %d\n", CODENAME, ID_OUTPUT_INTERVAL,  vals->output_interval);

	printf("%s %s: %e\n", CODENAME, ID_DENSITY_L,          vals->density_l);
	printf("%s %s: %e\n", CODENAME, ID_DENSITY_G,          vals->density_g);
	printf("%s %s: %e\n", CODENAME, ID_VISCOSITY_L,        vals->viscosity_l);
	printf("%s %s: %e\n", CODENAME, ID_VISCOSITY_G,        vals->viscosity_g);
	printf("%s %s: %e %e %e\n", CODENAME, ID_GRAVITY, vals->gravity[0], vals->gravity[1], vals->gravity[2]);
	printf("%s %s: %e\n", CODENAME, ID_SURF_TENSION_COEF, vals->surf_tension_coef);
	printf("%s %s: %e\n", CODENAME, ID_SIZE_INTERFACE, vals->size_interface);

	printf("%s %s: %e\n", CODENAME, ID_DT_REINIT, vals->dt_reinit);
	printf("%s %s: %e\n", CODENAME, ID_EPSILON_REINIT, vals->epsilon_reinit);
	printf("%s %s: %e\n", CODENAME, ID_DELTA_REINIT, vals->delta_reinit);

	printf("%s -------------------------------------------\n\n", CODENAME);
}


void read_calc_conditions(
		VALUES*     vals,
		const char* directory)
{
	printf("\n");

	assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading conditon file \"%s\".\n", CODENAME, filename);
		int num;
		num = BB_std_read_file_get_val_int_p(
				&(vals->num_ip_each_axis), filename, ID_NUM_IP_EACH_AXIS, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->mat_epsilon), filename, ID_MAT_EPSILON, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_int_p(
				&(vals->mat_max_iter), filename, ID_MAT_MAX_ITER, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->dt), filename, ID_DT, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->finish_time), filename, ID_FINISH_TIME, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_int_p(
				&(vals->output_interval), filename, ID_OUTPUT_INTERVAL, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->density_l), filename, ID_DENSITY_L, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->density_g), filename, ID_DENSITY_G, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->viscosity_l), filename, ID_VISCOSITY_L, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->viscosity_g), filename, ID_VISCOSITY_G, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				vals->gravity, filename, ID_GRAVITY, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->surf_tension_coef), filename, ID_SURF_TENSION_COEF, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->size_interface), filename, ID_SIZE_INTERFACE, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->dt_reinit), filename, ID_DT_REINIT, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->epsilon_reinit), filename, ID_EPSILON_REINIT, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->delta_reinit), filename, ID_DELTA_REINIT, BUFFER_SIZE, CODENAME);
		fclose(fp);
	}

	print_all_values(vals);


	printf("\n");
}

void read_levelset_file(
		BBFE_DATA*   fe,
		VALUES*      vals,
		const char*  filename,
		const char*  directory)
{
	FILE* fp;
	fp = BBFE_sys_read_fopen(fp, filename, directory);
	char* label[256];
	int num;

	// read label
	BB_std_scan_line(
			&fp, BUFFER_SIZE, "%s", &(label));
	printf("%s Label: %s\n", CODENAME, label);

	// read the number of nodes
	BB_std_scan_line(
			&fp, BUFFER_SIZE, "%d %d", &(fe->total_num_nodes), &(num));

	printf("%s Num. nodes: %d\n", CODENAME, fe->total_num_nodes);
	printf("%s Num. values per nodes: %d\n", CODENAME, num);

	// read levelset value of nodes
	for(int i=0; i<(fe->total_num_nodes); i++) {
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%lf", &(vals->levelset[i]));
	}

	fclose(fp);
}


void output_result_file_vtk(
		BBFE_DATA*       fe,
		VALUES*        vals,
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);
	printf("total_num_nodes: %d\n", fe->total_num_nodes);

	switch( fe->local_num_nodes ) {
		case 4:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
			break;

		case 8:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
			break;
	}

	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	BB_vtk_write_point_vals_vector(fp, vals->v, fe->total_num_nodes, "Velocity");
	BB_vtk_write_point_vals_scalar(fp, vals->p, fe->total_num_nodes, "Pressure");
	BB_vtk_write_point_vals_scalar(fp, vals->levelset, fe->total_num_nodes, "Levelset");
	BB_vtk_write_point_vals_scalar(fp, vals->heaviside, fe->total_num_nodes, "Heaviside");
	BB_vtk_write_point_vals_scalar(fp, vals->density, fe->total_num_nodes, "Density");
	BB_vtk_write_point_vals_scalar(fp, vals->viscosity, fe->total_num_nodes, "Viscosity");
	BB_vtk_write_point_vals_vector(fp, vals->surf_tension, fe->total_num_nodes, "SurfaceTension");
	BB_vtk_write_point_vals_vector(fp, vals->grad_phi, fe->total_num_nodes, "GradientLevelset");
	fclose(fp);

}


void output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	int myrank = monolis_mpi_get_global_my_rank();
	char fname_vtk[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, myrank, file_num);

	output_result_file_vtk(
			&(sys->fe), &(sys->vals), fname_vtk, sys->cond.directory, t);
}

void BBFE_fluid_sups_read_Dirichlet_bc(
		BBFE_BC*     bc,
		const char*  filename,
		const char*  directory,
		const int    total_num_nodes,
		const int    block_size)
{
	bc->total_num_nodes = total_num_nodes;
	bc->block_size      = block_size;

	BBFE_sys_memory_allocation_Dirichlet_bc(bc, total_num_nodes, bc->block_size);
	int n = total_num_nodes * bc->block_size;

	for(int i=0; i<n; i++) {
		bc->D_bc_exists[i]   = false;
		bc->imposed_D_val[i] = 0.0;
	}

	FILE* fp;
	fp = BBFE_sys_read_fopen_without_error(fp, filename, directory);
	if( fp == NULL ) {
		printf("%s WARNING: Dirichlet B.C. file, \"%s\", is not found.\n",
				CODENAME, filename);
		return;
	}

	int tmp;
	BB_std_scan_line(&fp, BUFFER_SIZE,
			"%d %d", &(bc->num_D_bcs), &(tmp));
	printf("%s Num. Dirichlet B.C.: %d, Num. block size: %d\n", CODENAME, bc->num_D_bcs, tmp);

	for(int i=0; i<(bc->num_D_bcs); i++) {
		int node_id;  int block_id;  double val;
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%d %d %lf", &node_id, &block_id, &val);

		int index = (bc->block_size)*node_id + block_id;
		bc->D_bc_exists[ index ]   = true;
		bc->imposed_D_val[ index ] = val;
	}

	fclose(fp);
}


void set_element_mat(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double*** val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_3d_double(val_ip     , 4 , 4, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double* local_viscosity;
	local_viscosity = BB_std_calloc_1d_double(local_viscosity, nl);
	double* local_density;
	local_density = BB_std_calloc_1d_double(local_density, nl);
	double* viscosity_ip;
	viscosity_ip = BB_std_calloc_1d_double(viscosity_ip, np);
	double* density_ip;
	density_ip = BB_std_calloc_1d_double(density_ip, np);

	double A[4][4];

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		BBFE_elemmat_set_local_array_scalar(local_viscosity, fe, vals->viscosity, e);
		BBFE_elemmat_set_local_array_scalar(local_density, fe, vals->density, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
			viscosity_ip[p] = BBFE_std_mapping_scalar(nl, local_viscosity, basis->N[p]);
			density_ip[p] = BBFE_std_mapping_scalar(nl, local_density, basis->N[p]);		
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {

					double tau = BBFE_elemmat_fluid_sups_coef(
							density_ip[p], viscosity_ip[p], v_ip[p], h_e, vals->dt);

					BBFE_elemmat_fluid_sups_mat(
							A, basis->N[p][i], basis->N[p][j], 
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], 
							v_ip[p], density_ip[p], viscosity_ip[p], tau, vals->dt);

					for(int a=0; a<4; a++){
						for(int b=0; b<4; b++) {
							val_ip[a][b][p] = A[a][b];
							A[a][b] = 0.0;
						}
					}
				}

				for(int a=0; a<4; a++){
					for(int b=0; b<4; b++) {
						double integ_val = BBFE_std_integ_calc(
								np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

						monolis_add_scalar_to_sparse_matrix_R(
								monolis, fe->conn[e][i], fe->conn[e][j], a, b, integ_val);
					}
				}
			}
		}
	}

	BB_std_free_3d_double(val_ip     , 4 , 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);

	// for multilayer flow
	BB_std_free_1d_double(local_density, nl);
	BB_std_free_1d_double(local_viscosity, nl);
	BB_std_free_1d_double(density_ip, np);
	BB_std_free_1d_double(viscosity_ip, np);
}


void set_element_vec(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double** val_ip;
	double*  Jacobian_ip;
	val_ip      = BB_std_calloc_2d_double(val_ip, 4, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double* local_viscosity;
	local_viscosity = BB_std_calloc_1d_double(local_viscosity, nl);
	double* local_density;
	local_density = BB_std_calloc_1d_double(local_density, nl);
	double* local_levelset;
	local_levelset = BB_std_calloc_1d_double(local_levelset, nl);
	double** local_grad_phi;
	local_grad_phi = BB_std_calloc_2d_double(local_grad_phi, nl, 3);

	double* viscosity_ip;
	viscosity_ip = BB_std_calloc_1d_double(viscosity_ip, np);
	double* density_ip;
	density_ip = BB_std_calloc_1d_double(density_ip, np);
	double* levelset_ip;
	levelset_ip = BB_std_calloc_1d_double(levelset_ip, np);
	double** grad_phi_ip;
	grad_phi_ip = BB_std_calloc_2d_double(grad_phi_ip, np, 3);

	double** surf_tension_ip; double** surf_tension_ip2;
	surf_tension_ip = BB_std_calloc_2d_double(val_ip, np, 3);
	surf_tension_ip2 = BB_std_calloc_2d_double(val_ip, 3, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(
				np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		BBFE_elemmat_set_local_array_scalar(local_viscosity, fe, vals->viscosity, e);
		BBFE_elemmat_set_local_array_scalar(local_density, fe, vals->density, e);
		BBFE_elemmat_set_local_array_scalar(local_levelset, fe, vals->levelset, e);
		BBFE_elemmat_set_local_array_vector(local_grad_phi, fe, vals->grad_phi, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
			BBFE_std_mapping_scalar_grad(grad_phi_ip[p], nl, local_levelset, fe->geo[e][p].grad_N);

			viscosity_ip[p] = BBFE_std_mapping_scalar(nl, local_viscosity, basis->N[p]);
			density_ip[p] = BBFE_std_mapping_scalar(nl, local_density, basis->N[p]);
			levelset_ip[p]  = BBFE_std_mapping_scalar(nl, local_levelset, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			double integ_val[4];

			for(int p=0; p<np; p++) {
				double tau = BBFE_elemmat_fluid_sups_coef(
						density_ip[p], viscosity_ip[p], v_ip[p], h_e, vals->dt);

				BBFE_elemmat_vec_surface_tension(
						basis->N[p][i], fe->geo[e][p].grad_N[i], levelset_ip[p], grad_phi_ip[p], 
						vals->surf_tension_coef, surf_tension_ip[p], vals->epsilon_reinit);

				double vec[4];
				BBFE_elemmat_fluid_sups_vec(
						vec, basis->N[p][i], fe->geo[e][p].grad_N[i],
						v_ip[p], density_ip[p], tau, vals->dt, vals->gravity, surf_tension_ip[p]);

				for(int d=0; d<4; d++) {
					val_ip[d][p] = vec[d];
				}
				for(int d=0; d<3; d++){
					surf_tension_ip2[d][p] = surf_tension_ip[p][d];
				}
			}

			for(int d=0; d<4; d++) {
				integ_val[d] = BBFE_std_integ_calc(
						np, val_ip[d], basis->integ_weight, Jacobian_ip);

				monolis->mat.R.B[ 4*fe->conn[e][i] + d ] += integ_val[d];
			}

			for(int d=0; d<3; d++) {
				vals->surf_tension[fe->conn[e][i]][d] += BBFE_std_integ_calc(np, surf_tension_ip2[d], basis->integ_weight, Jacobian_ip);
			}
		}
	}

	BB_std_free_2d_double(val_ip, 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);

	// for multilayer flow
	BB_std_free_1d_double(local_density, nl);
	BB_std_free_1d_double(local_viscosity, nl);
	BB_std_free_1d_double(density_ip, np);
	BB_std_free_1d_double(viscosity_ip, np);

	BB_std_free_1d_double(local_levelset, nl);
	BB_std_free_1d_double(levelset_ip, np);
	BB_std_free_2d_double(local_grad_phi, nl, 3);
	BB_std_free_2d_double(grad_phi_ip, np, 3);

	BB_std_free_2d_double(surf_tension_ip, np, 3);
	BB_std_free_2d_double(surf_tension_ip2, 3, np);
}

/**********************************************************
 * Matrix and Vector for Levelset 
 **********************************************************/
void set_element_mat_levelset(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double* local_viscosity;
	local_viscosity = BB_std_calloc_1d_double(local_viscosity, nl);
	double* local_density;
	local_density = BB_std_calloc_1d_double(local_density, nl);
	double* viscosity_ip;
	viscosity_ip = BB_std_calloc_1d_double(viscosity_ip, np);
	double* density_ip;
	density_ip = BB_std_calloc_1d_double(density_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		BBFE_elemmat_set_local_array_scalar(local_viscosity, fe, vals->viscosity, e);
		BBFE_elemmat_set_local_array_scalar(local_density, fe, vals->density, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
			viscosity_ip[p] = BBFE_std_mapping_scalar(nl, local_viscosity, basis->N[p]);
			density_ip[p] = BBFE_std_mapping_scalar(nl, local_density, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					double tau = elemmat_supg_coef_ml(v_ip[p], h_e, vals->dt);

					val_ip[p] = BBFE_elemmat_mat_levelset(
							basis->N[p][i], basis->N[p][j], fe->geo[e][p].grad_N[i], v_ip[p], tau);
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix_R(
						monolis, fe->conn[e][i], fe->conn[e][j], 0, 0, integ_val);
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);

	// for multilayer flow
	BB_std_free_1d_double(local_density, nl);
	BB_std_free_1d_double(local_viscosity, nl);
	BB_std_free_1d_double(density_ip, np);
	BB_std_free_1d_double(viscosity_ip, np);
}

void set_element_vec_levelset(
		MONOLIS*	monolis,
		BBFE_DATA*	fe,
		BBFE_BASIS* basis,
		VALUES*		vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double*  val_ip;
	double*  Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	// for multilayer flow
	double* local_levelset;
	local_levelset = BB_std_calloc_1d_double(local_levelset, nl);
	double* local_viscosity;
	local_viscosity = BB_std_calloc_1d_double(local_viscosity, nl);
	double* local_density;
	local_density = BB_std_calloc_1d_double(local_density, nl);
	double* levelset_ip;
	levelset_ip = BB_std_calloc_1d_double(levelset_ip, np);
	double* viscosity_ip;
	viscosity_ip = BB_std_calloc_1d_double(viscosity_ip, np);
	double* density_ip;
	density_ip = BB_std_calloc_1d_double(density_ip, np);
	double** grad_phi_ip;
	grad_phi_ip = BB_std_calloc_2d_double(grad_phi_ip, np, 3);

	double*** grad_v_ip;
	grad_v_ip = BB_std_calloc_3d_double(grad_v_ip, np, 3, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(
				np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
		// for multilayer flow
		BBFE_elemmat_set_local_array_scalar(local_levelset, fe, vals->levelset, e);
		BBFE_elemmat_set_local_array_scalar(local_viscosity, fe, vals->viscosity, e);
		BBFE_elemmat_set_local_array_scalar(local_density, fe, vals->density, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
			BBFE_std_mapping_vector3d_grad(grad_v_ip[p], nl, local_v, fe->geo[e][p].grad_N);
			BBFE_std_mapping_scalar_grad(grad_phi_ip[p], nl, local_levelset, fe->geo[e][p].grad_N);
			// for multilayer flow
			levelset_ip[p]  = BBFE_std_mapping_scalar(nl, local_levelset, basis->N[p]);
			viscosity_ip[p] = BBFE_std_mapping_scalar(nl, local_viscosity, basis->N[p]);
			density_ip[p]   = BBFE_std_mapping_scalar(nl, local_density, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				double tau_supg_ml = elemmat_supg_coef_ml(v_ip[p], h_e, vals->dt);
				double tau_lsic    = elemmat_lsic_coef(h_e, v_ip[p]);

				double vec[3];
				val_ip[p] = BBFE_elemmat_vec_levelset(
					vec, 
					basis->N[p][i], fe->geo[e][p].grad_N[i], 
					v_ip[p], grad_v_ip[p],
					levelset_ip[p], grad_phi_ip[p],
					density_ip[p], viscosity_ip[p], tau_supg_ml, tau_lsic, vals->dt);
			}
			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.R.B[ fe->conn[e][i] ] += integ_val;
		}
	}
	
	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_3d_double(grad_v_ip, np, 3, 3);
	// for multilayer flow
	BB_std_free_1d_double(local_density, nl);
	BB_std_free_1d_double(local_viscosity, nl);
	BB_std_free_1d_double(local_levelset, nl);
	BB_std_free_1d_double(density_ip, np);
	BB_std_free_1d_double(viscosity_ip, np);
	BB_std_free_1d_double(levelset_ip, np);
	BB_std_free_2d_double(grad_phi_ip, np, 3);
}

/**********************************************************
 * L2 projection of gradient of levelset
 **********************************************************/
void set_element_vec_L2_projection(
		MONOLIS*    monolis,
		BBFE_DATA*	fe,
		BBFE_BASIS* basis,
		VALUES*		vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double** val_ip;
	val_ip      = BB_std_calloc_2d_double(val_ip, 3, np);
	double*  Jacobian_ip;
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double* local_levelset;
	local_levelset = BB_std_calloc_1d_double(local_levelset, nl);
	double* levelset_ip;
	levelset_ip = BB_std_calloc_1d_double(levelset_ip, np);
	double** grad_phi_ip;
	grad_phi_ip = BB_std_calloc_2d_double(grad_phi_ip, np, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(
				np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_scalar(local_levelset, fe, vals->levelset, e);

		for(int p=0; p<np; p++) {
			levelset_ip[p]  = BBFE_std_mapping_scalar(nl, local_levelset, basis->N[p]);
			BBFE_std_mapping_scalar_grad(grad_phi_ip[p], nl, local_levelset, fe->geo[e][p].grad_N);

		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				double vec[3];
				
				BBFE_elemmat_vec_grad_phi_L2_projection(
						vec, basis->N[p][i], grad_phi_ip[p]);

				for(int d=0; d<3; d++) {
					val_ip[d][p] = vec[d];
				}
			}
			double integ_val[3];
			for(int d=0; d<3; d++) {
				integ_val[d] = BBFE_std_integ_calc(
						np, val_ip[d], basis->integ_weight, Jacobian_ip);

				monolis->mat.R.B[ 3*(fe->conn[e][i]) + d ] += integ_val[d];
			}
		}
	}
	
	BB_std_free_2d_double(val_ip, 3, np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_1d_double(local_levelset, nl);
	BB_std_free_1d_double(levelset_ip, np);
	BB_std_free_2d_double(grad_phi_ip, np, 3);
}

/**********************************************************
 * Reinitialization of Levelset
 **********************************************************/
void set_element_vec_levelset_reinitialize(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double* local_levelset;
	local_levelset = BB_std_calloc_1d_double(local_levelset, nl);
	double* levelset_ip;
	levelset_ip = BB_std_calloc_1d_double(levelset_ip, np);

	double* local_levelset_tmp;
	local_levelset_tmp = BB_std_calloc_1d_double(local_levelset_tmp, nl);
	double* levelset_ip_tmp;
	levelset_ip_tmp = BB_std_calloc_1d_double(levelset_ip_tmp, np);

	double** local_grad_phi;
	local_grad_phi = BB_std_calloc_2d_double(local_grad_phi, nl, 3);
	double** grad_phi_ip;
	grad_phi_ip = BB_std_calloc_2d_double(grad_phi_ip, np, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		BBFE_elemmat_set_local_array_scalar(local_levelset, fe, vals->levelset, e);
		BBFE_elemmat_set_local_array_scalar(local_levelset_tmp, fe, vals->levelset_tmp, e);
		BBFE_elemmat_set_local_array_vector(local_grad_phi, fe, vals->grad_phi, e, 3);

		for(int p=0; p<np; p++) {
			levelset_ip[p]  = BBFE_std_mapping_scalar(nl, local_levelset, basis->N[p]);
			levelset_ip_tmp[p]  = BBFE_std_mapping_scalar(nl, local_levelset_tmp, basis->N[p]);
			//BBFE_std_mapping_vector3d(grad_phi_ip[p], nl, local_grad_phi, basis->N[p]);
			BBFE_std_mapping_scalar_grad(grad_phi_ip[p], np, local_levelset, fe->geo[e][p].grad_N);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				val_ip[p] = BBFE_elemmat_vec_levelset_reinitialize(
					basis->N[p][i], levelset_ip[p], levelset_ip_tmp[p], grad_phi_ip[p], vals->dt_reinit, vals->epsilon_reinit);
			}
			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.R.B[ fe->conn[e][i] ] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_1d_double(local_levelset, nl);
	BB_std_free_1d_double(local_levelset_tmp, nl);
	BB_std_free_1d_double(levelset_ip, np);
	BB_std_free_1d_double(levelset_ip_tmp, np);

	BB_std_free_2d_double(local_grad_phi, nl, 3);
	BB_std_free_2d_double(grad_phi_ip, np, 3);
}

void reinit_levelset(
		FE_SYSTEM*  sys)
{
	double t = 0.0;
	int step = 0;

	for(int i=0; i<sys->fe.total_num_nodes; i++){
		sys->vals.levelset_tmp[i] = sys->vals.levelset[i];
	}

	while (step<5) {
		t += sys->vals.dt_reinit;
		step += 1;
		monolis_clear_mat_value_R(&(sys->mono_reinit));
		monolis_clear_mat_value_rhs_R(&(sys->mono_reinit));

		BBFE_elemmat_set_global_mat_cmass_const(
					&(sys->mono_reinit),
					&(sys->fe),
					&(sys->basis), 
					1.0, 
					1);
		set_element_vec_levelset_reinitialize(
					&(sys->mono_reinit),
					&(sys->fe),
					&(sys->basis),
					&(sys->vals));
		BBFE_sys_monowrap_solve(
					&(sys->mono_reinit),
					&(sys->mono_com),
					sys->mono_reinit.mat.R.X,
					MONOLIS_ITER_BICGSTAB,
					MONOLIS_PREC_DIAG,
					sys->vals.mat_max_iter,
					sys->vals.mat_epsilon);

		double error = 0.0;
		for(int i=0; i<(sys->fe.total_num_nodes); i++) {
			// Residual
			error += (sys->mono_reinit.mat.R.X[i] - sys->vals.levelset[i]) * (sys->mono_reinit.mat.R.X[i] - sys->vals.levelset[i]);
			// Update levelset value
			sys->vals.levelset[i] = sys->mono_reinit.mat.R.X[i];
		}
		printf("Reinitialization Step: %d\n", step);
		printf("Error: %f, Delta: %f\n", sqrt(error)/sys->vals.dt_reinit, sys->vals.delta_reinit);
		if(sqrt(error)/sys->vals.dt_reinit < sys->vals.delta_reinit){
			break;
		}
	}
}

/**********************************************************
 * Reinitialization of Conservative Levelset
 **********************************************************/
void set_element_mat_CLSM_reinitialize(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double* local_levelset;
	local_levelset = BB_std_calloc_1d_double(local_levelset, nl);
	double* levelset_ip;
	levelset_ip = BB_std_calloc_1d_double(levelset_ip, np);

	double** local_grad_phi;
	local_grad_phi = BB_std_calloc_2d_double(local_grad_phi, nl, 3);
	double** grad_phi_ip;
	grad_phi_ip = BB_std_calloc_2d_double(grad_phi_ip, np, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_scalar(local_levelset, fe, vals->levelset, e);
		BBFE_elemmat_set_local_array_vector(local_grad_phi, fe, vals->grad_phi, e, 3);

		for(int p=0; p<np; p++) {
			levelset_ip[p]  = BBFE_std_mapping_scalar(nl, local_levelset, basis->N[p]);
			BBFE_std_mapping_vector3d(grad_phi_ip[p], nl, local_grad_phi, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
					val_ip[p] = BBFE_elemmat_mat_CLSM_reinitialize(
							basis->N[p][i], basis->N[p][j], 
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j],
							levelset_ip[p], grad_phi_ip[p], vals->dt_reinit, vals->epsilon_reinit);
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				monolis_add_scalar_to_sparse_matrix_R(
						monolis, fe->conn[e][i], fe->conn[e][j], 0, 0, integ_val);
			}
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_1d_double(local_levelset, nl);
	BB_std_free_1d_double(levelset_ip, np);

	BB_std_free_2d_double(local_grad_phi, nl, 3);
	BB_std_free_2d_double(grad_phi_ip, np, 3);
}

void set_element_vec_CLSM_reinitialize(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double* local_levelset;
	local_levelset = BB_std_calloc_1d_double(local_levelset, nl);
	double* levelset_ip;
	levelset_ip = BB_std_calloc_1d_double(levelset_ip, np);

	double** local_normal_vec;
	local_normal_vec = BB_std_calloc_2d_double(local_normal_vec, nl, 3);
	double** normal_vec_ip;
	normal_vec_ip = BB_std_calloc_2d_double(normal_vec_ip, np, 3);

	double** grad_phi_ip;
	grad_phi_ip = BB_std_calloc_2d_double(grad_phi_ip, np, 3);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_scalar(local_levelset, fe, vals->levelset, e);
		BBFE_elemmat_set_local_array_vector(local_normal_vec, fe, vals->grad_phi, e, 3);

		for(int p=0; p<np; p++) {
			levelset_ip[p]  = BBFE_std_mapping_scalar(nl, local_levelset, basis->N[p]);
			BBFE_std_mapping_vector3d(normal_vec_ip[p], nl, local_normal_vec, basis->N[p]);
			BBFE_std_mapping_scalar_grad(grad_phi_ip[p], nl, local_levelset, fe->geo[e][p].grad_N);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				double vec[3];
				val_ip[p] = BBFE_elemmat_vec_CLSM_reinitialize(
					basis->N[p][i],
					fe->geo[e][p].grad_N[i],
					levelset_ip[p], normal_vec_ip[p], grad_phi_ip[p], vals->dt_reinit, vals->epsilon_reinit);
			}
			double integ_val = BBFE_std_integ_calc(
					np, val_ip, basis->integ_weight, Jacobian_ip);

			monolis->mat.R.B[ fe->conn[e][i] ] += integ_val;
		}
	}

	BB_std_free_1d_double(val_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_1d_double(local_levelset, nl);
	BB_std_free_1d_double(levelset_ip, np);

	BB_std_free_2d_double(local_normal_vec, nl, 3);
	BB_std_free_2d_double(normal_vec_ip, np, 3);

	//BB_std_free_2d_double(local_grad_phi, nl, 3);
	BB_std_free_2d_double(grad_phi_ip, np, 3);
}


void reinit_CLSM(
		FE_SYSTEM*  sys)
{
	double t = 0.0;
	int step = 0;
	double error = 1e10;
	while (sqrt(error)/sys->vals.dt_reinit > sys->vals.delta_reinit) {
		t += sys->vals.dt_reinit;
		step += 1;
		monolis_clear_mat_value_R(&(sys->mono_reinit));
		monolis_clear_mat_value_rhs_R(&(sys->mono_reinit));
		
		set_element_mat_CLSM_reinitialize(
					&(sys->mono_reinit),
					&(sys->fe),
					&(sys->basis),
					&(sys->vals));
		set_element_vec_CLSM_reinitialize(
					&(sys->mono_reinit),
					&(sys->fe),
					&(sys->basis),
					&(sys->vals));
		BBFE_sys_monowrap_solve(
					&(sys->mono_reinit),
					&(sys->mono_com),
					sys->mono_reinit.mat.R.X,
					MONOLIS_ITER_BICGSTAB,
					MONOLIS_PREC_DIAG,
					sys->vals.mat_max_iter,
					sys->vals.mat_epsilon);

		error = 0.0;
		for(int i=0; i<(sys->fe.total_num_nodes); i++) {
			// Residual
			error += (sys->mono_reinit.mat.R.X[i] - sys->vals.levelset[i]) * (sys->mono_reinit.mat.R.X[i] - sys->vals.levelset[i]);
			// Update levelset value
			sys->vals.levelset[i] = sys->mono_reinit.mat.R.X[i];
		}
		printf("Reinitialization Step: %d\n", step);
		printf("Error: %f, Delta: %f\n", sqrt(error)/sys->vals.dt_reinit, sys->vals.delta_reinit);
		if(sqrt(error)/sys->vals.dt_reinit < sys->vals.delta_reinit){
			break;
		}
	}
}

/**********************************************************
 * Element matrix for L2 projection
 **********************************************************/
void BBFE_elemmat_mat_L2projection(
		MONOLIS*    monolis,
		BBFE_DATA*  fe,
		BBFE_BASIS* basis,
		double      coef,
		int         block_size)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;

					val_ip[p] += coef * basis->N[p][i] * basis->N[p][j];
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				for(int b=0; b<block_size; b++) {
					if(fe->conn[e][i] == fe->conn[e][j]){ // only diagonal because I matrix is multiplied for L2 projection of liner element
						monolis_add_scalar_to_sparse_matrix_R(
								monolis,
								fe->conn[e][i], fe->conn[e][j], b, b,
								integ_val);
					}
				}
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);
}

/**********************************************************
 * Output Dambreak Data
 **********************************************************/
typedef struct {
	double key;
	double value;
}Pair;

int compare(const void *a, const void *b) {
    double diff = ((Pair*)a)->key - ((Pair*)b)->key;
    if (diff < 0) return -1;
    if (diff > 0) return 1;
    return 0;
}

void output_result_dambreak_data(
		FE_SYSTEM* sys,
		const char* directory,
		double time)
{
	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_DAMBREAK);

	FILE* fp;
	fp = BBFE_sys_write_add_fopen(fp, filename, directory);

	int cnt_bottom = 0;
	int cnt_height = 0;
	const double eps  = 1.0e-10;

	for(int i=0; i<sys->fe.total_num_nodes; i++){
		if(abs(sys->fe.x[i][2] - 0.0) < eps && abs(sys->fe.x[i][1] - 0.0) < eps){
			// (y,z)=(0,0)
			cnt_bottom++;
		}
		if(abs(sys->fe.x[i][0] - 0.0) < eps && abs(sys->fe.x[i][1] - 0.0) < eps){
			// (x,y)=(0,0)
			cnt_height++;
		}
	}

	double* x = BB_std_calloc_1d_double(x, cnt_bottom);
	double* phi_x = BB_std_calloc_1d_double(phi_x, cnt_bottom);
	double* z = BB_std_calloc_1d_double(z, cnt_height);
	double* phi_z = BB_std_calloc_1d_double(phi_z, cnt_height);
	Pair *pairs_x = (Pair*)malloc(cnt_bottom * sizeof(Pair));
	Pair *pairs_z = (Pair*)malloc(cnt_height * sizeof(Pair));
    if (pairs_x == NULL || pairs_z == NULL) {
        fprintf(stderr, "メモリ割り当てエラー\n");
        exit(1);
    }
	int id_x = 0;
	int id_z = 0;

	for(int i=0; i<sys->fe.total_num_nodes; i++){
		if(abs(sys->fe.x[i][2] - 0.0) < eps && abs(sys->fe.x[i][1] - 0.0) < eps){
			// (y,z)=(0,0)
			x[id_x] = sys->fe.x[i][0];
			phi_x[id_x] = sys->vals.levelset[i];
			pairs_x[id_x].key = x[id_x];
			pairs_x[id_x].value = phi_x[id_x];
			id_x++;
		}
		if(abs(sys->fe.x[i][0] - 0.0) < eps && abs(sys->fe.x[i][1] - 0.0) < eps){
			// (x,y)=(0,0)
			z[id_z] = sys->fe.x[i][2];
			phi_z[id_z] = sys->vals.levelset[i];
			pairs_z[id_z].key = z[id_z];
			pairs_z[id_z].value = phi_z[id_z];
			id_z++;
		}
	}

	qsort(pairs_x, cnt_bottom, sizeof(Pair), compare);
	qsort(pairs_z, cnt_height, sizeof(Pair), compare);

	if(abs(time-0.0)<eps){
		fprintf(fp, "%s, %s, %s, \n", "Time", "x", "z");
	}
	fprintf(fp, "%lf, ", time);
	for(int i=0; i<cnt_bottom-1; i++){
		//fprintf(fp, "%lf, %lf\n", pairs_x[i].key, pairs_x[i].value);
		if(pairs_x[i].value>0 && pairs_x[i+1].value<0){
			double x1 = pairs_x[i].key;
			double x2 = pairs_x[i+1].key;
			double p1 = pairs_x[i].value;
			double p2 = pairs_x[i+1].value;

			double x_zero = (p1*x2 - p2*x1)/(p1 - p2);
			
			//fprintf(fp, "%lf, %lf\n", x1, p1);
			fprintf(fp, "%lf, ", x_zero);
			//fprintf(fp, "%lf, %lf\n", x2, p2);
		}else if(abs(pairs_x[i].value)<eps){
			fprintf(fp, "%lf, ", pairs_x[i].key);
		}
	}

	for(int i=0; i<cnt_height; i++){
		if(pairs_z[i].value>0 && pairs_z[i+1].value<0){
			double z1 = pairs_z[i].key;
			double z2 = pairs_z[i+1].key;
			double p1 = pairs_z[i].value;
			double p2 = pairs_z[i+1].value;

			double z_zero = (p1*z2 - p2*z1)/(p1 - p2);
			
			//fprintf(fp, "%lf, %lf\n", z1, p1);
			fprintf(fp, "%lf\n", z_zero);
			//fprintf(fp, "%lf, %lf\n", z2, p2);
		}else if(abs(pairs_z[i].value)<eps){
			fprintf(fp, "%lf\n", pairs_z[i].key);
		}
	}

	BB_std_free_1d_double(x, cnt_bottom);
	BB_std_free_1d_double(phi_x, cnt_bottom);
	BB_std_free_1d_double(z, cnt_height);
	BB_std_free_1d_double(phi_z, cnt_height);
	free(pairs_x);
	free(pairs_z);

	fclose(fp);

}

/**********************************************************
 * Output Rising Bubble Data
 **********************************************************/
double calc_data_bubble(
			BBFE_DATA*   fe,
			BBFE_BASIS*  basis,
			VALUES*      vals,
			double*      data)
{
	double mean_vel_z = 0;
	double mean_pos_z = 0;
	double total_vol = 0;

	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double*  vel_ip;
	double*  pos_ip;
	double*  vol_ip;
	double*  Jacobian_ip;
	vel_ip      = BB_std_calloc_1d_double(vel_ip     , np);
	pos_ip      = BB_std_calloc_1d_double(pos_ip     , np);
	vol_ip      = BB_std_calloc_1d_double(vol_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);
	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double** local_x;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	double** x_ip; 
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);

	double* local_heaviside;
	local_heaviside = BB_std_calloc_1d_double(local_heaviside, nl);
	double* heaviside_ip;
	heaviside_ip = BB_std_calloc_1d_double(heaviside_ip, np);

	for(int e=0; e<(fe->total_num_elems); e++) {
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(
				np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);
		BBFE_elemmat_set_local_array_scalar(local_heaviside, fe, vals->heaviside, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);	
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);	
			heaviside_ip[p]  = BBFE_std_mapping_scalar(nl, local_heaviside, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int p=0; p<np; p++) {
				vel_ip[p] = basis->N[p][i] * (0.5 * (0+v_ip[p][2]) + heaviside_ip[p] * (0-v_ip[p][2]));
				pos_ip[p] = basis->N[p][i] * (0.5 * (0+x_ip[p][2]) + heaviside_ip[p] * (0-x_ip[p][2])) ;
				vol_ip[p] = basis->N[p][i] *  (0.5 * (0+1) + heaviside_ip[p] * (0-1));
			}
			double integ_vel = BBFE_std_integ_calc(
					np, vel_ip, basis->integ_weight, Jacobian_ip);

			double integ_pos = BBFE_std_integ_calc(
					np, pos_ip, basis->integ_weight, Jacobian_ip);

			double integ_vol = BBFE_std_integ_calc(
					np, vol_ip, basis->integ_weight, Jacobian_ip);

			mean_vel_z += integ_vel;
			mean_pos_z += integ_pos;
			total_vol += integ_vol;
		}
	}

	data[0] = mean_pos_z / total_vol;
	data[1] = mean_vel_z / total_vol;
	data[2] = total_vol / (0.25*0.25*0.25*4/3*M_PI);
	data[3] = 0;
	
	BB_std_free_1d_double(vel_ip,      np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);

	BB_std_free_1d_double(local_heaviside, nl);
	BB_std_free_1d_double(heaviside_ip, np);
}

void output_result_bubble_data(
			FE_SYSTEM* sys,
			const char* directory,
			double time)
{
	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_BUBBLE);

	FILE* fp;
	fp = BBFE_sys_write_add_fopen(fp, filename, directory);

	double* data = BB_std_calloc_1d_double(data, 4);

	calc_data_bubble(&(sys->fe), &(sys->basis),&(sys->vals), data);

	double eps = 1e-10;
	if(abs(time-0.0)<eps){
		fprintf(fp, "%s, %s, %s, %s, %s\n", "Time", "z", "vz", "sphericity", "size");
	}
	fprintf(fp, "%lf, %lf, %lf, %lf, %lf\n", time, data[0], data[1], data[2], data[3]);
	
	fclose(fp);

	BB_std_free_1d_double(data, 4);

}

int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();
	double t1 = monolis_get_time();

	sys.cond.directory = BBFE_fluid_get_directory_name(argc, argv, CODENAME);	
	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_fluid_pre(
			&(sys.fe), &(sys.basis),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis);

	const char* filename;
	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC_V);
	BBFE_fluid_sups_read_Dirichlet_bc(
			&(sys.bc),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			4);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);
	
	BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));
	
	BBFE_sys_monowrap_init_monomat(&(sys.monolis) , &(sys.mono_com), &(sys.fe), 4, sys.cond.directory);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_levelset)  , &(sys.mono_com), &(sys.fe), 1, sys.cond.directory);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_reinit)  , &(sys.mono_com), &(sys.fe), 1, sys.cond.directory);
	BBFE_sys_monowrap_init_monomat(&(sys.mono_L2)  , &(sys.mono_com), &(sys.fe), 3, sys.cond.directory);

	BBFE_elemmat_mat_L2projection(&(sys.mono_L2), &(sys.fe), &(sys.basis), 1.0, 3);

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_LEVELSET);
	read_levelset_file(&(sys.fe), &(sys.vals), filename, sys.cond.directory);
	//BBFE_fluid_convert_levelset2CLSM(sys.vals.levelset, sys.vals.size_interface, sys.fe.total_num_nodes);
	BBFE_fluid_convert_levelset2heaviside(sys.vals.heaviside, sys.vals.levelset, sys.vals.size_interface, sys.fe.total_num_nodes);

	/****************** solver ********************/
	double t = 0.0;
	int step = 0;
	int file_num = 0;
	output_files(&sys, file_num, t);
	output_result_dambreak_data(&sys, sys.cond.directory, t);
	output_result_bubble_data(&(sys), sys.cond.directory, t);

	while (t < sys.vals.finish_time) {
		t += sys.vals.dt;
		step += 1;

		printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);

		monolis_clear_mat_value_R(&(sys.monolis));
		
		monolis_clear_mat_value_R(&(sys.mono_levelset));
		monolis_clear_mat_value_rhs_R(&(sys.mono_levelset));

		monolis_clear_mat_value_rhs_R(&(sys.mono_L2));

		// clear surface tension vector
		for(int m=0;m<sys.fe.total_num_nodes;m++){
			sys.vals.surf_tension[m][0] = 0;
			sys.vals.surf_tension[m][1] = 0;
			sys.vals.surf_tension[m][2] = 0;
		}

		printf("%s --- update density and viscosity step ---\n", CODENAME);
		BBFE_fluid_convert_levelset2heaviside(sys.vals.heaviside, sys.vals.levelset, sys.vals.size_interface, sys.fe.total_num_nodes);
		//*
		BBFE_fluid_renew_vals_by_levelset(
				sys.vals.heaviside, 
				sys.vals.density,
				sys.vals.density_l,
				sys.vals.density_g,
				sys.fe.total_num_nodes);

		BBFE_fluid_renew_vals_by_levelset(
				sys.vals.heaviside, 
				sys.vals.viscosity,
				sys.vals.viscosity_l,
				sys.vals.viscosity_g,
				sys.fe.total_num_nodes); //*/

		/*
		BBFE_fluid_renew_vals_by_CLSM(
				sys.vals.levelset, 
				sys.vals.density,
				sys.vals.density_l,
				sys.vals.density_g,
				sys.fe.total_num_nodes);

		BBFE_fluid_renew_vals_by_levelset(
				sys.vals.levelset, 
				sys.vals.viscosity,
				sys.vals.viscosity_l,
				sys.vals.viscosity_g,
				sys.fe.total_num_nodes);
		//*/
		
		printf("%s --- directly solve velocity and pressure step ---\n", CODENAME);
		set_element_mat(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		set_element_vec(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));

		BBFE_sys_monowrap_set_Dirichlet_bc(
				&(sys.monolis),
				sys.fe.total_num_nodes,
				4,
				&(sys.bc),
				sys.monolis.mat.R.B);
		BBFE_sys_monowrap_solve(
				&(sys.monolis),
				&(sys.mono_com),
				sys.monolis.mat.R.X,
				MONOLIS_ITER_BICGSTAB,
				MONOLIS_PREC_DIAG,
				sys.vals.mat_max_iter,
				sys.vals.mat_epsilon);
		BBFE_fluid_sups_renew_velocity(
				sys.vals.v, 
				sys.monolis.mat.R.X,
				sys.fe.total_num_nodes);
		
		printf("%s --- levelset function convection step ---\n", CODENAME);
		set_element_mat_levelset(
				&(sys.mono_levelset),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		set_element_vec_levelset(
				&(sys.mono_levelset),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		BBFE_sys_monowrap_solve(
				&(sys.mono_levelset),
				&(sys.mono_com),
				sys.mono_levelset.mat.R.X,
				MONOLIS_ITER_BICGSTAB,
				MONOLIS_PREC_SOR,
				sys.vals.mat_max_iter,
				sys.vals.mat_epsilon);

		BBFE_fluid_renew_levelset(
				sys.vals.levelset, 
				sys.mono_levelset.mat.R.X,
				sys.fe.total_num_nodes);

		printf("%s --- L2 projection of Levelset step ---\n", CODENAME);
		set_element_vec_L2_projection(
				&(sys.mono_L2),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		BBFE_sys_monowrap_solve(
				&(sys.mono_L2),
				&(sys.mono_com),
				sys.mono_L2.mat.R.X,
				MONOLIS_ITER_BICGSTAB,
				MONOLIS_PREC_SOR,
				sys.vals.mat_max_iter,
				sys.vals.mat_epsilon);
		BBFE_fluid_renew_velocity(
				sys.vals.grad_phi,
				sys.mono_L2.mat.R.X,
				sys.fe.total_num_nodes);

		printf("%s --- Reinitialization step ---\n", CODENAME);
		reinit_levelset(&(sys));
		//reinit_CLSM(&(sys));

		/**********************************************/

		if(step%sys.vals.output_interval == 0) {

			BBFE_fluid_sups_renew_pressure(
				sys.vals.p,
				sys.monolis.mat.R.X,
				sys.vals.density,
				sys.fe.total_num_nodes);

			output_files(&sys, file_num+1, t);
			output_result_dambreak_data(&sys, sys.cond.directory, t);
			output_result_bubble_data(&(sys), sys.cond.directory, t);

			file_num += 1;
		}

	}

	BBFE_fluid_finalize(&(sys.fe), &(sys.basis));
	BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc), sys.fe.total_num_nodes, 4);
	monolis_finalize(&(sys.monolis));
	monolis_finalize(&(sys.mono_levelset));
	monolis_finalize(&(sys.mono_reinit));

	double t2 = monolis_get_time();
	int myrank = monolis_mpi_get_global_my_rank();

	if(myrank == 0) {
		printf("** Total time: %f\n", t2 - t1);
	}

	monolis_global_finalize();

	printf("\n");

	return 0;

}
