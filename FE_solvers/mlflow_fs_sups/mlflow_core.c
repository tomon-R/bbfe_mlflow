
#include "mlflow_core.h"


const char* BBFE_fluid_get_directory_name(
		int         argc,
		char*       argv[],
		const char* codename)
{
	const char* dir_name;

	if(argc < 2) { dir_name = "."; }
	else         { dir_name = argv[1]; }

	printf("%s Main directory: %s\n", codename, dir_name);

	return dir_name;
}


void BBFE_fluid_pre(
		BBFE_DATA*    fe,
		BBFE_BASIS*   basis,
		int           argc,
		char*         argv[],
		const char*   directory,
		int           num_integ_points_each_axis)
{
	BB_calc_void();

	int n_axis = num_integ_points_each_axis;
	const char* filename;

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_NODE);
	BBFE_sys_read_node(
			fe,
			filename,
			directory);

	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_ELEM);
	BBFE_sys_read_elem(
			fe,
			filename,
			directory,
			n_axis*n_axis*n_axis);

	BBFE_sys_memory_allocation_integ(
			basis,
			n_axis*n_axis*n_axis,
			3);
	BBFE_sys_memory_allocation_shapefunc(
			basis,
			fe->local_num_nodes,
			1,
			n_axis*n_axis*n_axis);

	BBFE_fluid_set_basis(
			basis,
			fe->local_num_nodes,
			n_axis);
}


void BBFE_fluid_set_basis(
		BBFE_BASIS*   basis,
		int           local_num_nodes,
		int           num_integ_points_each_axis)
{
	switch( local_num_nodes ) {
		case 4:
			basis->num_integ_points = 
				BBFE_std_integ_tet_set_arbitrary_points(
						num_integ_points_each_axis,
						basis->integ_point,
						basis->integ_weight);

			for(int i=0; i<(basis->num_integ_points); i++) {
				BBFE_std_shapefunc_tet1st_get_val(
						basis->integ_point[i],
						basis->N[i]);

				BBFE_std_shapefunc_tet1st_get_derivative(
						basis->integ_point[i],
						basis->dN_dxi[i],
						basis->dN_det[i],
						basis->dN_dze[i]);
			}
			printf("%s Element type: 1st-order tetrahedron.\n", CODENAME);
			break;

		case 8:
			basis->num_integ_points = 
				BBFE_std_integ_hex_set_arbitrary_points(
						num_integ_points_each_axis,
						basis->integ_point,
						basis->integ_weight);

			for(int i=0; i<(basis->num_integ_points); i++) {
				BBFE_std_shapefunc_hex1st_get_val(
						basis->integ_point[i],
						basis->N[i]);

				BBFE_std_shapefunc_hex1st_get_derivative(
						basis->integ_point[i],
						basis->dN_dxi[i],
						basis->dN_det[i],
						basis->dN_dze[i]);
			}
			printf("%s Element type: 1st-order hexahedron.\n", CODENAME);
			break;
	}
	printf("%s The number of integration points: %d\n", CODENAME, basis->num_integ_points);
}


void BBFE_fluid_renew_velocity(
		double**  v,
		double*   ans_vec,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		for(int d=0; d<3; d++) {
			v[i][d] = ans_vec[ 3*i + d ];
		}
	}
}

void BBFE_fluid_copy_velocity(
		double**  v_new,
		double**  v_pre,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		for(int d=0; d<3; d++) {
			v_new[i][d] = v_pre[i][d];
		}
	}
}

void BBFE_fluid_finalize(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis)
{
	BBFE_sys_memory_free_integ(basis, 3);
	BBFE_sys_memory_free_shapefunc(basis);

	BBFE_sys_memory_free_node(fe, 3);
	BBFE_sys_memory_free_elem(fe, basis->num_integ_points, 3);
}

void BBFE_fluid_renew_levelset(
		double* v,
		double* ans_vec,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		v[i] = ans_vec[i];
	}
}

void BBFE_fluid_renew_density(
		double* levelset,
		double* density,
		double density_l,
		double density_g,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++){
		density[i] = 0.5 * (density_l + density_g) + levelset[i] * (density_l - density_g);
	}
}

void BBFE_fluid_renew_viscosity(
		double* levelset,
		double* viscosity,
		double viscosity_l,
		double viscosity_g,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++){
		viscosity[i] = 0.5 * (viscosity_l + viscosity_g) + levelset[i] * (viscosity_l - viscosity_g);
	}
}

void BBFE_fluid_convert_levelset2heaviside(
		double* levelset,
		const double mesh_size,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		double h1 = levelset[i]/mesh_size+1/M_PI*sin(M_PI*levelset[i]/mesh_size);
		if(h1 > 1.0) h1 = 1.0;
		if(h1 < -1.0) h1 = -1.0;
		levelset[i] = 0.5 * h1;
	}
}

void BBFE_fluid_sups_renew_velocity(
		double**  v,
		double*   ans_vec,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		for(int d=0; d<3; d++) {
			v[i][d] = ans_vec[ 4*i + d ];
		}
	}
}


void BBFE_fluid_sups_renew_pressure(
		double*   p,
		double*   ans_vec,
		const int total_num_nodes)
{
	for(int i=0; i<total_num_nodes; i++) {
		p[i] = ans_vec[ 4*i + 3 ];
	}
}
