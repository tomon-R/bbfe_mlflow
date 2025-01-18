
#include "mlflow_utils.h"

static const char* OUTPUT_FILENAME_DAMBREAK = "result_dambreak.csv";
static const char* OUTPUT_FILENAME_BUBBLE   = "result_bubble.csv";
static const char* OUTPUT_FILENAME_SLOSHING = "result_sloshing.csv";

const int BUFFER_SIZE = 10000;
const double EPS = 1.0e-10;

/**********************************************************
 * Output Dambreak Data
 **********************************************************/
typedef struct {
	double key;
	double value;
}Pair;

int compare(
		const void *a, 
		const void *b) 
{
    double diff = ((Pair*)a)->key - ((Pair*)b)->key;
    if (diff < 0) return -1;
    if (diff > 0) return 1;
    return 0;
}

void output_result_dambreak_data(
		BBFE_DATA* fe,
		MONOLIS_COM* monolis_com,
		double* levelset,
		const char* directory,
		double time)
{
	int myrank = monolis_mpi_get_global_my_rank();
	bool* is_internal_elem = NULL;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem, fe->total_num_elems);
	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	FILE* fp = NULL;
	fp = BBFE_sys_write_add_fopen(fp, OUTPUT_FILENAME_DAMBREAK, directory);

	int cnt_bottom = 0;
	int cnt_height = 0;


	for(int e=0; e<(fe->total_num_elems); e++){
		if (!is_internal_elem[e]) continue;
		
		for(int j=0; j<(fe->local_num_nodes); j++) {
			int i = fe->conn[e][j];
			if(fabs(fe->x[i][2] - 0.0) < EPS && fabs(fe->x[i][1] - 0.0) < EPS){
				// (y,z)=(0,0)
				cnt_bottom++;
			}
			if(fabs(fe->x[i][0] - 0.0) < EPS && fabs(fe->x[i][1] - 0.0) < EPS){
				// (x,y)=(0,0)
				cnt_height++;
			}
		}
	}

	double* x = NULL;
	double* phi_x = NULL;
	double* z = NULL;
	double* phi_z  = NULL;
		x = BB_std_calloc_1d_double(x, cnt_bottom);
	phi_x = BB_std_calloc_1d_double(phi_x, cnt_bottom);
		z = BB_std_calloc_1d_double(z, cnt_height);
	phi_z = BB_std_calloc_1d_double(phi_z, cnt_height);
	
	Pair *pairs_x = (Pair*)malloc(cnt_bottom * sizeof(Pair));
	Pair *pairs_z = (Pair*)malloc(cnt_height * sizeof(Pair));
    if (pairs_x == NULL || pairs_z == NULL) {
        fprintf(stderr, "Error of memory allocation\n");
        exit(1);
    }
	int id_x = 0;
	int id_z = 0;

	for(int e=0; e<(fe->total_num_elems); e++){
		if (!is_internal_elem[e]) continue;
		
		for(int j=0; j<(fe->local_num_nodes); j++) {
			int i = fe->conn[e][j];
			if(fabs(fe->x[i][2] - 0.0) < EPS && fabs(fe->x[i][1] - 0.0) < EPS){
				// (y,z)=(0,0)
				x[id_x] = fe->x[i][0];
				phi_x[id_x] = levelset[i];
				pairs_x[id_x].key = x[id_x];
				pairs_x[id_x].value = phi_x[id_x];
				id_x++;
			}
			if(fabs(fe->x[i][0] - 0.0) < EPS && fabs(fe->x[i][1] - 0.0) < EPS){
				// (x,y)=(0,0)
				z[id_z] = fe->x[i][2];
				phi_z[id_z] = levelset[i];
				pairs_z[id_z].key = z[id_z];
				pairs_z[id_z].value = phi_z[id_z];
				id_z++;
			}
		}
	}

	qsort(pairs_x, cnt_bottom, sizeof(Pair), compare);
	qsort(pairs_z, cnt_height, sizeof(Pair), compare);

	if(myrank == 0){
		if(fabs(time-0.0)<EPS){
			fprintf(fp, "%s, %s, %s, \n", "Time", "x", "z");
		}
		fprintf(fp, "%lf, ", time);
	}

	double x_zero = 0;
	for(int i=0; i<cnt_bottom-1; i++){
		//fprintf(fp, "%lf, %lf\n", pairs_x[i].key, pairs_x[i].value);
		if(pairs_x[i].value>0 && pairs_x[i+1].value<0){
			double x1 = pairs_x[i].key;
			double x2 = pairs_x[i+1].key;
			double p1 = pairs_x[i].value;
			double p2 = pairs_x[i+1].value;

			x_zero = (p1*x2 - p2*x1)/(p1 - p2);
			
			//fprintf(fp, "%lf, %lf\n", x1, p1);
			//fprintf(fp, "%lf, %lf\n", x2, p2);
		}else if(fabs(pairs_x[i].value)<EPS){
			//fprintf(fp, "%lf, ", pairs_x[i].key);
		}
	}
	monolis_allreduce_R(1, &x_zero, MONOLIS_MPI_SUM, monolis_com->comm);
	if(myrank == 0){
		fprintf(fp, "%lf, ", x_zero);
	}

	double z_zero = 0;
	for(int i=0; i<cnt_height; i++){
		if(pairs_z[i].value>0 && pairs_z[i+1].value<0){
			double z1 = pairs_z[i].key;
			double z2 = pairs_z[i+1].key;
			double p1 = pairs_z[i].value;
			double p2 = pairs_z[i+1].value;

			z_zero = (p1*z2 - p2*z1)/(p1 - p2);
			
			//fprintf(fp, "%lf, %lf\n", z1, p1);
			//fprintf(fp, "%lf, %lf\n", z2, p2);
		}else if(fabs(pairs_z[i].value)<EPS){
			//fprintf(fp, "%lf\n", pairs_z[i].key);
		}
	}
	monolis_allreduce_R(1, &z_zero, MONOLIS_MPI_SUM, monolis_com->comm);
	if(myrank == 0){
		fprintf(fp, "%lf\n", z_zero);
	}

	BB_std_free_1d_double(x, cnt_bottom);
	BB_std_free_1d_double(phi_x, cnt_bottom);
	BB_std_free_1d_double(z, cnt_height);
	BB_std_free_1d_double(phi_z, cnt_height);
	free(pairs_x);
	free(pairs_z);
	BB_std_free_1d_bool(is_internal_elem, fe->total_num_elems);

	fclose(fp);

}

/**********************************************************
 * Output Rising Bubble Data
 **********************************************************/
void calc_data_bubble(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		MONOLIS_COM* monolis_com,
		double** v,
		double* heaviside,
		double*      data)
{
	double mean_vel_z = 0;
	double mean_pos_z = 0;
	double total_vol = 0;

	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double*  vel_ip = NULL;
	double*  pos_ip = NULL;
	double*  vol_ip = NULL;
	double*  Jacobian_ip = NULL;
	vel_ip      = BB_std_calloc_1d_double(vel_ip     , np);
	pos_ip      = BB_std_calloc_1d_double(pos_ip     , np);
	vol_ip      = BB_std_calloc_1d_double(vol_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v = NULL;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);
	double** v_ip = NULL; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double** local_x = NULL;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	double** x_ip = NULL; 
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);

	double* local_heaviside = NULL;
	local_heaviside = BB_std_calloc_1d_double(local_heaviside, nl);
	double* heaviside_ip = NULL;
	heaviside_ip = BB_std_calloc_1d_double(heaviside_ip, np);

	bool* is_internal_elem = NULL;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem, fe->total_num_elems);
	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
		if (!is_internal_elem[e]) continue;
		
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		BBFE_elemmat_set_local_array_vector(local_v, fe, v, e, 3);
		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);
		BBFE_elemmat_set_local_array_scalar(local_heaviside, fe, heaviside, e);

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

	monolis_allreduce_R(1, &mean_vel_z, MONOLIS_MPI_SUM, monolis_com->comm);
	monolis_allreduce_R(1, &mean_pos_z, MONOLIS_MPI_SUM, monolis_com->comm);
	monolis_allreduce_R(1, &total_vol, MONOLIS_MPI_SUM, monolis_com->comm);

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
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		MONOLIS_COM* monolis_com,
		double** v,
		double*  heaviside,
		const char* directory,
		double time)
{
	double* data = NULL;
	data = BB_std_calloc_1d_double(data, 4);
	int myrank = monolis_mpi_get_global_my_rank();

	calc_data_bubble(fe, basis, monolis_com, v, heaviside, data);

	if(myrank == 0){
		FILE* fp = NULL;
		fp = BBFE_sys_write_add_fopen(fp, OUTPUT_FILENAME_BUBBLE, directory);
		if(fabs(time-0.0)<EPS){
			fprintf(fp, "%s, %s, %s, %s, %s\n", "Time", "z", "vz", "sphericity", "size");
		}
		fprintf(fp, "%lf, %lf, %lf, %lf, %lf\n", time, data[0], data[1], data[2], data[3]);
		fclose(fp);
	}

	BB_std_free_1d_double(data, 4);

}

/**********************************************************
 * Output Sloshing Data
 **********************************************************/
int count_mlflow_measurement_node(
		BBFE_DATA* fe)
{
	int num_nodes = 0;
	for(int i=0; i<fe->total_num_nodes; i++){
		if(fabs(fe->x[i][0] - 0.0) < EPS && fabs(fe->x[i][1] - 0.0) < EPS){
			// (x,y)=(0,0)
			num_nodes += 1;
		}
	}
	return num_nodes;
}

void set_mlflow_measurement_node(
		BBFE_DATA* fe,
		int* measurement_node_id)
{
    int mnid = 0;
	for(int i=0; i<fe->total_num_nodes; i++){
		if(fabs(fe->x[i][0] - 0.0) < EPS && fabs(fe->x[i][1] - 0.0) < EPS){
			// (x,y)=(0,0)
			measurement_node_id[mnid] = i;
			mnid += 1;
		}
	}
}

void output_result_sloshing_data(
		BBFE_DATA* fe,
		double* levelset,
		const char* directory,
		int* measurement_node_id,
		int  num_nodes,
		double time)
{
	int myrank = monolis_mpi_get_global_my_rank();

	double* z     = NULL;
	double* phi_z = NULL;
	z = BB_std_calloc_1d_double(z, num_nodes);
	phi_z = BB_std_calloc_1d_double(phi_z, num_nodes);
	Pair *pairs_z = (Pair*)malloc(num_nodes * sizeof(Pair));
    if (pairs_z == NULL) {
        printf("Memory allocation error\n");
        exit(1);
    }
	for(int i=0; i<num_nodes; i++){
		int mnid = measurement_node_id[i];
		//printf("i, mnid: %d %d\n", i, mnid);
		// (x,y)=(0,0)
		z[i] = fe->x[mnid][2];
		phi_z[i] = levelset[mnid];
		pairs_z[i].key = z[i];
		pairs_z[i].value = phi_z[i];
	}

	qsort(pairs_z, num_nodes, sizeof(Pair), compare);

	FILE* fp = NULL;
	fp = BBFE_sys_write_add_fopen(fp, OUTPUT_FILENAME_SLOSHING, directory);

	if(myrank == 0){
		if(fabs(time-0.0)<EPS){
			fprintf(fp, "%s, %s, \n", "Time", "z");
		}
	}
	for(int i=0; i<num_nodes-1;i++){
		if(pairs_z[i].value>0 && pairs_z[i+1].value<0){
			double z1 = pairs_z[i].key;
			double z2 = pairs_z[i+1].key;
			double p1 = pairs_z[i].value;
			double p2 = pairs_z[i+1].value;

			double z_zero = (p1*z2 - p2*z1)/(p1 - p2);
			
			fprintf(fp, "%lf, ", time);
			fprintf(fp, "%lf\n", z_zero);

		}else if(fabs(pairs_z[i].value)<EPS){
			//fprintf(fp, "%lf, ", time);
			//fprintf(fp, "%lf\n", pairs_z[i].key);
		}
	}
	BB_std_free_1d_double(z, num_nodes);
	BB_std_free_1d_double(phi_z, num_nodes);
	free(pairs_z);

	fclose(fp);
}