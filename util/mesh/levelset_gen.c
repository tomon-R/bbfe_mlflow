#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "BB/std.h"
#include "BB/vtk.h"
#include "BBFE/sys/memory.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"
#include "mesh_core.h"

const double MIN_NUM = 1e-4;
static const char* CODENAME = "levelset_gen >";
static const char* VOIDNAME = "           ";
static int SIMULATION_ID = -1;

void args_manager_levelset(
		SETTINGS*   set,
		int         argc,
		char*       argv[],
		bool        surface_mesh)
{
	if(argc < 2) {
		printf("%s Please input parameters.\n", CODENAME);
		printf("%s Format: \n", CODENAME);
		printf("%s    %s 1: Dambreak, 2: Rising Bubble (2D), 3: Rising  Bubble (3D)\n", 
				VOIDNAME, argv[0]);

		printf("%s Options: \n", CODENAME);
		printf("%s     %s [input/output directory]\n", 
				VOIDNAME, OPTION_DIRECTORY);
		printf("%s     %s [input filename (nodes)]\n", 
				VOIDNAME, OPTION_INFILE_NODE);
		printf("%s     %s [input filename (elements)]\n", 
				VOIDNAME, OPTION_INFILE_ELEM);
		printf("%s     %s [output filename (nodes)]\n", 
				VOIDNAME, OPTION_OUTFILE_NODE);
		printf("%s     %s [output filename (elements)]\n", 
				VOIDNAME, OPTION_OUTFILE_ELEM);
		printf("%s     %s [output filename (vtk)]\n", 
				VOIDNAME, OPTION_OUTFILE_VTK);

		printf("\n");
		exit(0);
	}

	SIMULATION_ID = atoi(argv[1]);

	if(SIMULATION_ID == 1 && argc < 8){
		printf("%s Please input parameters.\n", CODENAME);
		printf("%s Format: \n", CODENAME);
		printf("%s    %s 1 [x_min] [y_min] [z_min] [x_max] [y_max] [z_max]\n", 
				VOIDNAME, argv[0]);
		exit(0);
	}
	set->x_min[0] = atof(argv[2]);
	set->x_min[1] = atof(argv[3]);
	set->x_min[2] = atof(argv[4]);

	set->x_max[0] = atof(argv[5]);
	set->x_max[1] = atof(argv[6]);
	set->x_max[2] = atof(argv[7]);

	int num;
	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_DIRECTORY);
	if(num == -1) {
		set->directory = DEF_DIRECTORY;
	}
	else {
		set->directory = argv[num+1];
	}

	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_INFILE_NODE);
	if(num == -1) {
		set->infile_node = DEF_FILENAME_NODE;
	}
	else {
		set->infile_node = argv[num+1];
	}

	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_INFILE_ELEM);
	if(num == -1) {
		if(surface_mesh) {
			set->infile_elem = DEF_FILENAME_SURF;
		}
		else{
			set->infile_elem = DEF_FILENAME_ELEM;
		}
	}
	else {
		set->infile_elem = argv[num+1];
	}

	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_OUTFILE_NODE);
	if(num == -1) {
		set->outfile_node = DEF_FILENAME_NODE;
	}
	else {
		set->outfile_node = argv[num+1];
	}

	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_OUTFILE_ELEM);
	if(num == -1) {
		if(surface_mesh) {
			set->outfile_elem = DEF_FILENAME_SURF;
		}
		else{
			set->outfile_elem = DEF_FILENAME_ELEM;
		}
	}
	else {
		set->outfile_elem = argv[num+1];
	}

	num = BB_std_read_args_return_char_num(
				argc, argv, OPTION_OUTFILE_VTK);
	if(num == -1) {
		if(surface_mesh) {
			set->outfile_vtk = DEF_FILENAME_SVTK;
		}
		else{
			set->outfile_vtk = DEF_FILENAME_VTK;
		}
	}
	else {
		set->outfile_vtk = argv[num+1];
	}
	
	printf("%s SIMULATION_ID: %e\n", CODENAME, SIMULATION_ID);
	printf("%s Input/output directory    : %s\n", CODENAME, set->directory);
	printf("%s Input filename (nodes)    : %s\n", CODENAME, set->infile_node);
	printf("%s Input filename (elements) : %s\n", CODENAME, set->infile_elem);
	printf("%s Output filename (nodes)   : %s\n", CODENAME, set->outfile_node);
	printf("%s Output filename (elements): %s\n", CODENAME, set->outfile_elem);
	printf("%s Output filename (vtk)     : %s\n", CODENAME, set->outfile_vtk);
	
	printf("\n");
}


void write_levelset_file(
		SETTINGS* set,
		BBFE_DATA*  fe,
		double* d,
		const char* filename,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);
	
	fprintf(fp, "%s\n", "#levelset");
	fprintf(fp, "%d %d\n", fe->total_num_nodes, 1);
	
	switch(SIMULATION_ID){
		case 1: 
			//* dambreak
			for(int i=0; i<(fe->total_num_nodes); i++) {
				//printf("%d\n", i);
				if( ((fabs(fe->x[i][0]-set->x_max[0]) < MIN_NUM) && (fe->x[i][2] <= set->x_max[2])) ||
				    ((fabs(fe->x[i][2]-set->x_max[2]) < MIN_NUM) && (fe->x[i][0] <= set->x_max[0])) ){
					fprintf(fp, "%.15e\n", 0.0);
				}else if(fe->x[i][0] >= set->x_min[0] && fe->x[i][0] <= set->x_max[0] &&
				   fe->x[i][1] >= set->x_min[1] && fe->x[i][1] <= set->x_max[1] &&
				   fe->x[i][2] >= set->x_min[2] && fe->x[i][2] <= set->x_max[2]){
				   	fprintf(fp, "%.15e\n", d[i]);
				}else{
					fprintf(fp, "%.15e\n", -d[i]);
				}
			}
			break;
		case 2:
			//* 3d-cavity 
			for(int i=0; i<(fe->total_num_nodes); i++) {
				//printf("%d\n", i);
				if(fabs(fe->x[i][2]-0.5) < MIN_NUM){
					fprintf(fp, "%.15e\n", 0.0);
				}else if(fe->x[i][2] <= set->x_max[2]){
				   	fprintf(fp, "%.15e\n", d[i]);
				}else{
					fprintf(fp, "%.15e\n", -d[i]);
				}
			}
			break;
		case 3:
			//* 2d-bubble
			for(int i=0; i<(fe->total_num_nodes); i++) {
				double dx = fe->x[i][0] - 0.5;
				double dy = fe->x[i][2] - 0.5;
				double dd = sqrt(dx*dx+dy*dy);

				fprintf(fp, "%.15e\n", dd-0.25);
			}
			break;
		case 4:
			//* 3d-bubble
			for(int i=0; i<(fe->total_num_nodes); i++) {
				double dx = fe->x[i][0] - 0.5;
				double dy = fe->x[i][1] - 0.5;
				double dz = fe->x[i][2] - 0.5;
				double dd = sqrt(dx*dx+dy*dy+dz*dz);

				fprintf(fp, "%.15e\n", dd-0.25);
			}
			break;
		defalut:
			printf("SIMULATION_ID is not from 1 to 4");
	}

	fclose(fp);
}

void calc_distance(
		SETTINGS* set,
		BBFE_DATA* fe,
		double* dist)
{
	int* surface;
	surface = (int*)calloc(fe->total_num_nodes, sizeof(int));

	// identify surface for 1: Dambreak and 2: Cavity 
	if(SIMULATION_ID == 1){
		//* dambreak 
		for(int i=0; i<(fe->total_num_nodes); i++) {
			if( ((fabs(fe->x[i][0] - set->x_max[0]) < MIN_NUM) && (fe->x[i][2] <= set->x_max[2])) ||
			    ((fabs(fe->x[i][2] - set->x_max[2]) < MIN_NUM) && (fe->x[i][0] <= set->x_max[0])) ){
				dist[i] = 0.0;
				surface[i] = 1;
			}else{
				dist[i] = 1e7;
				surface[i] = 0;
			}
		}
	}else if(SIMULATION_ID == 2){
		//* 3d-cavity 
		for(int i=0; i<(fe->total_num_nodes); i++) {
			if(fabs(fe->x[i][2]-0.5) < MIN_NUM){
				dist[i] = 0.0;
				surface[i] = 1;
			}else{
				dist[i] = 1e7;
				surface[i] = 0;
			}
		}
	}else if(SIMULATION_ID == 3){
		//* Rising Bubble (2D)
		//  Do nothing
	}else if(SIMULATION_ID == 4){
		//* Rising Bubble (3D)
		//  Do nothing
	}

	// calc distance for all cases 
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if(surface[i] == 0){
			for(int j=0; j<(fe->total_num_nodes); j++){
				if(surface[j] == 1){
					double x1 = fe->x[i][0];
					double y1 = fe->x[i][1];
					double z1 = fe->x[i][2];
					double x2 = fe->x[j][0];
					double y2 = fe->x[j][1];
					double z2 = fe->x[j][2];
					double d = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
					if(dist[i]>d)dist[i] = d;
				}
			}
			//printf("dist[i] = %f\n", dist[i]); 
		}
	}

	free(surface);
	surface = NULL;
	printf("calc_distance()\n");
}

double* init_dist(
		BBFE_DATA* fe, 
		double* dist)
{
	return (double*)calloc(fe->total_num_nodes, sizeof(double));
}

int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	SETTINGS  set;
	BBFE_DATA fe;

	args_manager_levelset(&set, argc, argv, false);

	open_fe_files(
			&fe, 
			set.infile_node,
			set.infile_elem,
			set.directory);

	double* dist = init_dist(&fe, dist);
	calc_distance(&set, &fe, dist);
	write_levelset_file(&set, &fe, dist, "levelset.dat", set.directory);

	printf("\n");

	free(dist);
	dist = NULL;

	return 0;
}
