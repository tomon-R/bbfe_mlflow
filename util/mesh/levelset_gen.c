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


void write_levelset_file(
		SETTINGS* set,
		BBFE_DATA*  fe,
		double* d,
		const char* filename,
		const char* directory)
{
	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, filename, directory);
	
	fprintf(fp, "%d\n", fe->total_num_nodes);
	
	/* dambreak 
	for(int i=0; i<(fe->total_num_nodes); i++) {
		//printf("%d\n", i);
		if( ((fabs(fe->x[i][0]-0.146) < MIN_NUM) && (fe->x[i][2] <= 0.292)) ||
		    ((fabs(fe->x[i][2]-0.292) < MIN_NUM) && (fe->x[i][0] <=0.146)) ){
			fprintf(fp, "%.15e\n", 0.0);
		}else if(fe->x[i][0] >= set->x_min[0] && fe->x[i][0] <= set->x_max[0] &&
		   fe->x[i][1] >= set->x_min[1] && fe->x[i][1] <= set->x_max[1] &&
		   fe->x[i][2] >= set->x_min[2] && fe->x[i][2] <= set->x_max[2]){
		   	fprintf(fp, "%.15e\n", d[i]);
		}else{
			fprintf(fp, "%.15e\n", -d[i]);
		}
	}
	//*/
	/* 3d-cavity 
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
	//*/

	//* 2d-bubble
	for(int i=0; i<(fe->total_num_nodes); i++) {
		double dx = fe->x[i][0] - 0.5;
		double dy = fe->x[i][2] - 0.5;
		double dd = sqrt(dx*dx+dy*dy);

		fprintf(fp, "%.15e\n", dd-0.25);
	}
	//*/
	/* 3d-bubble
	for(int i=0; i<(fe->total_num_nodes); i++) {
		double dx = fe->x[i][0] - 0.5;
		double dy = fe->x[i][1] - 0.5;
		double dz = fe->x[i][2] - 0.5;
		double dd = sqrt(dx*dx+dy*dy+dz*dz);

		fprintf(fp, "%.15e\n", dd-0.25);
	}
	//*/

	fclose(fp);
}

void calc_distance(BBFE_DATA* fe, double* dist){
	int* surface;
	surface = (int*)calloc(fe->total_num_nodes, sizeof(int));

	/* dambreak 
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if( ((fabs(fe->x[i][0]-0.146) < MIN_NUM) && (fe->x[i][2] <= 0.292)) ||
		    ((fabs(fe->x[i][2]-0.292) < MIN_NUM) && (fe->x[i][0] <=0.45)) ){
			dist[i] = 0.0;
			surface[i] = 1;
		}else{
			dist[i] = 1e7;
			surface[i] = 0;
		}
	}
	//*/
	/* 3d-cavity 
	for(int i=0; i<(fe->total_num_nodes); i++) {
		if(fabs(fe->x[i][2]-0.5) < MIN_NUM){
			dist[i] = 0.0;
			surface[i] = 1;
		}else{
			dist[i] = 1e7;
			surface[i] = 0;
		}
	}
	//*/

	/*
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
			printf("dist[i] = %f\n", dist[i]); 
		}
	}
	//*/

	free(surface);
	surface = NULL;
	printf("calc_distance()\n");
}

double* init_dist(BBFE_DATA* fe, double* dist){
	return (double*)calloc(fe->total_num_nodes, sizeof(double));
}

int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	SETTINGS  set;
	BBFE_DATA fe;

	args_manager(&set, argc, argv, false);

	open_fe_files(
			&fe, 
			set.infile_node,
			set.infile_elem,
			set.directory);

	double* dist = init_dist(&fe, dist);
	calc_distance(&fe, dist);
	write_levelset_file(&set, &fe, dist, "levelset.dat", set.directory);

	printf("\n");

	free(dist);
	dist = NULL;

	return 0;
}
