#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"

#include <stdio.h>
#include <stdlib.h>

static const char* CODENAME    = "util/vtk/vtk_wireframe_hex >";
static const int   BUFFER_SIZE = 10000;
static const int   DIMENSION   = 3;
static const char* OUTPUT_FILENAME = "wireframe.vtk";

void BBFE_sfem_write_wireframe_vtk(
    int         argc,
    char*       argv[],
    const int   cell_type,
    const int   dimension)
{
    int total_num_nodes;
    int total_num_elems;
    int local_num_nodes;
    double** x_nodes;
    int**    connectivity;

    if(argc != 3){
        printf("Please input 2 filename : node.dat & elem.dat\n");
        return ;
    }

    char filename[BUFFER_SIZE];
    FILE* fp;

    //read node file
    snprintf(filename, BUFFER_SIZE, "%s", argv[1]);
    fp  = fopen(filename,"r");
    if(fp == NULL){
        printf("%s WARNING : we cannot open node file\n", CODENAME);
        return ;
    }
    BB_std_scan_line(&fp, BUFFER_SIZE, "%d", &(total_num_nodes));
    printf("%s Num. Wireframe nodes : %d\n", CODENAME, total_num_nodes);
    x_nodes = BB_std_calloc_2d_double(x_nodes, total_num_nodes, dimension);
    for(int i=0; i<(total_num_nodes); i++) {
        BB_std_scan_line(&fp, BUFFER_SIZE, "%lf %lf %lf",
                        &(x_nodes[i][0]), &(x_nodes[i][1]), &(x_nodes[i][2]));
    }
    fclose(fp);

    //read elem file
    snprintf(filename, BUFFER_SIZE, "%s", argv[2]);
    fp  = fopen(filename,"r");
    if(fp == NULL){
        printf("%s WARNING : we cannot open elem file\n", CODENAME);
        return ;
    }
    BB_std_scan_line(&fp, BUFFER_SIZE, "%d %d",&(total_num_elems), &(local_num_nodes));
    printf("%s Num. Wireframe elems: %d\n", CODENAME, total_num_elems);
    connectivity = BB_std_calloc_2d_int(connectivity, total_num_elems, local_num_nodes);
    for(int e=0; e<(total_num_elems); e++) {
        for(int i=0; i<(local_num_nodes); i++) {
            fscanf(fp, "%d", &(connectivity[e][i]));
        }
    }
    fclose(fp);

    //output wireframe vtk
    snprintf(filename, BUFFER_SIZE, "%s", OUTPUT_FILENAME);
    fp  = fopen(filename,"w");
    if(fp == NULL){
        printf("%s WARNING : we cannot open vtk file\n", CODENAME);
        return ;
    }
    BB_vtk_write_header(fp);
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    BB_vtk_write_points_3d(fp, total_num_nodes, x_nodes);
    BB_vtk_write_cells(fp, total_num_elems, local_num_nodes, connectivity);
    BB_vtk_write_cell_types(fp, total_num_elems, cell_type);

    fclose(fp);

    printf("write wireframe.vtk\n");
}


int main(
    int     argc,
    char*   argv[])
{
    BBFE_sfem_write_wireframe_vtk(argc, argv, TYPE_VTK_HEXAHEDRON, DIMENSION);

    return 0;
}