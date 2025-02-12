#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"
#include "BBFE/sys/FE_dataset.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"

#include <stdio.h>
#include <stdlib.h>

static const char* CODENAME               = "util/vtkddm/vtkddm_hex >";
static const char* VOIDNAME               = "                        ";
static const int   BUFFER_SIZE            = 10000;
static const int   FEM_LOCAL_NUM_NODE_HEX = 8;

static const char* OUTPUT_FILENAME_VTK_LABEL = "num_parted";

static const char* DEF_DIRECTORY            = ".";
static const char* DEF_FILENAME_NODE        = "node.dat";
static const char* DEF_FILENAME_ELEM        = "elem.dat";
static const char* DEF_FILENAME_PART_DIR    = "parted.0";
static const char* DEF_FILENAME_CONN_PAR_ID = "connectivity.dat.id";
static const char* DEF_FILENAME_VTK         = "vtkddm_hex.vtk";

static const char* OPTION_INFILE_NODE        = "-in";
static const char* OPTION_INFILE_ELEM        = "-ie";
static const char* OPTION_INFILE_PART_DIR    = "-id";
static const char* OPTION_INFILE_CONN_PAR_ID = "-im";
static const char* OPTION_OUTFILE_VTK        = "-ov";

static const char* OPTION_DIRECTORY = "-d";

typedef struct
{
    int         total_num_nodes;
    int         total_num_elems;
    double**    node;
    int**       conn;
} FE_DATA;

typedef struct {
    const char* infile_node;
    const char* infile_elem;
    const char* infile_part_dir;
    const char* infile_conn_par_id;
    const char* outfile_vtk;
    const char* directory;
} SETTINGS;

typedef struct
{
    int     num_parallel;
    int*    connectivity_id;
} CONDITION;


int args_manager(
    SETTINGS*   set,
    int         argc,
    char*       argv[])
{
    if(argc < 2) {
        printf("%s Please input parameters.\n", CODENAME);
        printf("%s Format: \n", CODENAME);
        printf("%s    %s [the num. of parallel]\n", 
                VOIDNAME, argv[0]);

        printf("%s Options: \n", CODENAME);
        printf("%s     %s [input/output directory]\n", 
                VOIDNAME, OPTION_DIRECTORY);
        printf("%s     %s [input filename (nodes)]\n", 
                VOIDNAME, OPTION_INFILE_NODE);
        printf("%s     %s [input filename (elements)]\n", 
                VOIDNAME, OPTION_INFILE_ELEM);
        printf("%s     %s [input filename (parted directory name)]\n", 
                VOIDNAME, OPTION_INFILE_PART_DIR);
        printf("%s     %s [input filename (connectivity id)]\n", 
                VOIDNAME, OPTION_INFILE_CONN_PAR_ID);
        printf("%s     %s [output filename (vtk)]\n", 
                VOIDNAME, OPTION_OUTFILE_VTK);

        printf("\n");
        exit(0);
    }

    int num_parallel = atoi(argv[1]);

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
        set->infile_elem = DEF_FILENAME_ELEM;
    }
    else {
        set->infile_elem = argv[num+1];
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_INFILE_PART_DIR);
    if(num == -1) {
        set->infile_part_dir = DEF_FILENAME_PART_DIR;
    }
    else {
        set->infile_part_dir = argv[num+1];
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_INFILE_CONN_PAR_ID);
    if(num == -1) {
        set->infile_conn_par_id = DEF_FILENAME_CONN_PAR_ID;
    }
    else {
        set->infile_conn_par_id = argv[num+1];
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_OUTFILE_VTK);
    if(num == -1) {
        set->outfile_vtk = DEF_FILENAME_VTK;
    }
    else {
        set->outfile_vtk = argv[num+1];
    }

    printf("%s The num. of parallel                   : %d\n", CODENAME, num_parallel);
    printf("%s Input/output directory                 : %s\n", CODENAME, set->directory);
    printf("%s Input filename (nodes)                 : %s\n", CODENAME, set->infile_node);
    printf("%s Input filename (elements)              : %s\n", CODENAME, set->infile_elem);
    printf("%s Input filename (parted directory name) : %s\n", CODENAME, set->infile_part_dir);
    printf("%s Input filename (connectivity id)       : %s\n", CODENAME, set->infile_conn_par_id);
    printf("%s Output filename (vtk)                  : %s\n", CODENAME, set->outfile_vtk);

    return num_parallel;
}


void BB_vtk_cell_data(
    const int       cell_type,
    const int       num_nodes,
    const int       total_num_elems,
    const int       num_local_nodes,
    double**        node,
    int**           conn,
    const char*     filename,
    const char*     directory,
    double*         cell_vals)
{
    FILE* fp;
    fp = BBFE_sys_write_fopen(fp, filename, directory);

    BB_vtk_write_header(fp);
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    BB_vtk_write_points_3d(fp, num_nodes, node);
    BB_vtk_write_cells(fp, total_num_elems, num_local_nodes, conn);
    BB_vtk_write_cell_types(fp, total_num_elems, cell_type);
    fprintf(fp, "CELL_DATA %d\n", total_num_elems);
    BB_vtk_write_elem_vals_scalar(fp, cell_vals, total_num_elems, OUTPUT_FILENAME_VTK_LABEL);

    fclose(fp);
}


void read_fem_mesh_files(
    FE_DATA*    fe,
    SETTINGS*   set)
{
    char f_name[BUFFER_SIZE];
    FILE* fp;

    //節点座標の読み込み
    snprintf(f_name, BUFFER_SIZE, "%s", set->infile_node);
    fp = BBFE_sys_read_fopen(fp, f_name, set->directory);
    fe->total_num_nodes = BB_std_read_file_pointer_get_num_nodes(fp, BUFFER_SIZE);
    fe->node = BB_std_calloc_2d_double(fe->node, fe->total_num_nodes, 3);
    BB_std_read_file_pointer_get_val_2d_double(fp, fe->total_num_nodes, 3, fe->node);

    //コネクティビティの読み込み
    snprintf(f_name, BUFFER_SIZE, "%s", set->infile_elem);
    fp = BBFE_sys_read_fopen(fp, f_name, set->directory);
    fe->total_num_elems = BB_std_read_file_pointer_get_num_elems(fp, FEM_LOCAL_NUM_NODE_HEX, BUFFER_SIZE);
    fe->conn = BB_std_calloc_2d_int(fe->conn, fe->total_num_elems, FEM_LOCAL_NUM_NODE_HEX);
    BB_std_read_file_pointer_get_val_2d_int(fp, fe->total_num_elems, FEM_LOCAL_NUM_NODE_HEX, fe->conn);

    fclose(fp);
}


char* get_input_filename_mpi(
    int             myrank,
    const char*     filename_body)
{
    char ctmp[BUFFER_SIZE];
    snprintf(ctmp, BUFFER_SIZE, "%s.%d", filename_body, myrank);

    char* filename;
    filename = &(ctmp[0]);
    return filename;
}


//MPI connectivity_id.{rank_id} を読込
int read_mpi_id(
    CONDITION*      cond,
    const char*     filename,
    const char*     directory)
{
    FILE* fp;
    char label[BUFFER_SIZE];
    int num_points = 0;
    int each_point = 0;

    fp = BBFE_sys_read_fopen(fp, filename, directory);
    BB_std_scan_line(&fp, BUFFER_SIZE, "%s", label);
    BB_std_scan_line(&fp, BUFFER_SIZE, "%d %d", &num_points, &each_point);
    if( num_points == 0){
        printf("%s File \"%s\", Label \"%s\", num_points = 0.\n", CODENAME, filename, label);
        return 0;
    }
    if( each_point != 1){
        printf("%s File \"%s\", Label \"%s\", each_point != 1.\n", CODENAME, filename, label);
        return 0;
    }

    cond->connectivity_id = BB_std_calloc_1d_int(cond->connectivity_id, num_points);
    for(int i=0; i<num_points; i++){
        BB_std_scan_line(&fp, BUFFER_SIZE, "%d", &(cond->connectivity_id[i]));
    }

    fclose(fp);
    return num_points;
}


void build_num_parted(
    SETTINGS*       set,
    CONDITION*      cond,
    double*         e_global,
    const char*     directory)
{
    char filename_body[BUFFER_SIZE];
    snprintf(filename_body, BUFFER_SIZE, "%s/%s", set->infile_part_dir ,set->infile_conn_par_id);

    char* filename;
    for(int i=0; i<(cond->num_parallel); i++){
        filename = get_input_filename_mpi(i, filename_body);
        int num_elem_par = read_mpi_id(cond, filename, set->directory);

        for(int j=0; j<(num_elem_par); j++){
            int g = cond->connectivity_id[j];
            double num_parted = e_global[g];

            if(num_parted == 0){
                e_global[g] = i + 1;
            }
            else if(num_parted != 0){
                e_global[g] = -1;
            }    
        }
    }
}


int main(
    int     argc,
    char*   argv[])
{
    printf("\n");

    FE_DATA     fe;
    SETTINGS    set;
    CONDITION   cond;

    BB_calc_void();

    cond.num_parallel = args_manager(&set, argc, argv);

    //FEMメッシュデータ読込
    read_fem_mesh_files(&(fe), &(set));

    double* e_global;
    e_global = BB_std_calloc_1d_double(e_global, fe.total_num_elems);

    build_num_parted(&(set), &(cond), e_global, set.directory);

    //vtkファイルへの書出
    BB_vtk_cell_data(
                TYPE_VTK_HEXAHEDRON,
                fe.total_num_nodes,
                fe.total_num_elems,
                FEM_LOCAL_NUM_NODE_HEX,
                fe.node,
                fe.conn,
                set.outfile_vtk,
                set.directory,
                e_global);

    return 0;
}
