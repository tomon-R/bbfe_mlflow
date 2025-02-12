#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"
#include "BBFE/sys/FE_dataset.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"

#include <stdio.h>
#include <stdlib.h>

static const char* CODENAME               = "util/vtk/txt2vtk_hex_nodalvals >";
static const char* VOIDNAME               = "                                ";
static const int   BUFFER_SIZE            = 10000;
static const int   FEM_LOCAL_NUM_NODE_HEX = 8;

static const char* DEF_DIRECTORY            = ".";
static const char* DEF_FILENAME_NODE        = "node.dat";
static const char* DEF_FILENAME_ELEM        = "elem.dat";
static const char* DEF_FILENAME_PART_DIR    = "parted.0";
static const char* DEF_FILENAME_RESULT_VALS = "result_%d_%06d.txt";
static const char* DEF_FILENAME_VTK         = "result_%06d.vtk";

static const char* OPTION_INFILE_NODE        = "-in";
static const char* OPTION_INFILE_ELEM        = "-ie";
static const char* OPTION_INFILE_PART_DIR    = "-id";
static const char* OPTION_INFILE_RESULT_VALS = "-ir";
static const char* OPTION_OUTFILE_VTK        = "-ov";

static const char* OPTION_DIRECTORY = "-d";


typedef struct
{
    int         total_num_nodes;
    int         total_num_elems;
    double**    node;
    int**       conn;
    double**    result_v;
    double*     result_p;
} FE_DATA;


typedef struct {
    const char* infile_node;
    const char* infile_elem;
    const char* infile_part_dir;
    const char* infile_resultvals_txt;
    const char* outfile_vtk;
    const char* directory;
    int         file_num;
} SETTINGS;


int args_manager(
    SETTINGS*   set,
    int         argc,
    char*       argv[])
{
    if(argc < 3) {
        printf("%s Please input parameters.\n", CODENAME);
        printf("%s Format: \n", CODENAME);
        printf("%s    %s [the num. of parallel] [the num. of file]\n", 
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
                VOIDNAME, OPTION_INFILE_RESULT_VALS);
        printf("%s     %s [output filename (vtk)]\n", 
                VOIDNAME, OPTION_OUTFILE_VTK);

        printf("\n");
        exit(0);
    }

    int num_parallel = atoi(argv[1]);

    set->file_num = atoi(argv[2]);

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
                argc, argv, OPTION_INFILE_RESULT_VALS);
    if(num == -1) {
        set->infile_resultvals_txt = DEF_FILENAME_RESULT_VALS;
    }
    else {
        set->infile_resultvals_txt = argv[num+1];
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
    printf("%s The num. of file                       : %d\n", CODENAME, set->file_num);
    printf("%s Input/output directory                 : %s\n", CODENAME, set->directory);
    printf("%s Input filename (nodes)                 : %s\n", CODENAME, set->infile_node);
    printf("%s Input filename (elements)              : %s\n", CODENAME, set->infile_elem);
    printf("%s Input filename (parted directory name) : %s\n", CODENAME, set->infile_part_dir);
    printf("%s Input filename (result values txt)     : %s\n", CODENAME, set->infile_resultvals_txt);
    printf("%s Output filename (result values vtk)    : %s\n", CODENAME, set->outfile_vtk);

    return num_parallel;
}


void read_FEM_mesh_files(
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


//結果txtファイルの読込
void read_resultvals_txt_fluid(
    FE_DATA*        fe,
    const char*     filename,
    const char*     directory)
{
    FILE* fp;
    int num_vals_par=0;
    int global_id;
    double val[4];

    fp = BBFE_sys_read_fopen(fp, filename, directory);
    BB_std_scan_line(&fp, BUFFER_SIZE, "%d", &num_vals_par);

    if( num_vals_par == 0){
        printf("%s File \"%s\", num_vals_par = 0.\n", CODENAME, filename);
        return;
    }
    for(int i=0; i<num_vals_par; i++){
        BB_std_scan_line(&fp, BUFFER_SIZE, "%d\t%le\t%le\t%le\t%le", &global_id, &val[0], &val[1], &val[2], &val[3]);

        fe->result_v[global_id][0] = val[0];
        fe->result_v[global_id][1] = val[1];
        fe->result_v[global_id][2] = val[2];
        fe->result_p[global_id]    = val[3];
    }

    fclose(fp);
    return;
}


//結果vtkファイルの書き出し
void output_resultvals_vtk_fluid(
    FE_DATA*        fe,
    const char*     filename,
    const char*     directory)
{
    FILE* fp;
    fp = BBFE_sys_write_fopen(fp, filename, directory);

    BB_vtk_write_header(fp);
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
    BB_vtk_write_points_3d(fp, fe->total_num_nodes, fe->node);
    BB_vtk_write_cells(fp, fe->total_num_elems, FEM_LOCAL_NUM_NODE_HEX, fe->conn);
    BB_vtk_write_cell_types(fp, fe->total_num_elems, TYPE_VTK_HEXAHEDRON);

    fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
    BB_vtk_write_point_vals_vector(fp, fe->result_v, fe->total_num_nodes, "Velocity");
    BB_vtk_write_point_vals_scalar(fp, fe->result_p, fe->total_num_nodes, "Pressure");

    fclose(fp);
}


//txt読み込みとvtk書き出し
void read_output_nodalvals(
    FE_DATA*        fe,
    SETTINGS*       set,
    const int       num_parallel)
{
    char filename_read[BUFFER_SIZE];
    char filename_write[BUFFER_SIZE];
    for(int i=0; i<(set->file_num); i++){
        for(int j=0; j<(num_parallel); j++){
            snprintf(filename_read, BUFFER_SIZE, set->infile_resultvals_txt, j, i);
            read_resultvals_txt_fluid(fe, filename_read, set->directory);
        }

        snprintf(filename_write, BUFFER_SIZE, set->outfile_vtk, i);
        output_resultvals_vtk_fluid(fe, filename_write, set->directory);
    }
}


int main(
    int     argc,
    char*   argv[])
{
    printf("\n");

    FE_DATA      fe;
    SETTINGS     set;

    BB_calc_void();

    int num_parallel = args_manager(&set, argc, argv);

    //FEMメッシュデータ読込
    read_FEM_mesh_files(&(fe), &(set));
    fe.result_v = BB_std_calloc_2d_double(fe.result_v, fe.total_num_nodes, 3);
    fe.result_p = BB_std_calloc_1d_double(fe.result_p, fe.total_num_nodes);

    //txtファイル読み込み&vtkファイル書き出し
    read_output_nodalvals(&(fe), &(set), num_parallel);

    BB_std_free_2d_double(fe.result_v, fe.total_num_nodes, 3);
    BB_std_free_1d_double(fe.result_p, fe.total_num_nodes);

    return 0;
}
