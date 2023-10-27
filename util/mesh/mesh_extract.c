#include <stdio.h>
#include <stdbool.h>

#include "BB/std.h"
#include "BB/vtk.h"
#include "BBFE/sys/memory.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/write.h"
#include "mesh_core.h"

static const char* CODENAME = "mesh_extract >";


int main(
		int   argc,
		char* argv[])
{
	printf("\n");

	SETTINGS  set;
	BBFE_DATA fe_orig;
	BBFE_DATA fe;

	args_manager(&set, argc, argv, false);

	open_fe_files(
			&fe_orig, 
			set.infile_node,
			set.infile_elem,
			set.directory);

	extract_elements(&fe, &fe_orig, set.x_min, set.x_max);
	remove_floating_nodes(&fe, &fe_orig);

	write_node_file(   &fe, set.outfile_node, set.directory);
	write_element_file(&fe, set.outfile_elem, set.directory);

	int cell_type = 0;
	cell_type = get_cell_type_vtk_3d(fe.local_num_nodes);
	write_vtk_shape(&fe, set.outfile_vtk, set.directory, cell_type);

	printf("\n");

	return 0;
}
