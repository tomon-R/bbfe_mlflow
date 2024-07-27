#pragma once

#include "fluid_core.h"

int compare(
		const void *a, 
		const void *b);

void output_result_dambreak_data(
		BBFE_DATA* fe,
		double* levelset,
		const char* directory,
		double time);

double calc_data_bubble(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		double** v,
		double*  heaviside,
		double*      data);

void output_result_bubble_data(
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		double** v,
		double*  heaviside,
		const char* directory,
		double time);

void output_result_sloshing_data(
		BBFE_DATA* fe,
		double* levelset,
		const char* directory,
		double time);