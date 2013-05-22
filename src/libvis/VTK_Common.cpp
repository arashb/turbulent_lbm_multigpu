
#include <stdio.h>
#include "VTK_Common.hpp"

/**
 * Method for writing header information in vtk format.
 * @author Arash Bakhtiari
 */
void write_vtkHeader( FILE *fp )
{
    if( fp == NULL )
    {
        fprintf(stderr, "Null pointer in write_vtkHeader" );
        return;
    }

    fprintf(fp,"# vtk DataFile Version 3.1\n");
    fprintf(fp,"Turbulent Fluid Simulation on MultiGPU.\n");
    fprintf(fp,"ASCII\n");
    fprintf(fp,"\n");
}


/**
 * Append the values of a double variable to a VTK file
 *
 * @param out_file_name target output file
 * @param var_name name of the variable
 * @param start_index index of the first element to include in the output file
 * @param end_index index of the last element to include in the output file
 * @param values array with the values of var_name at the grid points
 */
//void vtk_append_scalar_double(const char *out_file_name, const char *var_name, int start_index,
//                       int end_index, double *values) {
//    int i;
//    FILE *fp = NULL;
//
//    if ( (fp = fopen(out_file_name, "a")) == NULL ) {
//        fprintf(stderr, "Failed to open %s", out_file_name);
//        return;
//    }
//
//    /*
//     * The first line gives the name of the dataset (variable) and its type (here "DOUBLE")
//     * The second line selects the color table to use and it is usually "LOOKUP_TABLE default".
//     * The following lines contain a value of the dataset per line
//     */
//    fprintf(fp, "SCALARS %s DOUBLE\n", var_name);
//    fprintf(fp, "LOOKUP_TABLE default\n");
//    for ( i = start_index; i <= end_index; i++ )
//        fprintf(fp, "%f\n", values[i]);
//
//    fprintf(fp, "\n");
//
//    if ( fclose(fp) ) fprintf(stderr, "Failed to close %s", out_file_name);
//}
//
/**
 * Append the values of an integer variable to a VTK file
 *
 * @param out_file_name target output file
 * @param var_name name of the variable
 * @param start_index index of the first element to include in the output file
 * @param end_index index of the last element to include in the output file
 * @param values array with the values of var_name at the grid points
 */
//void vtk_append_scalar_integer(const char *out_file_name, const char *var_name, int start_index,
//                        int end_index, int *values) {
//    int i;
//    FILE *fp = NULL;
//
//    fp = fopen(out_file_name, "a");
//    if ( fp == NULL ) {
//        fprintf(stderr, "Failed to open %s", out_file_name);
//        return;
//    }
//
//    /*
//     * The first line gives the name of the dataset (variable) and its type (here "INT")
//     * The second line selects the color table to use and it is usually "LOOKUP_TABLE default".
//     * The following lines contain a value of the dataset per line
//     */
//    fprintf(fp, "SCALARS %s INT\n", var_name);
//    fprintf(fp, "LOOKUP_TABLE default\n");
//    for ( i = start_index; i <= end_index; i++ )
//        fprintf(fp, "%d\n", values[i]);
//
//    fprintf(fp, "\n");
//
//    if ( fclose(fp) ) fprintf(stderr, "Failed to close %s", out_file_name);
//}

//void vtk_append_vector_float(FILE *fp, const char *var_name, int start_index,
//        int end_index, float **values) {
//
//	fprintf(fp, "VECTORS %s float\n", var_name);
//	for(j = 0; j < jmax+1; j++) {
//		for(i = 0; i < imax+1; i++) {
//		fprintf(fp, "%f %f 0\n", (U[i][j] + U[i][j+1]) * 0.5, (V[i][j] + V[i+1][j]) * 0.5 );
//		}
//	}
//}

