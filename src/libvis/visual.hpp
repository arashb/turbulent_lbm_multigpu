#include <stdio.h>
#include "visual.h"


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


void write_vtkPointCoordinates( FILE *fp, int imax, int jmax, int kmax,
                      double dx, double dy, double dz)
{
	double originX = 0.0;
	double originY = 0.0;
	double originZ = 0.0;

	int i = 0;
	int j = 0;
	int k = 0;

	for (k = 0; k < kmax + 1; k++) {
		for (j = 0; j < jmax + 1; j++) {
			for (i = 0; i < imax + 1; i++) {
				fprintf(fp, "%f %f %f\n", originX + (i * dx),
						originY + (j * dy), originZ + (k * dz));
			}
		}
	}
}

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
