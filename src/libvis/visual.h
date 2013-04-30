#ifndef __VISUAL_H__
#define __VISUAL_H__


/**
 * Method for writing header information in vtk format. 
 * @author Arash Bakhtiari
 */
void write_vtkHeader( FILE *fp );

/**
 * Method for writing grid coordinate information in vtk format. 
 * 
 * @param fp      File pointer for writing info.  
 * @param imax    Maximum number of entries (minus 2) in x-direction
 * @param jmax    Maximum number of entries (minus 2) in y-direction
 * @param kmax	  Maximum number of entries (minus 2) in z-direction
 * @param dx      mesh size dx
 * @param dy      mesh size dy
 * @param dz      mesh size dz
 * 
 * @author Arash Bakhtiari
 */
void write_vtkPointCoordinates( FILE *fp, int imax, int jmax, int kmax,
                                double dx, double dy, double dz);

/**
 * Append the values of a double variable to a VTK file
 *
 * @param out_file_name target output file
 * @param var_name name of the variable
 * @param start_index index of the first element to include in the output file
 * @param end_index index of the last element to include in the output file
 * @param values array with the values of var_name at the grid points
 */
//void vtk_append_double(const char *out_file_name, const char *var_name, int start_index,
//                       int end_index, double *values);

/**
 * Append the values of an integer variable to a VTK file
 *
 * @param out_file_name target output file
 * @param var_name name of the variable
 * @param start_index index of the first element to include in the output file
 * @param end_index index of the last element to include in the output file
 * @param values array with the values of var_name at the grid points
 */
//void vtk_append_integer(const char *out_file_name, const char *var_name, int start_index,
//                        int end_index, int *values);
#endif
