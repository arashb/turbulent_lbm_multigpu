/*
 * VTK_Common.hpp
 *
 *  Created on: May 21, 2013
 *      Author: Arash
 */

#ifndef VTK_COMMON_HPP_
#define VTK_COMMON_HPP_

#include "../libmath/CVector.hpp"

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
template<typename T>
void write_vtkPointCoordinates( FILE *fp, int imax, int jmax, int kmax,
                      T dx, T dy, T dz, CVector<3,int> origin_cell)
{
	T originX = origin_cell[0]*dx;
	T originY = origin_cell[1]*dy;
	T originZ = origin_cell[2]*dz;

	for (int k = 0; k < kmax + 1; k++) {
		for (int j = 0; j < jmax + 1; j++) {
			for (int i = 0; i < imax + 1; i++) {
				fprintf(fp, "%f %f %f\n", originX + (i * dx),
						originY + (j * dy), originZ + (k * dz));
			}
		}
	}
}

#endif /* VTK_COMMON_HPP_ */
