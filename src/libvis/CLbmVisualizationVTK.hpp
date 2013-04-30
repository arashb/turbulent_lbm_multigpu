#ifndef CLBMVISUALIZATIONVTK_HPP
#define CLBMVISUALIZATIONVTK_HPP

#include "ILbmVisualization.hpp"
#include "../libmath/CVector.hpp"
#include "visual.hpp"

template <typename T>
class CLbmVisualizationVTK : virtual public ILbmVisualization<T>
{
	std::string  _file_name;
	int _timeStepNumber;

public:
	CLbmVisualizationVTK(): ILbmVisualization<T>()
	{
		_file_name = "OUTPUT";
		_timeStepNumber = -1;
	}

	void setup(CLbmOpenCl<T> &p_cLbmOpencl, std::string file_name, int timeStepNumber = -1) {
		ILbmVisualization<T>::setup(p_cLbmOpencl);
		_file_name = file_name;
		_timeStepNumber = timeStepNumber;
	}

	void render(int increment = -1)
	{
		// using the increment in VTK visualization as the time step value
		_timeStepNumber = increment;
		std::cout << " I am making vtk file: " << _file_name << std::endl;

		// load array from GPU
		this->cLbmOpencl->storeVelocity(this->velocity);
		this->cLbmOpencl->storeDensity(this->density);
		this->cLbmOpencl->storeFlags(this->flags);

		// number of cells in each dimension (don't confuse it with the number of points!)
//		int imax = params->imax;
//		int jmax = params->jmax;
//		int kmax = params->kmax;
//		double dx = params->dx;
//		double dy = params->dy;
//		double dz = params->dz;
//
		// TODO: implement the render function for VTK visualization
		char szFileName[80];
		FILE *fp = NULL;
		sprintf( szFileName, "%s.%i.vtk", _file_name.c_str(), _timeStepNumber );
		fp = fopen( szFileName, "w");
		if( fp == NULL )
		{
			fprintf(stderr, "Failed to open %s", szFileName );
			return;
		}

		write_vtkHeader( fp );
//	    fprintf(fp,"DATASET STRUCTURED_GRID\n");
//	    fprintf(fp,"DIMENSIONS  %i %i %i \n", imax+1, jmax+1, kmax+1);
//	    fprintf(fp,"POINTS %i float\n", (imax+1)*(jmax+1)*(kmax+1) );
//	    fprintf(fp,"\n");
//		write_vtkPointCoordinates(fp, imax, jmax, kmax, dx, dy, dz);
//
//		fprintf(fp,"POINT_DATA %i \n", (imax+1)*(jmax+1)*(kmax+1) );
//
//
//	    fprintf(fp,"\n");
//	    fprintf(fp,"CELL_DATA %i \n", ((imax)*(jmax)*(kmax)) );

		if( fclose(fp) )
		{
			fprintf(stderr, "Failed to close %s", szFileName);
		}

	}
};
#endif
