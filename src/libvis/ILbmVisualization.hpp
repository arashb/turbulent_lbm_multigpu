

#ifndef CVISUALIZER_HPP
#define CVISUALIZER_HPP

#include "../CLbmSolver.hpp"
#include "../libmath/CVector.hpp"

/**
 * \brief Visualization Abstraction Class
 */
template <typename T>
class ILbmVisualization
{
protected:
	T *velocity;
	T *density;
	int *flags;
	CLbmSolver<T> *cLbmOpencl;

public:
	ILbmVisualization()
		: velocity(NULL),
		  density(NULL),
		  flags(NULL)
	{
	}

    virtual ~ILbmVisualization() {
		if (velocity)
			delete[] velocity;

		if (density)
			delete[] density;

		if (flags)
			delete[] flags;
    };

	virtual void setup(CLbmSolver<T> *p_cLbmOpencl) {
		cLbmOpencl = p_cLbmOpencl;

		CVector<3,int> domain_cells = cLbmOpencl->domain_cells;

		delete [] velocity;
		velocity = new T[domain_cells.elements()*3];

		delete [] density;
		density = new T[domain_cells.elements()];

		delete [] flags;
		flags = new int[domain_cells.elements()];

	}

	virtual void render(int increment = -1) = 0;
};

#endif
