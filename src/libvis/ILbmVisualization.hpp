

#ifndef CVISUALIZER_HPP
#define CVISUALIZER_HPP

#include "../CLbmOpenCl.hpp"
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
	CLbmOpenCl<T> *cLbmOpencl;

public:
    virtual ~ILbmVisualization() {
		if (velocity)
			delete velocity;

		if (density)
			delete density;

		if (flags)
			delete flags;
    };

	virtual void setup(CLbmOpenCl<T> &p_cLbmOpencl) {
		cLbmOpencl = &p_cLbmOpencl;

		CVector<3,int> domain_cells = p_cLbmOpencl.domain_cells;

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
