#ifndef CLBMVISUALIZATIONVTK_HPP
#define CLBMVISUALIZATIONVTK_HPP

#include "ILbmVisualization.hpp"
#include "../libmath/CVector.hpp"

template <typename T>
class CLbmVisualizationVTK : virtual public ILbmVisualization<T>
{
public:
	CLbmVisualizationVTK(): ILbmVisualization<T>()
	{

	}

	void render(int increment = -1)
	{
		// TODO
	}
};
#endif
