/*
 * Copyright 2010 Martin Schreiber
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


/*
 * CFlatTextureLayout.hpp
 *
 *  Created on: Jan 3, 2010
 *      Author: martin
 */

#ifndef CFLATTEXTURELAYOUT_HPP_
#define CFLATTEXTURELAYOUT_HPP_

#include "libmath/CVector.hpp"

/**
 * \brief	create and handle information about flat textures
 */
class CFlatTextureLayout
{
public:
	/**
	 * ft_a_height, ft_a_width:     width and size of flat-texture in direction a
	 * width and size of flat textures
	 */
	int ft_x_height;	///< height of flat texture in view direction x
	int ft_x_width;		///< width of flat texture in view direction x

	int ft_y_height;	///< height of flat texture in view direction y
	int ft_y_width;		///< width of flat texture in view direction y

	int ft_z_height;	///< height of flat texture in view direction z
	int ft_z_width;		///< width of flat texture in view direction z

	int ft_z_elements;	///< ft_z_height*ft_z_width

	int ft_z_mod;	///< number of textures in flat textures aligned in texture direction "u" if view is aligned along x
	int ft_y_mod;	///< number of textures in flat textures aligned in texture direction "u" if view is aligned along y
	int ft_x_mod;	///< number of textures in flat textures aligned in texture direction "u" if view is aligned along z

	/**
	 * initialize flat textures
	 * \param domain_cells	domain cells of volume texture
	 */
	inline void init(CVector<3,int> &domain_cells)
	{
		init(domain_cells[0], domain_cells[1], domain_cells[2]);
	}

	/**
	 * initialize flat texture
	 */
	void init(	int res_x,		///< resolution in x direction
				int res_y,		///< resolution in y direction
				int res_z		///< resolution in z direction
	)
	{
		/*
		 * initialize flat texture measurement data (resolution, etc.)
		 * ft_a_mod: specifies the amount of textures along axis a which are aligned together in the flat texture
		 * ft_a_height, ft_a_width:     width and size of flat-texture
		 */
		int tmp;        // mod in y direction to init variables

		tmp = ft_z_mod = (int)sqrt(res_z);
		while (ft_z_mod*tmp < res_z)            // search for the best fitting texture sizes
		{
				ft_z_mod++;
				if (ft_z_mod*tmp >= res_z)
						break;
				tmp++;
		}
		ft_z_height = res_y * tmp;
		ft_z_width = res_x * ft_z_mod;
		ft_z_elements = ft_z_height*ft_z_width;

		tmp = ft_y_mod = (int)sqrt(res_y);
		while (ft_y_mod*tmp < res_y)
		{
				ft_y_mod++;
				if (ft_y_mod*tmp >= res_y)
						break;
				tmp++;
		}
		ft_y_height = res_z * tmp;
		ft_y_width = res_x * ft_y_mod;

		tmp = ft_x_mod = (int)sqrt(res_x);
		while (ft_x_mod*tmp < res_x)
		{
				ft_x_mod++;
				if (ft_x_mod*tmp >= res_x)
						break;
				tmp++;
		}
		ft_x_height = res_z * tmp;
		ft_x_width = res_y * ft_x_mod;
	}
};

#endif /* CFLATTEXTURELAYOUT_HPP_ */
