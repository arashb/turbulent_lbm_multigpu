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


#ifndef CLBMSKELETON_HPP
#define CLBMSKELETON_HPP

/*! \file
 *
 * \brief Lattice Boltzmann parameter handler and parametrization class
 * \author Martin Schreiber
 * \version 0.1
 * \date 2009-05-20
 *
 */

#include "libmath/CMath.hpp"
#include "libmath/CVector.hpp"
#include "lib/CError.hpp"
#include "CDomain.hpp"

/**
 * \brief skeleton for lattice boltzmann simulation to handle parameters and to do the parametrization
 */
template <typename T>
class CLbmSkeleton
{
public:
	CError error;		///< error hanlder
	bool debug;			///< true, if debug mode is on

	CDomain<T> domain;
	CVector<3,int> domain_cells;	///< available simulation cells in each dimension
	int domain_cells_count;			///< overall number of simulation cells
	CVector<4, T> d_drivenCavityVelocity;
//protected:
	// input values
	T d_domain_x_length;		/**< domain length in x direction (meters)
								 * the length in y and z direction can be computated by domain_cells because the cells have unit lengths
								 * maybe it's possible to specify the domain length in each dimension by a modification of the
								 * equilibrium and collision operators
								 * (with dimension)
								 */

	CVector<3,T> d_gravitation;	///< gravitation vector (with dimension)
	T d_viscosity;				///< viscosity of fluid (with dimension)
	T mass_exchange_factor;		///< accelleration for mass exchange

	// computed values
	T d_cell_length;			///< length of cell
	T d_timestep;				///< timestep

	T d_reynolds;
	// simulation values (parametrized dimensionless)
//	T viscosity;				///< dimensionless viscosity of fluid
	T tau;						///< dimensionless tau for collision operator
	T inv_tau;					///< inverse tau
	T inv_trt_tau;				///< inverse two time relaxation model tau

	CVector<4, T> drivenCavityVelocity;
	CVector<3,T> gravitation;	///< gravitation vector

	T max_sim_gravitation_length;	///< maximum length of dimensionless gravitation vector to restrict maximum force acting on fluid

	/**
	 * update the values for lbm simulation (parametrization)
	 */
	void updateValues(bool info_output = false)
	{
		d_timestep = (d_cell_length*d_cell_length)*((T)2.0*tau - (T)1.0)/((T)6.0 * d_viscosity * CMath<T>::sqrt(mass_exchange_factor));
		gravitation = d_gravitation*((d_timestep*d_timestep)/d_cell_length);
/*
		std::cout << "=================================" << std::endl;
		std::cout << "tau: " << tau << std::endl;
		std::cout << "gravitation: " << d_gravitation << std::endl;
		std::cout << "=================================" << std::endl;
*/
		/*
		 * limit the gravitation parameter for the simulation
		 * to avoid large velocities and thus an unstable simulation
		 */

		if (gravitation.length() >= max_sim_gravitation_length)
		{
			if (info_output)
			{
				std::cout << "limiting timestep (gravitation: " << gravitation << ")" << std::endl;
			}

			d_timestep = CMath<T>::sqrt((max_sim_gravitation_length * d_cell_length) / d_gravitation.length());
			gravitation = d_gravitation*((d_timestep*d_timestep)/d_cell_length);
			tau = (T)0.5*(d_timestep*d_viscosity*CMath<T>::sqrt(mass_exchange_factor)*(T)6.0)/(d_cell_length*d_cell_length)+(T)0.5;
		}

		if (tau < 0.51 || tau > 2.5)
		{
			error << "tau has to be within the boundary [0.51; 2.5]" << std::endl;
			error << "otherwise the simulation becomes unstable! current value: " << tau << std::endl;
		}

		inv_tau = (T)1.0/tau;
		inv_trt_tau = (T)1.0/((T)0.5 + (T)3.0/((T)16.0*tau - (T)8.0));
	}


	/**
	 * set new gravitation force
	 */
	void setGravitation(CVector<3,T> p_d_gravitation)
	{
		d_gravitation = p_d_gravitation;
		updateValues();
	}


	/**
	 * initialize skeleton
	 */
	void init(	//CVector<3,int> &p_domain_cells,
				//T p_d_domain_x_length,
				//CDomain<T> &domain,
				CVector<3,T> &p_d_gravitation,
				T p_d_viscosity,
				T p_mass_exchange_factor,
				T p_max_sim_gravitation_length = 0.0001,
//				T p_tau = (T)1.0/(T)1.95
//				T p_tau = (T)1.0/(T)1.95
//				T p_tau = (T)1.0/(T)1.5
				T p_tau = (T)0.953575
	    )
	{
		domain_cells = domain.getSize();//p_domain_cells;
		d_domain_x_length = domain.getLength()[0];//p_d_domain_x_length;

		d_gravitation = p_d_gravitation;

		domain_cells_count = domain_cells.elements();


		d_viscosity = p_d_viscosity;
		mass_exchange_factor = p_mass_exchange_factor;

		// compute sidelength of a lattice cell
		d_cell_length = d_domain_x_length / (T)domain_cells[0];

		max_sim_gravitation_length = p_max_sim_gravitation_length;
		tau = p_tau;

		updateValues(true);
		drivenCavityVelocity = d_drivenCavityVelocity * d_timestep;// / d_cell_length ;// / d_domain_x_length;
		d_reynolds = d_domain_x_length*d_drivenCavityVelocity[0] / d_viscosity;

#if DEBUG
		{
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << "dim domain cells: " << domain_cells << std::endl;
			std::cout << "dim domain cells max: " << domain_cells.max() << std::endl;
			std::cout << "dim domain x length: " << d_domain_x_length << std::endl;
			std::cout << "dim cell length: " << d_cell_length << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
			std::cout << "dim viscocity: " << d_viscosity << std::endl;
			std::cout << "dim gravitation: " << d_gravitation << std::endl;
			std::cout << "dim timestep: " << d_timestep << std::endl;
#endif
			std::cout << "dim reynolds number: " << d_reynolds << std::endl;
#if DEBUG
			std::cout << "-------------------------------------------" << std::endl;

//			std::cout << "viscocity: " << viscosity << std::endl;
			std::cout << "tau: " << tau << std::endl;
			std::cout << "inv_trt_tau: " << inv_trt_tau << std::endl;
			std::cout << "inv_tau: " << inv_tau << std::endl;
			std::cout << "gravitation: " << gravitation << std::endl;
			std::cout << "-------------------------------------------" << std::endl;
		}
#endif
	}

	/**
	 * default constructor: initialize only debug variable, use init(...) for further initialization
	 */
	//CLbmSkeleton(CDomain<T> _domain, bool p_debug = false): domain(_domain)
	CLbmSkeleton(CDomain<T> _domain, CVector<4, T> _drivenCavityVelocity) :
			domain(_domain), d_drivenCavityVelocity(_drivenCavityVelocity) {
	}
};

#endif
