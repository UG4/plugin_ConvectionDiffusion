/*
 * Copyright (c) 2012-2013:  G-CSC, Goethe University Frankfurt
 * Authors: Daniel Gonzalez
 * Based on the modules by Dmitry Logashenko, Andreas Vogel
 *
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FV1__BND__OUTFLOW_BASE__
#define __H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FV1__BND__OUTFLOW_BASE__
// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../convection_diffusion_base.h"

namespace ug{
namespace ConvectionDiffusionPlugin{


/// \ingroup lib_disc_elem_disc
/// @{

/// The zero-stress (neutral) outflow boundary condition for the incompressible NS equation
/**
 * This class implements the so-called neutral boundary condition for outflow
 * boundaries. Note that this class can be used only with the stabilized
 * vertex-centered discretization of the Navier-Stokes equations.
 *
 * This boundary condition imposes two equations on the unknown functions
 * at the outflow boundary \f$ F \f$:
 * <ul>
 * <li> \f$ \int_F p ds = 0 \f$, and
 * <li> \f$ \sigma \vec{n} \cdot \vec{n} = 0 \f$ on \f$ F \f$.
 * </ul>
 * where \f$ \vec{n} \f$ is the outer normal at \f$ F \f$ and
 * \f$ \sigma = \mu (\nabla \vec{u} + (\nabla \vec{u})^T) \f$ the stress tensor.
 */
template<	typename TDomain>
class ConvectionDiffusionOutflowBase
	: public IElemDisc<TDomain>
{
	protected:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	///	own type
		typedef ConvectionDiffusionOutflowBase<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
    ConvectionDiffusionOutflowBase(SmartPtr< ConvectionDiffusionBase<TDomain> > spMaster);
	
	///	adds a boundary segment
		void add(const char* subsets);
	
	protected:
	///	sets the kinematic viscosity
		virtual void set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > data) = 0;

	public:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return true;}

	protected:
	/// The master discretization:
		SmartPtr< ConvectionDiffusionBase<TDomain> > m_spMaster;
	
	/// The boundary subsets:
		std::vector<std::string> m_vScheduledBndSubSets; // names
		std::vector<int> m_vBndSubSetIndex; // indices
	
		void extract_scheduled_data(); // convert m_vScheduledBndSubSets -> m_vBndSubSetIndex
};

/// @}

} // namespace NavierStokes
} // end namespace ug

#endif /*_H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FV1__BND__OUTFLOW_BASE__*/
