/*
 * Copyright (c) 2012-2014:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FV1__BND__OUTFLOW_FV1__
#define __H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FV1__BND__OUTFLOW_FV1__


// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"

#include "../../bnd/outflow_base.h"

namespace ug{
namespace ConvectionDiffusionPlugin{


/// \ingroup lib_disc_elem_disc
/// @{

/// Outflow boundary condition for the  Transport equation
/**
 * This class implements the so-called neutral boundary condition for outflow
 * boundaries. Note that this class can be used only with the stabilized
 * vertex-centered discretization of the Navier-Stokes equations.
 *
 * This boundary condition imposes the  outflow Q at boundary \f$ F \f$:
 * <ul>
 * <li> \f$ \int_Q  = c \vec{u} \cdot \vec{ds} \f$,
 * </ul>
 * where \f$ \vec{ds} \f$ is the outer normal at \f$ F \f$, and
 * \f$ c \f$ the concentration.
 */
template<	typename TDomain>
class ConvectionDiffusionOutflowFV1
	: public ConvectionDiffusionOutflowBase<TDomain>
{
	private:
	///	Base class type
		typedef ConvectionDiffusionOutflowBase<TDomain> base_type;

	///	own type
		typedef ConvectionDiffusionOutflowFV1<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor (setting default values)
        ConvectionDiffusionOutflowFV1(SmartPtr< ConvectionDiffusionBase<TDomain> > spMaster);

	protected:
	///	sets the velocity
		virtual void set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
			{m_imVelocity.set_data(data);}


	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	public:
	///	prepares the element loop
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for evaluation
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	///	finishes the element loop
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

	///	adds the stiffness part to the local jacobian
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	adds the stiffness part to the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	public:
	///	dummy implementations
	///	\{
		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]){}
	/// \}

    
	private:

	/// adds the convective part of the local Jacobian of the momentum equation
		template <typename BF>
		inline void convective_flux_Jac
		(
			const size_t ip,
			const BF& bf,
			LocalMatrix& J,
			const LocalVector& u,
            const number& stdValue
		);
	/// adds the convective part of the local defect of the momentum equation
		template <typename BF>
		inline void convective_flux_defect
		(
			const size_t ip,
			const BF& bf,
			LocalVector& d,
			const LocalVector& u,
			const number& stdValue
		);



    
	protected:
	/// abbreviation for pressure
		static const size_t _C_ = 0;

		using base_type::m_spMaster;
		using base_type::m_vBndSubSetIndex;

    ///    Data import for the Velocity field
        DataImport<MathVector<dim>, dim > m_imVelocity;

	/// Boundary integration points of the velocity
		std::vector<MathVector<dim> > m_vLocIP;
		std::vector<MathVector<dim> > m_vGloIP;

	protected:
		void register_all_funcs(bool bHang);
		template<typename TElem, typename TFVGeom>
		void register_func();

};

/// @}

} // namespace ConvectionDiffusionPlugin
} // end namespace ug

#endif /*__H__UG__PLUGINS__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FV1__BND__OUTFLOW_FV1__*/
