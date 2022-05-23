/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_BASE__
#define __H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_BASE__

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"


#include "lib_disc/spatial_disc/user_data/linker/scale_add_linker.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

/// \addtogroup convection_diffusion
/// \{

/// Discretization for the Convection-Diffusion Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the convection diffusion equation.
 * The Equation has the form
 * \f[
 * 	\partial_t (m1*c + m2) - \nabla \cdot \left( D \nabla c - \vec{v} c \right)
 *      + \nabla \cdot F + r1 \cdot c + r2 = f + f2
 * \f]
 * with
 * <ul>
 * <li>	\f$ c \f$ is the unknown solution
 * <li>	\f$ m1 \equiv m(\vec{x},t) \f$ is the Mass Scaling Term
 * <li>	\f$ m2 \equiv m(\vec{x},t) \f$ is the Mass Term
 * <li>	\f$ D \equiv D(\vec{x},t) \f$ is the Diffusion Tensor
 * <li>	\f$ v \equiv \vec{v}(\vec{x},t) \f$ is the Velocity Field
 * <li>	\f$ r1 \equiv r(\vec{x},t) \f$ is the Reaction Rate
 * <li>	\f$ r2 \equiv r(\vec{x},t) \f$ is a Reaction Term
 * <li>	\f$ f \equiv f(\vec{x},t) \f$ is a Source Term
 * <li> \f$ f2 \equiv f_2(\vec{x},t) \f$ is a Vector Source Term
 * </ul>
 *
 * \tparam	TDomain		Domain
 * \tparam	TAlgebra	Algebra
 */
template<typename TDomain>
class ConvectionDiffusionUserDataFactory;

template<	typename TDomain>
class ConvectionDiffusionBase
: public IElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IElemDisc<TDomain> base_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor
		ConvectionDiffusionBase(const char* functions, const char* subsets);
	protected:
		void init_imports();
	public:
	///	sets the diffusion tensor
	/**
	 * This method sets the Diffusion tensor used in computations. If no
	 * Tensor is set, a zero value is assumed.
	 */
	///	\{
		void set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > user);
		void set_diffusion(number val);
#ifdef UG_FOR_LUA
		void set_diffusion(const char* fctName);
		void set_diffusion(LuaFunctionHandle fct);
#endif
	///	\}

	///	sets the velocity field
	/**
	 * This method sets the Velocity field. If no field is provided a zero
	 * value is assumed.
	 */
	/// \{
		void set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_velocity(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_velocity(const char* fctName);
		void set_velocity(LuaFunctionHandle fct);
#endif
	/// \}

	///	sets the flux
	/**
	 * This method sets the Flux. If no field is provided a zero
	 * value is assumed.
	 */
	/// \{
		void set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_flux(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_flux(const char* fctName);
		void set_flux(LuaFunctionHandle fct);
#endif
	/// \}

	///	sets the reaction rate
	/**
	 * This method sets the Reaction Rate. A zero value is assumed as default.
	 */
	///	\{
		void set_reaction_rate(SmartPtr<CplUserData<number, dim> > user);
		void set_reaction_rate(number val);
#ifdef UG_FOR_LUA
		void set_reaction_rate(const char* fctName);
		void set_reaction_rate(LuaFunctionHandle fct);
#endif
	///	\}

	///	sets the reaction
	/**
	 * This method sets the Reaction. A zero value is assumed as default.
	 */
	///	\{
		void set_reaction(SmartPtr<CplUserData<number, dim> > user);
		void set_reaction(number val);
#ifdef UG_FOR_LUA
		void set_reaction(const char* fctName);
		void set_reaction(LuaFunctionHandle fct);
#endif
	///	\}

		void set_reaction_rate_explicit(SmartPtr<CplUserData<number, dim> > user);
		void set_reaction_rate_explicit(number val);
#ifdef UG_FOR_LUA
		void set_reaction_rate_explicit(const char* fctName);
#endif

		void set_reaction_explicit(SmartPtr<CplUserData<number, dim> > user);
		void set_reaction_explicit(number val);
#ifdef UG_FOR_LUA
		void set_reaction_explicit(const char* fctName);
#endif

		void set_source_explicit(SmartPtr<CplUserData<number, dim> > user);
		void set_source_explicit(number val);
#ifdef UG_FOR_LUA
		void set_source_explicit(const char* fctName);
#endif

	///	sets the source / sink term
	/**
	 * This method sets the source/sink value. A zero value is assumed as
	 * default.
	 */
	///	\{
		void set_source(SmartPtr<CplUserData<number, dim> > user);
		void set_source(number val);
#ifdef UG_FOR_LUA
		void set_source(const char* fctName);
		void set_source(LuaFunctionHandle fct);
#endif
	///	\}

	///	sets the vector source term
	/**
	 * This method sets the divergence of the source as an effect of an
	 * external field. A zero value is assumed as default, thus this term is
	 * ignored then.
	 */
	///	\{
		void set_vector_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user);
		void set_vector_source(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_vector_source(const char* fctName);
		void set_vector_source(LuaFunctionHandle fct);
#endif
	///	\}

	///	sets mass scale
	/**
	 * This method sets the mass scale value. The default value is 1.0.
	 */
	///	\{
		void set_mass_scale(SmartPtr<CplUserData<number, dim> > user);
		void set_mass_scale(number val);
#ifdef UG_FOR_LUA
		void set_mass_scale(const char* fctName);
		void set_mass_scale(LuaFunctionHandle fct);
#endif
	///	\}

	///	sets mass
	/**
	 * This method sets the mass value. The default value is 0.0.
	 */
	///	\{
		void set_mass(SmartPtr<CplUserData<number, dim> > user);
		void set_mass(number val);
#ifdef UG_FOR_LUA
		void set_mass(const char* fctName);
		void set_mass(LuaFunctionHandle fct);
#endif
	///	\}

	protected:
	///	Data import for Diffusion
		DataImport<MathMatrix<dim,dim>, dim> m_imDiffusion;

	///	Data import for the Velocity field
		DataImport<MathVector<dim>, dim > m_imVelocity;

	///	Data import for the Flux
		DataImport<MathVector<dim>, dim > m_imFlux;

	///	Data import for the reaction term
		DataImport<number, dim> m_imReactionRate;

	///	Data import for the reaction term
		DataImport<number, dim> m_imReaction;

	///	Data import for the reaction_rate term explicit
		DataImport<number, dim> m_imReactionRateExpl;

	///	Data import for the reaction term explicit
		DataImport<number, dim> m_imReactionExpl;

	///	Data import for the source term explicit
		DataImport<number, dim> m_imSourceExpl;

	///	Data import for the right-hand side (volume)
		DataImport<number, dim> m_imSource;

	///	Data import for the right-hand side (vector)
		DataImport<MathVector<dim>, dim > m_imVectorSource;

	///	Data import for the mass scale
		DataImport<number, dim> m_imMassScale;

	///	Data import for the mass scale
		DataImport<number, dim> m_imMass;

	private:
	///	returns if local time series is needed
		virtual bool requests_local_time_series() {return false;}

	public:
		typedef SmartPtr<CplUserData<number, dim> > NumberExport;
		typedef SmartPtr<CplUserData<MathVector<dim>, dim> > GradExport;

	///	returns the export of the value of associated unknown function
		virtual NumberExport value();

	///	returns the export of the gradient of associated unknown function
		virtual GradExport gradient();

	protected:
	///	Export for the concentration
		SmartPtr<DataExport<number, dim> > m_exValue;

	///	Export for the gradient of concentration
		SmartPtr<DataExport<MathVector<dim>, dim> > m_exGrad;



	public:
		void set_partial_velocity(int mask) {m_partialAssMask_Conv = mask;}
		void set_partial_flux(int mask) {m_partialAssMask_Flux = mask;}
		void set_partial_mass(int mask) {m_partialAssMask_Mass = mask;}

	protected:
		// bit 1: do not assemble contribution of derivative
		// bit 2: do not assemble convection
		int m_partialAssMask_Conv;
		int m_partialAssMask_Flux;
		int m_partialAssMask_Mass;

		friend class ConvectionDiffusionUserDataFactory<TDomain>;
};

// end group convection_diffusion
/// \}

/// This class provides linker for mass and flux.
template<typename TDomain>
class ConvectionDiffusionUserDataFactory
{
public:
	static const int dim = TDomain::dim;

protected:
	typedef ConvectionDiffusionBase<TDomain> TConvectionDiffusion;

	typedef ScaleAddLinker<number, dim, number> TMassLinker;
	typedef ScaleAddLinker<MathVector<dim>, dim, number> TFluxLinker;

public:
	ConvectionDiffusionUserDataFactory(){}

	SmartPtr<TMassLinker> create_mass_linker(SmartPtr<TConvectionDiffusion> spConvDiff) const
	{
		SmartPtr<TMassLinker> mass_linker = make_sp<TMassLinker> (new TMassLinker());

		if (spConvDiff->m_imMass.user_data().valid())
		{ mass_linker->add(1.0, spConvDiff->m_imMass.user_data()); }

		if (spConvDiff->m_imMassScale.user_data().valid())
		{ mass_linker->add(spConvDiff->value(), spConvDiff->m_imMassScale.user_data());}

		return mass_linker;
	}

	SmartPtr<TFluxLinker> create_flux_linker(SmartPtr<TConvectionDiffusion> spConvDiff) const
	{
		SmartPtr<TFluxLinker> flux_linker = make_sp<TFluxLinker> (new TFluxLinker());

		if (spConvDiff->m_imFlux.user_data().valid())
		{ flux_linker->add(1.0, spConvDiff->m_imFlux.user_data()); }

		if (spConvDiff->m_imDiffusion.user_data().valid())
		{
			typedef MathVector<dim> TData;
			typedef MathMatrix<dim,dim> TDataScale;
			typedef ScaleAddLinker<TData, dim, TDataScale> T;
			SmartPtr<T> DgradC = make_sp<T> (new T());
			DgradC->add(spConvDiff->m_imDiffusion.user_data(), spConvDiff->gradient());

			flux_linker->add(1.0, DgradC);
		}

		if (spConvDiff->m_imVelocity.user_data().valid())
		flux_linker->add(spConvDiff->value(), spConvDiff->m_imVelocity.user_data());


		return flux_linker;
	}

};

} // end ConvectionDiffusionPlugin
} // end namespace ug


#endif /*__H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_BASE__*/
