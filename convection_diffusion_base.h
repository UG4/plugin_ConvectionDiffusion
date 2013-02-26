/*
 * convection_diffusion.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
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
#include "lib_disc/spatial_disc/user_data/std_data_export.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

/// \ingroup lib_disc_elem_disc
/// @{

/// Discretization for the Convection-Diffusion Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the convection diffusion equation.
 * The Equation has the form
 * \f[
 * 	\partial_t (m1*c + m2) - \nabla \left( D \nabla c - \vec{v} c \right) +
 * 		r1 \cdot c + r2 = f + f2
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
template<	typename TDomain>
class ConvectionDiffusionBase
: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor
		ConvectionDiffusionBase(const char* functions, const char* subsets);

	///	sets the diffusion tensor
	/**
	 * This method sets the Diffusion tensor used in computations. If no
	 * Tensor is set, a zero value is assumed.
	 */
	///	\{
		void set_diffusion(SmartPtr<UserData<MathMatrix<dim, dim>, dim> > user);
		void set_diffusion(number val);
#ifdef UG_FOR_LUA
		void set_diffusion(const char* fctName);
#endif
	///	\}

	///	sets the velocity field
	/**
	 * This method sets the Velocity field. If no field is provided a zero
	 * value is assumed.
	 */
	/// \{
		void set_velocity(SmartPtr<UserData<MathVector<dim>, dim> > user);
		void set_velocity(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_velocity(const char* fctName);
#endif
	/// \}

	///	sets the reaction rate
	/**
	 * This method sets the Reaction Rate. A zero value is assumed as default.
	 */
	///	\{
		void set_reaction_rate(SmartPtr<UserData<number, dim> > user);
		void set_reaction_rate(number val);
#ifdef UG_FOR_LUA
		void set_reaction_rate(const char* fctName);
#endif
	///	\}

	///	sets the reaction
	/**
	 * This method sets the Reaction. A zero value is assumed as default.
	 */
	///	\{
		void set_reaction(SmartPtr<UserData<number, dim> > user);
		void set_reaction(number val);
#ifdef UG_FOR_LUA
		void set_reaction(const char* fctName);
#endif
	///	\}

		void set_reaction_rate_explicit(SmartPtr<UserData<number, dim> > user);
		void set_reaction_rate_explicit(number val);
#ifdef UG_FOR_LUA
		void set_reaction_rate_explicit(const char* fctName);
#endif

		void set_reaction_explicit(SmartPtr<UserData<number, dim> > user);
		void set_reaction_explicit(number val);
#ifdef UG_FOR_LUA
		void set_reaction_explicit(const char* fctName);
#endif

		void set_source_explicit(SmartPtr<UserData<number, dim> > user);
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
		void set_source(SmartPtr<UserData<number, dim> > user);
		void set_source(number val);
#ifdef UG_FOR_LUA
		void set_source(const char* fctName);
#endif
	///	\}

	///	sets the vector source term
	/**
	 * This method sets the divergence of the source as an effect of an
	 * external field. A zero value is assumed as default, thus this term is
	 * ignored then.
	 */
	///	\{
		void set_vector_source(SmartPtr<UserData<MathVector<dim>, dim> > user);
		void set_vector_source(const std::vector<number>& vVel);
#ifdef UG_FOR_LUA
		void set_vector_source(const char* fctName);
#endif
	///	\}

	///	sets mass scale
	/**
	 * This method sets the mass scale value. A value of 1.0 is assumed as
	 * default.
	 */
	///	\{
		void set_mass_scale(SmartPtr<UserData<number, dim> > user);
		void set_mass_scale(number val);
#ifdef UG_FOR_LUA
		void set_mass_scale(const char* fctName);
#endif
	///	\}

	///	sets mass
	/**
	 * This method sets the mass value. A value of 0.0 is assumed as
	 * default.
	 */
	///	\{
		void set_mass(SmartPtr<UserData<number, dim> > user);
		void set_mass(number val);
#ifdef UG_FOR_LUA
		void set_mass(const char* fctName);
#endif
	///	\}

	protected:
	///	Data import for Diffusion
		DataImport<MathMatrix<dim,dim>, dim> m_imDiffusion;

	///	Data import for the Velocity field
		DataImport<MathVector<dim>, dim > m_imVelocity;

	///	Data import for the reaction term
		DataImport<number, dim> m_imReactionRate;

	///	Data import for the reaction term
		DataImport<number, dim> m_imReaction;

	///	Data import for the reaction_rate term explicit
		DataImport<number, dim> m_imReactionRate_explicit;

	///	Data import for the reaction term explicit
		DataImport<number, dim> m_imReaction_explicit;

	///	Data import for the source term explicit
		DataImport<number, dim> m_imSource_explicit;

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
		typedef SmartPtr<UserData<number, dim> > NumberExport;
		typedef SmartPtr<UserData<MathVector<dim>, dim> > GradExport;

	///	returns the export of the value of associated unknown function
		virtual SmartPtr<UserData<number, dim> > value();

	///	returns the export of the gradient of associated unknown function
		virtual SmartPtr<UserData<MathVector<dim>, dim> > gradient();

	protected:
	///	Export for the concentration
		SmartPtr<ValueDataExport<dim> > m_exValue;

	///	Export for the gradient of concentration
		SmartPtr<GradientDataExport<dim> > m_exGrad;
};

/// @}

} // end ConvectionDiffusionPlugin
} // end namespace ug


#endif /*__H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_BASE__*/
