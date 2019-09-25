/*
 * Copyright (c) 2013-2018:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
 * Based on the modules by Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FRACT_FV1__
#define __H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FRACT_FV1__

// ug4 headers
#include "lib_grid/algorithms/deg_layer_mngr.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape_interface.h"

// plugin's internal headers
#include "../convection_diffusion_base.h"
#include "../convection_diffusion_sss.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

/// Discretization for the Convection-Diffusion Equation in fractures
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the convection diffusion equation in the fractures. This
 * module is complementary to convection_diffusion_fv1.
 * The PDE has the form
 * \f[
 * 	\partial_t (m1*c + m2) - \nabla \left( D \nabla c - \vec{v} c - \vec{F} \right )
 * 		+ r1 \cdot c + r2 = f + f2
 * \f]
 * along the fracture. We assume the corresponding fluxes also in the orthogonal
 * direction. Here:
 * <ul>
 * <li>	\f$ c \f$ is the unknown solution
 * <li>	\f$ m1 \equiv m(\vec{x},t) \f$ is the Mass Scaling Term
 * <li>	\f$ m2 \equiv m(\vec{x},t) \f$ is the Mass Term
 * <li>	\f$ D \equiv D(\vec{x},t) \f$ is the Diffusion Tensor
 * <li>	\f$ v \equiv \vec{v}(\vec{x},t) \f$ is the Velocity Field
 * <li>	\f$ F \equiv \vec{F}(\vec{x},t) \f$ is the Flux
 * <li>	\f$ r1 \equiv r(\vec{x},t) \f$ is the Reaction Rate
 * <li>	\f$ r2 \equiv r(\vec{x},t) \f$ is a Reaction Term
 * <li>	\f$ f \equiv f(\vec{x},t) \f$ is a Source Term
 * <li> \f$ f2 \equiv f_2(\vec{x},t) \f$ is a Vector Source Term
 * </ul>
 *
 * \tparam	TDomain		Domain
 */
template<typename TDomain>
class ConvectionDiffusionFractFV1 : public ConvectionDiffusionBase<TDomain>
{
	///	Own type
		typedef ConvectionDiffusionFractFV1<TDomain> this_type;

	public:
	
	///	Base class type
		typedef ConvectionDiffusionBase<TDomain> base_type;

	///	domain type
		typedef typename base_type::domain_type domain_type;

	///	position type
		typedef typename base_type::position_type position_type;
	
	///	World ('full') dimension
		static const int dim = base_type::dim;
		
	///	Manifold ('low') dimension
		static const int low_dim = dim - 1;

	///	fracture manager type
		typedef DegeneratedLayerManager<dim> fract_manager_type;
		
	private:
	
	///	FV geometry for the fractures
		typedef DimFV1Geometry<low_dim, dim> TFractFVGeom;
	
	///	type of the sides of elements (als low-dimensional fracture elements)
		typedef typename fract_manager_type::side_type side_type;
	
	///	max. number of corners of non-degenerated sides
		static const size_t maxFractSideCorners = fract_manager_type::maxLayerSideCorners;
	
	///	abbreviation for local function: brine mass fraction
		static const size_t _C_ = 0;
	
	/// convection shapes type (for the upwind)
		typedef IConvectionShapes<dim> conv_shape_type;
	
	public:
	///	Constructor
		ConvectionDiffusionFractFV1(const char* functions, const char* subsets);

	///	sets the fracture manager
	/**
	 * This function sets the fracture manager object to recognize the fractures.
	 */
		void set_fract_manager (SmartPtr<fract_manager_type> fract_manager) {m_spFractManager = fract_manager;}
	
	///	set the upwind method along the fracture
	/**
	 * This method sets the upwind method used to upwind the convection
	 * along the fracture.
	 *
	 * \param	shapes		upwind method
	 */
		void set_upwind (SmartPtr<conv_shape_type> shapes) {m_spConvShape = shapes;}
	
	//	set Specific parameters for low-dimensional subdomains (fractures):

		void set_aperture(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imAperture.set_data(user);
		}
		void set_aperture(number val)
		{
			set_aperture(make_sp(new ConstUserNumber<dim>(val)));
		}
	#ifdef UG_FOR_LUA
		void set_aperture(const char* fctName)
		{
			set_aperture(LuaUserDataFactory<number, dim>::create(fctName));
		}
		void set_aperture(LuaFunctionHandle fct)
		{
			set_aperture(make_sp(new LuaUserData<number,dim>(fct)));
		}
	#endif

		void set_ortho_velocity(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imOrthoVelocity.set_data(user);
		}
		void set_ortho_velocity(number val)
		{
			set_ortho_velocity(make_sp(new ConstUserNumber<dim>(val)));
		}
    #ifdef UG_FOR_LUA
		void set_ortho_velocity(const char* fctName)
		{
			set_ortho_velocity(LuaUserDataFactory<number,dim>::create(fctName));
		}
		void set_ortho_velocity(LuaFunctionHandle fct)
		{
			set_ortho_velocity(make_sp(new LuaUserData<number,dim>(fct)));
		}
    #endif

		void set_ortho_diffusion(SmartPtr<CplUserData<number, dim> > user)
		{
			m_imOrthoDiffusion.set_data(user);
		}
		void set_ortho_diffusion(number val)
		{
			set_ortho_diffusion(make_sp(new ConstUserNumber<dim>(val)));
		}
    #ifdef UG_FOR_LUA
		void set_ortho_diffusion(const char* fctName)
		{
			set_ortho_diffusion(LuaUserDataFactory<number,dim>::create(fctName));
		}
		void set_ortho_diffusion(LuaFunctionHandle fct)
		{
			set_ortho_diffusion(make_sp(new LuaUserData<number,dim>(fct)));
		}
    #endif

	/// set singular sources and sinks
		void set_sss_manager(SmartPtr<CDSingularSourcesAndSinks<dim> > sss_mngr) {m_sss_mngr = sss_mngr;}

	/// get singular sources and sinks
		SmartPtr<CDSingularSourcesAndSinks<dim> > sss_manager() {return m_sss_mngr;}

	///	hanging nodes are not allowed in this discretization
		virtual bool use_hanging() const {return false;}
	
	private:
	
	/// registers the interface function
		void register_all_funcs ();
		
	///	auxiliary class for registering functions
		struct RegisterLocalDiscr
		{
			this_type * m_pThis;
		
			RegisterLocalDiscr(this_type * pThis) : m_pThis(pThis) {}
		
			template<typename TElem> void operator () (TElem &) {m_pThis->register_func<TElem> ();}
		};

	///	registers the interface function for one element type
		template <typename TElem> void register_func ();
	
	private:
	
//---- General interface assembling functions: ----

	///	check the grid and the type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template <typename TElem>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	finishes the loop over all elements
		template <typename TElem>
		void fsh_elem_loop ();

	///	prepares the element for assembling
		template <typename TElem>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const position_type vCornerCoords[]);

	///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	///	assembles the stiffness part of the local defect
		template <typename TElem>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	///	assembles the stiffness part of the local defect explicit reaction, reaction_rate and source
		template <typename TElem>
		void add_def_A_expl_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	///	assembles the mass part of the local defect
		template <typename TElem>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const position_type vCornerCoords[]);

	///	assembles the local right hand side
		template <typename TElem>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const position_type vCornerCoords[]);

	///	computes the linearized defect w.r.t. to the fracture velocity
		template <typename TElem>
		void lin_def_velocity(const LocalVector& u,
		                      std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                      const size_t nip);

	///	computes the linearized defect w.r.t. to the orthogonal velocity
		template <typename TElem>
		void lin_def_ortho_velocity(const LocalVector& u,
		                      std::vector<std::vector<number> > vvvLinDef[],
		                      const size_t nip);

	///	computes the linearized defect w.r.t. to the fracture diffusion
		template <typename TElem>
		void lin_def_diffusion(const LocalVector& u,
		                       std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
		                       const size_t nip);

	///	computes the linearized defect w.r.t. to the orthogonal diffusion
		template <typename TElem>
		void lin_def_ortho_diffusion(const LocalVector& u,
		                       std::vector<std::vector<number> > vvvLinDef[],
		                       const size_t nip);

	///	computes the linearized defect w.r.t. to the fracture flux
		template <typename TElem>
		void lin_def_flux(const LocalVector& u,
		                  std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                  const size_t nip);

	///	computes the linearized defect w.r.t. to the orthogonal flux
		template <typename TElem>
		void lin_def_ortho_flux(const LocalVector& u,
		                  std::vector<std::vector<number> > vvvLinDef[],
		                  const size_t nip);

	///	computes the linearized defect w.r.t. to the reaction
		template <typename TElem>
		void lin_def_reaction(const LocalVector& u,
		                      std::vector<std::vector<number> > vvvLinDef[],
		                      const size_t nip);

	///	computes the linearized defect w.r.t. to the reaction
		template <typename TElem>
		void lin_def_reaction_rate(const LocalVector& u,
		                           std::vector<std::vector<number> > vvvLinDef[],
		                           const size_t nip);

	///	computes the linearized defect w.r.t to the source term
		template <typename TElem>
		void lin_def_source(const LocalVector& u,
		                    std::vector<std::vector<number> > vvvLinDef[],
		                    const size_t nip);

	///	computes the linearized defect w.r.t to the vector source term
		template <typename TElem>
		void lin_def_vector_source(const LocalVector& u,
		                           std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                           const size_t nip);

	///	computes the linearized defect w.r.t to the vector source term
		template <typename TElem>
		void lin_def_ortho_vector_source(const LocalVector& u,
		                           std::vector<std::vector<number> > vvvLinDef[],
		                           const size_t nip);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem>
		void lin_def_mass_scale(const LocalVector& u,
		                        std::vector<std::vector<number> > vvvLinDef[],
		                        const size_t nip);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem>
		void lin_def_mass(const LocalVector& u,
		                  std::vector<std::vector<number> > vvvLinDef[],
		                  const size_t nip);

//---- Assembling functions for fractures: ----

		template <typename TElem>
		inline void fract_add_jac_A_elem(LocalMatrix& J, const LocalVector& u, TElem* elem, const position_type vCornerCoords[]);

		template <typename TElem>
		inline void fract_bulk_add_jac_A_elem(LocalMatrix& J, const LocalVector& u, TElem* elem, const position_type vCornerCoords[]);

		template <typename TElem>
		inline void fract_add_def_A_elem(LocalVector& d, const LocalVector& u, TElem* elem, const position_type vCornerCoords[]);

		template <typename TElem>
		inline void fract_bulk_add_def_A_elem(LocalVector& d, const LocalVector& u, TElem* elem, const position_type vCornerCoords[]);

		template <typename TElem>
		inline void fract_add_def_A_expl_elem(LocalVector& d, const LocalVector& u, TElem* elem, const position_type vCornerCoords[]);

		template <typename TElem>
		inline void fract_bulk_add_def_A_expl_elem(LocalVector& d, const LocalVector& u, TElem* elem, const position_type vCornerCoords[]);

		template <typename TElem>
		inline void fract_add_rhs_elem(LocalVector& b, TElem* elem, const position_type vCornerCoords[]);

		template <typename TElem>
		inline void fract_bulk_add_rhs_elem(LocalVector& b, TElem* elem, const position_type vCornerCoords[]);

		template <typename TElem>
		inline void add_sss_def_elem(LocalVector& d, const LocalVector& u, TElem* pElem, size_t co, number flux);
	
		template <typename TElem>
		inline void add_sss_jac_elem(LocalMatrix& J, const LocalVector& u, TElem* pElem, size_t co, number flux);
	
	private:
	///	returns the updated convection shapes
		const conv_shape_type& get_updated_conv_shapes (bool computeDeriv);
		
	protected:
	
//	Standard parameters of the convection-diffusion plugin
		using base_type::m_imDiffusion;
		using base_type::m_imVelocity;
		using base_type::m_imFlux;
		using base_type::m_imSource;
		using base_type::m_imSourceExpl;
		using base_type::m_imVectorSource;
		using base_type::m_imReactionRate;
		using base_type::m_imReactionRateExpl;
		using base_type::m_imReaction;
		using base_type::m_imReactionExpl;
		using base_type::m_imMassScale;
		using base_type::m_imMass;

		using base_type::m_exGrad;
		using base_type::m_exValue;
	
//	Specific parameters for low-dimensional subdomains (fractures):
		DataImport<number, dim> m_imAperture; ///< the fracture width (constant per element)
		DataImport<number, dim> m_imOrthoDiffusion; ///< diffusion through the fracture-bulk interface
		DataImport<number, dim> m_imOrthoVelocity; ///< convection velocity through the fracture-bulk interface
		DataImport<number, dim> m_imOrthoFlux; ///< flux through the fracture-bulk interface
		DataImport<number, dim> m_imOrthoVectorSource; ///< vector through the fracture-bulk interface

    private:

// 	singular sources and sinks manager
		SmartPtr<CDSingularSourcesAndSinks<dim> > m_sss_mngr;		
		
	protected:
	
//	Parameters of the discretization:
		SmartPtr<fract_manager_type> m_spFractManager; ///< degenerated fracture manager (may be SPNULL)
		SmartPtr<conv_shape_type> m_spConvShape; ///< method to compute the upwind shapes

	private:
	
//	Temporary data used in the assembling
	TFractFVGeom * m_pFractGeo; ///< FV geometry object of fracture elements
	size_t m_numFractCo; ///< number of corners of the fracture side
	side_type * m_innerFractSide; ///< inner side of the fracture element
	side_type * m_outerFractSide; ///< outer side of the fracture element
	size_t m_innerFractSideIdx; ///< index of the inner side in the ref. elem.
	size_t m_outerFractSideIdx; ///< index of the outer side in the ref. elem.
	size_t m_innerSideCo[maxFractSideCorners]; ///< inner side corner idx -> elem. corner idx
	size_t m_outerSideCo[maxFractSideCorners]; ///< outer side corner idx -> elem. corner idx
	size_t m_assCo [2 * maxFractSideCorners]; ///< correspondence of the corners of the sides
};

} // end ConvectionDiffusionPlugin
} // end namespace ug

#endif /*__H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FRACT_FV1__*/

/* End of File */
