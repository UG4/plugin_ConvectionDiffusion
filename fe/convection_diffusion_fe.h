/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FE__
#define __H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FE__

// library intern headers
#include "../convection_diffusion_base.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape_interface.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

// \ingroup lib_disc_elem_disc
/// \addtogroup convection_diffusion
/// \{

/// Discretization for the Convection-Diffusion Equation
/**
 * This class implements the IElemDisc interface to provide element local
 * assemblings for the convection diffusion equation.
 * The Equation has the form
 * \f[
 * 	\partial_t (m1*c + m2) - \nabla \left( D \nabla c - \vec{v} c \right + \vec{F}) +
 * 		r1 \cdot c + r2 = f + f2
 * \f]
 * with
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
 * \tparam	TAlgebra	Algebra
 */
template<	typename TDomain>
class ConvectionDiffusionFE : public ConvectionDiffusionBase<TDomain>
{
	private:
	///	Base class type
		typedef ConvectionDiffusionBase<TDomain> base_type;

	///	Own type
		typedef ConvectionDiffusionFE<TDomain> this_type;

	/// error estimator type
		typedef SideAndElemErrEstData<TDomain> err_est_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor
		ConvectionDiffusionFE(const char* functions, const char* subsets);

	///	sets the quad order
		void set_quad_order(size_t order);

	private:
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template <typename TElem, typename TFEGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template <typename TElem, typename TFEGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	///	finishes the loop over all elements
		template <typename TElem, typename TFEGeom>
		void fsh_elem_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem, typename TFEGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem, typename TFEGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the stiffness part of the local defect
		template <typename TElem, typename TFEGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the mass part of the local defect
		template <typename TElem, typename TFEGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local right hand side
		template <typename TElem, typename TFEGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	prepares the loop over all elements of one type for the computation of the error estimator
		template <typename TElem, typename TFEGeom>
		void prep_err_est_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for assembling the error estimator
		template <typename TElem, typename TFEGeom>
		void prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	computes the error estimator contribution for one element
		template <typename TElem, typename TFEGeom>
		void compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);

	///	computes the error estimator contribution for one element
		template <typename TElem, typename TFEGeom>
		void compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);

	///	computes the error estimator contribution for one element
		template <typename TElem, typename TFEGeom>
		void compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);

	///	postprocesses the loop over all elements of one type in the computation of the error estimator
		template <typename TElem, typename TFEGeom>
		void fsh_err_est_elem_loop();

	protected:
	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TFEGeom>
		void lin_def_velocity(const LocalVector& u,
		                      std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                      const size_t nip);

	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TFEGeom>
		void lin_def_diffusion(const LocalVector& u,
		                       std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
		                       const size_t nip);

	///	computes the linearized defect w.r.t to the flux
		template <typename TElem, typename TFEGeom>
		void lin_def_flux(const LocalVector& u,
		                  std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                  const size_t nip);

	///	computes the linearized defect w.r.t to the reaction
		template <typename TElem, typename TFEGeom>
		void lin_def_reaction(const LocalVector& u,
		                      std::vector<std::vector<number> > vvvLinDef[],
		                      const size_t nip);

	///	computes the linearized defect w.r.t to the reaction
		template <typename TElem, typename TFEGeom>
		void lin_def_reaction_rate(const LocalVector& u,
		                           std::vector<std::vector<number> > vvvLinDef[],
		                           const size_t nip);

	///	computes the linearized defect w.r.t to the source term
		template <typename TElem, typename TFEGeom>
		void lin_def_source(const LocalVector& u,
		                    std::vector<std::vector<number> > vvvLinDef[],
		                    const size_t nip);

	///	computes the linearized defect w.r.t to the vector source term
		template <typename TElem, typename TFEGeom>
		void lin_def_vector_source(const LocalVector& u,
		                           std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                           const size_t nip);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem, typename TFEGeom>
		void lin_def_mass_scale(const LocalVector& u,
		                        std::vector<std::vector<number> > vvvLinDef[],
		                        const size_t nip);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem, typename TFEGeom>
		void lin_def_mass(const LocalVector& u,
		                  std::vector<std::vector<number> > vvvLinDef[],
		                  const size_t nip);
							  
	private:
	///	abbreviation for the local solution
		static const size_t _C_ = 0;

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

	protected:
	///	computes the concentration
		template <typename TElem, typename TFEGeom>
		void ex_value(number vValue[],
		              const MathVector<dim> vGlobIP[],
		              number time, int si,
		              const LocalVector& u,
		              GridObject* elem,
		              const MathVector<dim> vCornerCoords[],
		              const MathVector<TFEGeom::dim> vLocIP[],
		              const size_t nip,
		              bool bDeriv,
		              std::vector<std::vector<number> > vvvDeriv[]);

	///	computes the gradient of the concentration
		template <typename TElem, typename TFEGeom>
		void ex_grad(MathVector<dim> vValue[],
		             const MathVector<dim> vGlobIP[],
		             number time, int si,
		             const LocalVector& u,
		             GridObject* elem,
		             const MathVector<dim> vCornerCoords[],
		             const MathVector<TFEGeom::dim> vLocIP[],
		             const size_t nip,
		             bool bDeriv,
		             std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);

	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	returns if hanging nodes are needed
		virtual bool use_hanging() const;

	protected:
	///	current integration order
		bool m_bQuadOrderUserDef;
		int m_quadOrder;

	///	current shape function set
		LFEID m_lfeID;

	///	register utils
	///	\{
		void register_all_funcs(const LFEID& lfeid, const int quadOrder);
		template <typename TElem, typename TFEGeom> void register_func();
	/// \}

	private:
		/// struct holding values of shape functions in IPs
		struct ShapeValues
		{
			public:
				void resize(std::size_t nEip, std::size_t nSip, std::size_t _nSh)
				{
					nSh = _nSh;
					elemVals.resize(nEip);
					sideVals.resize(nSip);
					for (std::size_t i = 0; i < nEip; i++) elemVals[i].resize(nSh);
					for (std::size_t i = 0; i < nSip; i++) sideVals[i].resize(nSh);
				}
				number& shapeAtElemIP(std::size_t sh, std::size_t ip) {return elemVals[ip][sh];}
				number& shapeAtSideIP(std::size_t sh, std::size_t ip) {return sideVals[ip][sh];}
				number* shapesAtElemIP(std::size_t ip) {return &elemVals[ip][0];}
				number* shapesAtSideIP(std::size_t ip) {return &sideVals[ip][0];}
				std::size_t num_sh() {return nSh;}
			private:
				std::size_t nSh;
				std::vector<std::vector<number> > elemVals;
				std::vector<std::vector<number> > sideVals;
		} m_shapeValues;

};

// end group convection_diffusion
/// \}

} // end ConvectionDiffusionPlugin
} // end namespace ug


#endif /*__H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FE__*/
