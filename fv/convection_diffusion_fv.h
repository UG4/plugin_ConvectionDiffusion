
#ifndef __H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FV__
#define __H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FV__

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
class ConvectionDiffusionFV : public ConvectionDiffusionBase<TDomain>
{
	private:
	///	Base class type
		typedef ConvectionDiffusionBase<TDomain> base_type;

	///	Own type
		typedef ConvectionDiffusionFV<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	Constructor
		ConvectionDiffusionFV(const char* functions, const char* subsets);

	private:
	///	prepares the loop over all elements
	/**
	 * This method prepares the loop over all elements. It resizes the Position
	 * array for the corner coordinates and schedules the local ip positions
	 * at the data imports.
	 */
		template <typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	prepares the element for assembling
	/**
	 * This methods prepares an element for the assembling. The Positions of
	 * the Element Corners are read and the Finite Volume Geometry is updated.
	 * The global ip positions are scheduled at the data imports.
	 */
		template <typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	///	finishes the loop over all elements
		template <typename TElem, typename TFVGeom>
		void fsh_elem_loop();

	///	assembles the local stiffness matrix using a finite volume scheme
		template <typename TElem, typename TFVGeom>
		void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local mass matrix using a finite volume scheme
		template <typename TElem, typename TFVGeom>
		void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the stiffness part of the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the mass part of the local defect
		template <typename TElem, typename TFVGeom>
		void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	assembles the local right hand side
		template <typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);
		
	protected:
	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TFVGeom>
		void lin_def_velocity(const LocalVector& u,
		                      std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                      const size_t nip);

	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TFVGeom>
		void lin_def_diffusion(const LocalVector& u,
		                       std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
		                       const size_t nip);

	///	computes the linearized defect w.r.t to the velocity
		template <typename TElem, typename TFVGeom>
		void lin_def_flux(const LocalVector& u,
		                  std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                  const size_t nip);

	///	computes the linearized defect w.r.t to the reaction
		template <typename TElem, typename TFVGeom>
		void lin_def_reaction(const LocalVector& u,
		                      std::vector<std::vector<number> > vvvLinDef[],
		                      const size_t nip);

	///	computes the linearized defect w.r.t to the reaction
		template <typename TElem, typename TFVGeom>
		void lin_def_reaction_rate(const LocalVector& u,
		                           std::vector<std::vector<number> > vvvLinDef[],
		                           const size_t nip);

	///	computes the linearized defect w.r.t to the source term
		template <typename TElem, typename TFVGeom>
		void lin_def_source(const LocalVector& u,
		                    std::vector<std::vector<number> > vvvLinDef[],
		                    const size_t nip);

	///	computes the linearized defect w.r.t to the vector source term
		template <typename TElem, typename TFVGeom>
		void lin_def_vector_source(const LocalVector& u,
		                           std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
		                           const size_t nip);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem, typename TFVGeom>
		void lin_def_mass_scale(const LocalVector& u,
		                        std::vector<std::vector<number> > vvvLinDef[],
		                        const size_t nip);

	///	computes the linearized defect w.r.t to the mass scale term
		template <typename TElem, typename TFVGeom>
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
		template <typename TElem, typename TFVGeom>
		void ex_value(number vValue[],
		              const MathVector<dim> vGlobIP[],
		              number time, int si,
		              const LocalVector& u,
		              GridObject* elem,
		              const MathVector<dim> vCornerCoords[],
		              const MathVector<TFVGeom::dim> vLocIP[],
		              const size_t nip,
		              bool bDeriv,
		              std::vector<std::vector<number> > vvvDeriv[]);

	///	computes the gradient of the concentration
		template <typename TElem, typename TFVGeom>
		void ex_grad(MathVector<dim> vValue[],
		             const MathVector<dim> vGlobIP[],
		             number time, int si,
		             const LocalVector& u,
		             GridObject* elem,
		             const MathVector<dim> vCornerCoords[],
		             const MathVector<TFVGeom::dim> vLocIP[],
		             const size_t nip,
		             bool bDeriv,
		             std::vector<std::vector<MathVector<dim> > > vvvDeriv[]);

	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	///	returns if hanging nodes are needed
		virtual bool use_hanging() const;

	///	sets the quad order
		void set_quad_order(size_t order);

	protected:
	///	current shape function set
		LFEID m_lfeID;

	///	current integration order
		bool m_bQuadOrderUserDef;
		int m_quadOrder;

	///	register utils
	///	\{
		void register_all_funcs(const LFEID& lfeID, const int quadOrder);
		template <typename TElem, typename TFVGeom> void register_func();
	/// \}
};

// end group convection_diffusion
/// \}

} // end ConvectionDiffusionPlugin
} // end namespace ug


#endif /*__H__UG__LIB_DISC__CONVECTION_DIFFUSION__CONVECTION_DIFFUSION_FV__*/
