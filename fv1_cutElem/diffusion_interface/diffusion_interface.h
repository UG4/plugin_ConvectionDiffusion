/*
 * diffusion_interface.h
 *
 *  Created on: 24.08.2017
 *      Author: suze
 */

#ifndef IMMERSED_INTERFACE_DIFFUSION_H_
#define IMMERSED_INTERFACE_DIFFUSION_H_


#ifdef UG_PARALLEL
 	#include "lib_grid/parallelization/load_balancer_util.h"
#endif

#include "../convection_diffusion_fv1_cutElem.h"
#include "../../convection_diffusion_base.h"
#include "lib_disc/spatial_disc/immersed_util/interface_handler/interface_handler_two_sided_cut/interface_handler_diffusion.h"
#include "loc_to_glob_mapper_diffusion.h"

namespace ug{
namespace ConvectionDiffusionPlugin{



template <	typename TDomain, typename TAlgebra>
class ImmersedInterfaceDiffusion
		: public IImmersedInterface<TDomain, TAlgebra>
{
	public:
	///	world Dimension
		static const int dim = TDomain::dim;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

		typedef typename domain_traits<dim>::grid_base_object grid_base_object;

 		ImmersedInterfaceDiffusion(
 					   SmartPtr<IAssemble<TAlgebra> > ass,
 					   SmartPtr<ConvectionDiffusionPlugin::ConvectionDiffusionFV1_cutElem<TDomain> > spMaster,
 					   SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider,
 					   SmartPtr<CutElementHandler_TwoSided<dim> > cutElementHandler);

 	// destructor
		virtual ~ImmersedInterfaceDiffusion(){};

    //////////////////////////////////////////////////////////////////////////////////////////
    // main methods
    //////////////////////////////////////////////////////////////////////////////////////////

    // general initialisation of set up data;
    // most important: call of
    //      (1) 'initialize_interface()' and
    //      (2) 'update_interface_data()'  to mark the cut elements and interface vertices
        void init(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel,
              const int topLevel, bool bScaleDoFs);
    
    // mainly updates the BoolMarker for elements and vertices
        void update(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                const int baseLevel, const int topLevel, const number time);
    
    //////////////////////////////////////////////////////////////////////////////////////////
    // helper methods for init() and update()
    //////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////
	/// Info - 'initialize_interface()': Main infos: see CutElementHandler
	/// Important for call here:
    ///     --> counts number of cut elements!
	/// 	--> called during init()
	//////////////////////////////////////////////////////////////////////////////////////////

	    const size_t initialize_interface(ConstSmartPtr<DoFDistribution> dd)
            { return m_spCutElementHandler->initialize_interface(dd); }
	    const size_t initialize_interface_Nitsche(vector_type& u, ConstSmartPtr<DoFDistribution> dd);

    //////////////////////////////////////////////////////////////////////////////////////////
    /// methods for writing the source term and boundary conditions on the interface
    //////////////////////////////////////////////////////////////////////////////////////////
    
        void set_source_data_lua(const number interfaceSourceVal)
            { m_spInterfaceHandlerLocal->set_source_data_lua(interfaceSourceVal); }
        void set_jump_data_lua(const number interfaceJumpVal)
            { m_spInterfaceHandlerLocal->set_jump_data_lua(interfaceJumpVal); }
        void set_jump_grad_data_lua(const MathVector<2>& interfaceJumpGradVal)
            { m_spInterfaceHandlerLocal->set_jump_grad_data_lua(interfaceJumpGradVal); }
        void set_diffusion_data_lua(const MathVector<2>& diffusionCoeffs)
            { m_spInterfaceHandlerLocal->set_diffusion_coeff_data_lua(diffusionCoeffs); }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    /// lua-methods for set up:
    //////////////////////////////////////////////////////////////////////////////////////////
    
    // the 'threshold' defines the bandwidth around the immersed interface, in which a node
    //  counts as 'OUTSIDE' or 'ON_INTERFACE' during call of 'CutElementHandler::is_outside()',
    //  'CutElementHandler::is_nearInterface()
        void initialize_threshold(TDomain& domain, const int baseLevel, const int topLevel);
        void set_threshold(size_t level, const number threshold)
            { m_spCutElementHandler->set_threshold(level, threshold); }
        number MeanElementDiameter(TDomain& domain, int level);
    
    // If StdFV-assembling is ON, NO new 'vCornerCoords' will be computed on the cut elements.
    // The original nodes ACCROSS the interface (and ON the euclidian mesh) will be chosen for the
    // computation of the solution of the interface
    //  ==> the standard shape functions as in common 'ficticious domain' methods will be used
    //  ==> the shape functions will NOT be 1 ON the interface and the gradient will NOT
    //          point normal to the interface
        void set_StdFV_assembling(bool bValue) { m_spInterfaceHandlerLocal->set_StdFV_assembling(bValue);}
        bool StdFV_assembling() { return m_spInterfaceHandlerLocal->StdFV_assembling(); }

    
    //////////////////////////////////////////////////////////////////////////////////////////
    /// lua-methods for output
    //////////////////////////////////////////////////////////////////////////////////////////
    
		void set_analytic_solution(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                                   SmartPtr<MultiGrid> mg, const int topLevel);
    
	// adjustment of solution vector in order to compute the error WITHOUT nodes near the interface:
    //   ==> (1) remove additional nodes ON the interface, not lying on the original grid
    //       (2) set solution in nodes lying NEAR the interface, but IN the original grid to zero
         void adjust_for_error(vector_type& u, vector_type& uCopy, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                              SmartPtr<MultiGrid> mg, const int topLevel);

		double compute_solution_value(const MathVector<dim>& vrtPos);
    
    // returns the interface adapted l2 error, computed during 'ConvectionDiffusionFV1_cutElem::add_def_A_elem()':
		number get_L2Error() { return m_spInterfaceHandlerLocal->get_L2Error(); }

    // returns the number of DoFs on the original grid
        size_t get_numDoFs() { return m_spInterfaceMapper->get_numDoFs(); }

    // boolian used for perform the Nitsche-method (i.e. CutElem-method) for the treatment of the immersed interface
		void set_Nitsche(bool bNitsche) { return m_spInterfaceHandlerLocal->set_Nitsche(bNitsche); }

    // setting and getting flag for printing of cut-element data into file:
        void set_print_cutElemData(bool bValue) { m_spInterfaceHandlerLocal->set_print_cutElemData(bValue); }

    // returns the number of DoFs on the original grid
    size_t get_numCutElements(const int gridlevel, const size_t prtIndex)
    {
        ConstSmartPtr<DoFDistribution> dd = m_spApproxSpace->dof_distribution(GridLevel(gridlevel, GridLevel::LEVEL));
        const int levIndex = m_spCutElementHandler->get_Index(gridlevel, dd);
        
        return m_spCutElementHandler->get_numCutElements(levIndex);
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    /// class member
    //////////////////////////////////////////////////////////////////////////////////////////

	private:
	///	current ApproxSpace
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

		SmartPtr<DiffusionInterfaceProvider<dim> > m_spInterfaceProvider;
		SmartPtr<CutElementHandler_TwoSided<dim> > m_spCutElementHandler;
		SmartPtr<InterfaceHandlerLocalDiffusion<dim> > m_spInterfaceHandlerLocal;

		SmartPtr<DiffusionInterfaceMapper<TDomain, TAlgebra> > m_spInterfaceMapper;

};

} // end namespace ConvectionDiffusionPlugin
} // end namespace ug


#include "diffusion_interface_impl.h"

#endif /* IMMERSED_INTERFACE_DIFFUSION_H_ */




