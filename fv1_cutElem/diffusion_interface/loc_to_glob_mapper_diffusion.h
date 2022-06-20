/*
 * diffusion_interface_handler_local.h
 *
 *  Created on: 19.01.2015
 *      Author: suze
 */

#ifndef DIFFUSION_INTERFACE_MAPPER_H_
#define DIFFUSION_INTERFACE_MAPPER_H_

#include "lib_disc/spatial_disc/immersed_util/immersed_interface_base.h"

namespace ug{
namespace ConvectionDiffusionPlugin{
    

template <typename TDomain, typename TAlgebra>
class DiffusionInterfaceMapper : public IInterfaceMapper<TAlgebra>
{

	public:
 	///	World dimension
		static const int dim = TDomain::dim;
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of geometric base object
		typedef typename domain_traits<dim>::grid_base_object grid_base_object;

        DiffusionInterfaceMapper(){};

        DiffusionInterfaceMapper(SmartPtr<InterfaceHandlerLocalDiffusion<dim> > localHandler)
        	: m_numGridNodes(0),
              m_numNewDoFs(0),
        	  m_resized(false),
        	  m_resized_defect(false),
        	  m_scaleDoFs(false),
              m_spInterfaceHandlerLocal(localHandler)
        {}

    	virtual ~DiffusionInterfaceMapper()	{}

    ///////////////////////////////////////////////////////////////////////////////
    ///
    ///  base class methods not needed for 'DiffusionInterfaceMapper' class
    ///
    ///////////////////////////////////////////////////////////////////////////////

    ///	modifies local solution vector for adapted defect computation
        void modify_LocalData(LocalMatrix& locJ, LocalVector& locU,
                          ConstSmartPtr<DoFDistribution> dd){};
        void modify_LocalData(LocalVectorTimeSeries& uT, LocalMatrix& locJ, LocalVector& locU,
                          ConstSmartPtr<DoFDistribution> dd){};
    
        void modify_LocalData(LocalVector& locD, LocalVector& tmpLocD, LocalVector& locU,
                          ConstSmartPtr<DoFDistribution> dd){};
        void modify_LocalData(LocalVectorTimeSeries& uT, LocalVector& locD, LocalVector& tmpLocD,
                          LocalVector& locU, ConstSmartPtr<DoFDistribution> dd, size_t t){};

    ///////////////////////////////////////////////////////////////////////////////
    ///
    ///  base class methods and helper methods called by them
    ///
    ///////////////////////////////////////////////////////////////////////////////

    	///	base method: send local entries to global rhs
    		void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec,
    				ConstSmartPtr<DoFDistribution> dd);
        /// methods called by base method for cut element case: special mapping due to new DoFs!
    		void add_local_vec_to_global_interface(vector_type& vec, const LocalVector& lvec,
    				ConstSmartPtr<DoFDistribution> dd);
        /// methods called by base method for Nitsche-treatment on cut elements
    		void add_local_vec_to_global_Nitsche(vector_type& vec, const LocalVector& lvec,
    				ConstSmartPtr<DoFDistribution> dd);

       	///	base method: send local entries to global matrix
    		void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat,
    				ConstSmartPtr<DoFDistribution> dd);
        /// methods called by base method for cut element case: special mapping due to new DoFs!
    		void add_local_mat_to_global_interface(matrix_type& mat, const LocalMatrix& lmat,
    				ConstSmartPtr<DoFDistribution> dd);
    		void add_local_mat_to_global_Nitsche(matrix_type& mat, const LocalMatrix& lmat,
    				ConstSmartPtr<DoFDistribution> dd);
    
        /// sets all non-DoFs to identity rows
        ///   --> or diffusion, iff only INSIDE circle-computation only! ==> all other DoFs set to Dirichlet!
    		void adjust_mat_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd);


    ///////////////////////////////////////////////////////////////////////////////
    /// REMARK:
    ///	During DomainDiscretization::assemble_jacobian:
    /// calling
    /// 		--->  m_spAssTuner->modify_GlobalSol(pModifyMemory, vSol, dd);
    ///////////////////////////////////////////////////////////////////////////////

    /// instead of modifying global data: the computed values at the interface DoFs get
    /// written/stored into data of class 'InterfaceHandlerLocalDiffusion::m_verticesValue'
    /// via call of 'set_interface_values()' (for each call of domainDisc, i.e. newton step)
    
        void modify_GlobalSol(vector_type& vecMod, const vector_type& vec,
    				ConstSmartPtr<DoFDistribution> dd);

        void modify_GlobalSol(SmartPtr<VectorTimeSeries<vector_type> > vSolMod,
    				ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
    				ConstSmartPtr<DoFDistribution> dd){};

   ///////////////////////////////////////////////////////////////
   /// new methods:
   ///////////////////////////////////////////////////////////////
    
    /// sets dirichlet rows for non-DoFs (not yet needed here)
        void set_identity_mat_constraint(matrix_type& mat, const LocalMatrix& lmat,
                                     ConstSmartPtr<DoFDistribution> dd);

    // access methods to local data stored in class 'InterfaceHandlerLocalDiffusion'
    //   REMARK: the contribution of the boundary condition on the immersed interface
    //           were assembled within the ElemDisc 'ConvectionDiffusionFV1_cutElem'
    //           and stored in 'InterfaceHandlerLocalDiffusion' to access it here for
    //           performing the local-to-global mapping
        LocalMatrix& get_local_jacobian(const ReferenceObjectID roid)
            { return m_spInterfaceHandlerLocal->get_local_jacobian(roid); }
        LocalVector& get_local_defect(const ReferenceObjectID roid)
            { return m_spInterfaceHandlerLocal->get_local_defect(roid); }
 
    
    
    // writes the solution in the interface nodes, i.e. the NEW DoFs,
    //  into data of 'InterfaceHandlerLocalDiffusion'
		void set_interface_values(const std::vector<double > verticesValues)
		{ m_spInterfaceHandlerLocal->set_interface_values(verticesValues); }

    // called during 'InterfaceHandlerLocalDiffusion::init()
		void set_numDoFs(const size_t numDoFs)          { m_numGridNodes = numDoFs;}
        size_t get_numDoFs()                            { return m_numGridNodes;}
        void set_numNewDoFs(const size_t numNewDoFs)    { m_numNewDoFs = numNewDoFs;}
    
    // returns the value, by which the global matrix and vector needs to be resized
    //          REMARK: for a jump in the solution, 2 new DoFs will be placed for each interface node
        const size_t get_resize_measure(const size_t numDoFs);

		void set_bScaleDoFs(bool bScaleDoF) { m_scaleDoFs = bScaleDoF; }

	private:
		size_t m_numGridNodes;      // number of grid nodes, i.e. of DoFs WITHOUT interface nodes
		size_t m_numNewDoFs;        // number of interface nodes, i.e. additional DoFs
		bool m_resized;
		bool m_resized_defect;
        bool m_scaleDoFs;           // if m_scaleDoFs = true: 2 new DoFs will be placed for
                                    //  each interface node (for jump in value)

        SmartPtr<InterfaceHandlerLocalDiffusion<dim> > m_spInterfaceHandlerLocal;

};

}// namespace ConvectionDiffusionPlugin
} // end ug namespace

#include "loc_to_glob_mapper_diffusion_impl.h"


#endif /* DIFFUSION_INTERFACE_MAPPER_H_ */
