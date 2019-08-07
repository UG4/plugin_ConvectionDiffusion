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
        	: m_spInterfaceHandlerLocal(localHandler),
        	  m_numDoFs(0),
        	  m_resized(false),
        	  m_resized_defect(false),
        	  m_scaleDoFs(false)
        {}

    	virtual ~DiffusionInterfaceMapper()	{}


    	///	send local entries to global rhs
    		void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec,
    				ConstSmartPtr<DoFDistribution> dd);
    		void add_local_vec_to_global_interface(vector_type& vec, const LocalVector& lvec,
    				ConstSmartPtr<DoFDistribution> dd);
    		void add_local_vec_to_global_Nitsche(vector_type& vec, const LocalVector& lvec,
    				ConstSmartPtr<DoFDistribution> dd);

       	///	send local entries to global matrix
    		void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat,
    				ConstSmartPtr<DoFDistribution> dd);
    		void add_local_mat_to_global_interface(matrix_type& mat, const LocalMatrix& lmat,
    				ConstSmartPtr<DoFDistribution> dd);
    		void add_local_mat_to_global_Nitsche(matrix_type& mat, const LocalMatrix& lmat,
    				ConstSmartPtr<DoFDistribution> dd);

    		void adjust_mat_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd);

    	///	modifies local solution vector for adapted defect computation
     		void modify_LocalData(LocalMatrix& locJ, LocalVector& locU,
     				ConstSmartPtr<DoFDistribution> dd){};
     		void modify_LocalData(LocalVectorTimeSeries& uT, LocalMatrix& locJ, LocalVector& locU,
     				ConstSmartPtr<DoFDistribution> dd){};

     		void modify_LocalData(LocalVector& locD, LocalVector& tmpLocD, LocalVector& locU,
     				ConstSmartPtr<DoFDistribution> dd){};
     		void modify_LocalData(LocalVectorTimeSeries& uT, LocalVector& locD, LocalVector& tmpLocD, LocalVector& locU,
     				ConstSmartPtr<DoFDistribution> dd, size_t t){};

    	///	modifies global solution vector for adapted defect computation
     		void modify_GlobalSol(vector_type& vecMod, const vector_type& vec,
    				ConstSmartPtr<DoFDistribution> dd);

    		void modify_GlobalSol(SmartPtr<VectorTimeSeries<vector_type> > vSolMod,
    				ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
    				ConstSmartPtr<DoFDistribution> dd){};

   ///////////////////////////////////////////////////////////////
   /// new methods:
   ///////////////////////////////////////////////////////////////
    	void set_identity_mat_constraint(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd);

		LocalMatrix& get_local_jacobian_tri()
		{ return m_spInterfaceHandlerLocal->get_local_jacobian_tri(); }
		LocalMatrix& get_local_jacobian_quad()
		{ return m_spInterfaceHandlerLocal->get_local_jacobian_quad(); }

		LocalVector& get_local_defect_tri()
		{ return m_spInterfaceHandlerLocal->get_local_defect_tri(); }
		LocalVector& get_local_defect_quad()
		{ return m_spInterfaceHandlerLocal->get_local_defect_quad(); }

		void write_solution(const std::vector<double > verticesValues)
		{ m_spInterfaceHandlerLocal->write_solution(verticesValues); }

		// called during init() of diffusionInterface:
		void set_numDoFs(const size_t numDoFs)
		{ m_numDoFs = numDoFs;}
		void set_numNewDoFs(const size_t numNewDoFs)
		{ m_numNewDoFs = numNewDoFs;}

		void set_bScaleDoFs(bool bScaleDoF) { m_scaleDoFs = bScaleDoF; }

	private:
		SmartPtr<InterfaceHandlerLocalDiffusion<dim> > m_spInterfaceHandlerLocal;
	// number of DoFs in global matrix
		size_t m_numDoFs;
		size_t m_numNewDoFs;
		bool m_resized;
		bool m_resized_defect;
		bool m_scaleDoFs;


};

}// namespace ConvectionDiffusionPlugin
} // end ug namespace

#include "loc_to_glob_mapper_diffusion_impl.h"


#endif /* DIFFUSION_INTERFACE_MAPPER_H_ */
