/*
 * diffusion_interface.h
 *
 *  Created on: 24.08.2017
 *      Author: suze
 */

#ifndef DIFFUSION_INTERFACE_DIFFUSION_H_
#define DIFFUSION_INTERFACE_DIFFUSION_H_


#ifdef UG_PARALLEL
 	#include "lib_grid/parallelization/load_balancer_util.h"
#endif

#include "../../../../ConvectionDiffusion/fv1/convection_diffusion_fv1.h"
#include "../../../../ConvectionDiffusion/convection_diffusion_base.h"
#include "interface_handler_diffusion.h"
#include "loc_to_glob_mapper_diffusion.h"
#include "../immersed_interface_base/immersed_interface_base.h"

namespace ug{
namespace MovingInterfaceDiffusion{



template <	typename TDomain, typename TAlgebra>
class MovingInterfaceDiffusion
		: public MovingInterfaceBase::IMovingInterface<TDomain, TAlgebra>
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

 		MovingInterfaceDiffusion(
 					   SmartPtr<IAssemble<TAlgebra> > ass,
 					   SmartPtr<ConvectionDiffusionPlugin::ConvectionDiffusionFV1<TDomain> > spMaster,
 					   SmartPtr<MovingInterfaceBase::DiffusionInterfaceProvider<dim> > interfaceProvider,
 					   SmartPtr<MovingInterfaceBase::CutElementHandlerImmersed<dim> > cutElementHandler);

		void set_StdFV_assembling(bool bValue) { m_spInterfaceHandlerLocal->set_StdFV_assembling(bValue);}
     	bool StdFV_assembling() { return m_spInterfaceHandlerLocal->StdFV_assembling(); }


 	// destructor
		~MovingInterfaceDiffusion(){};

	/// called via .lua:
		void initialize_threshold(TDomain& domain, const int baseLevel, const int topLevel);
	    void set_threshold(size_t level, const number threshold)
	    { m_spCutElementHandler->set_threshold(level, threshold); }

	//////////////////////////////////////////////////////////////////////////////////
	/// Info - 'initialize_interface()':
	///
	/// computes vertices on intersection of cut element edges and interface:
	/// for 2d: instead of computing intersections: count number of cut elements!
	/// 	-> #cutElements == #m_vertices
	/// 	-> called during init()
	//////////////////////////////////////////////////////////////////////////////////

	    const size_t initialize_interface(vector_type& u, ConstSmartPtr<DoFDistribution> dd);
	    const size_t initialize_interface_Nitsche(vector_type& u, ConstSmartPtr<DoFDistribution> dd);

		number MeanElementDiameter(TDomain& domain, int level);

		//////////////////////////////////////////////////////////////////////////////////
		// ---> .lua:  update(u, deltaT, u:grid_level()): update global indices (transInd...)
		// 				=> A. copy_solution(topLev)
		// 				   B. update(baseLev-topLev)
		// 				   C. update_solution(topLev)
 		//////////////////////////////////////////////////////////////////////////////////
		// write solution to nodes outside fluid with particle velocities
		//		--> call method vie .lua BEFORE 'solTimeSeries:push_discard_oldest(oldestSol, time)':
		//			=> in case that outside nodes are inside AFTER update_prtCoords: NO solution defined here!


		/// call of the method via lua to set the real velocity values within the particle domain
		void adjust_global_solution(vector_type& u, const int topLevel);
		void fill_particle_solution(vector_type& u, const int topLevel, const number time);
		void update(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel, const int topLevel, const number time)
		{
			int topLev = spApproxSpace->num_levels()-1;
			if ( topLev != topLevel )
				UG_THROW("update: parameter 'topLevel' = " << topLevel << " != "
								 << topLev << "current top leven! \n");

		// fill particle nodes with their real solution
 			fill_particle_solution(u, topLevel, time);

		// update data: bool_marker
			ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
			m_spCutElementHandler->template init<TDomain>(dd, baseLevel, topLevel);

  		}

		void set_interface_values(vector_type& u, const int numDoFs, const int num_newDoFs)
		{
			double value = 20.0*0.4*0.4*0.4*0.4;
			DoFIndex index;

			for (size_t i = 0; i < num_newDoFs; ++i)
			{
				index = DoFIndex(numDoFs + i,0);
		//		DoFRef(u, index) = value;
			}
		}
		void set_analytic_solution(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, SmartPtr<MultiGrid> mg, const int topLevel);
		void adjust_for_error(vector_type& u, vector_type& uCopy, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, SmartPtr<MultiGrid> mg, const int topLevel);

		double compute_solution_value(const MathVector<dim>& vrtPos);

		number get_integral()
		{ return m_spInterfaceHandlerLocal->get_integral(); }

		void set_Nitsche(bool bNitsche)
		{ return m_spInterfaceHandlerLocal->set_Nitsche(bNitsche); }

 		void init(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel, const int topLevel, bool bScaleDoFs)
		{
   			m_spApproxSpace = spApproxSpace;

			ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));

			m_spInterfaceMapper->set_numDoFs(u.size());

			size_t numDoFs = u.size();
			UG_LOG("domain disc: numDoFs = " << numDoFs << "\n");
			size_t num_newDoFs = initialize_interface(u, dd);
			UG_LOG("________________ num_newDoFs = " << initialize_interface(u, dd) << "\n");
			m_spInterfaceMapper->set_numNewDoFs(num_newDoFs);


			if ( m_spInterfaceHandlerLocal->get_Nitsche() )
			{
				const size_t buffer = initialize_interface_Nitsche(u, dd);
				UG_LOG(" buffer = " << buffer << "\n");

				const size_t num_NitscheDoFs = m_spInterfaceHandlerLocal->get_num_NitscheDoFs();
				UG_LOG(" num_NitscheDoFs = " << num_NitscheDoFs << "\n");

				num_newDoFs = num_NitscheDoFs;
				m_spInterfaceMapper->set_numNewDoFs(num_newDoFs);

				for ( size_t i = 0; i < num_NitscheDoFs; ++i )
				{
					size_t global_index = m_spInterfaceHandlerLocal->get_global_index_Nitsche(i);
					UG_LOG(" global_index = " << global_index << "\n");
				}
			}

		// values for new DoFs are set to 0.0 by the 'resize()'-method (see vector.h):
			if ( bScaleDoFs )
				u.resize(numDoFs + 2*num_newDoFs);
			else
				u.resize(numDoFs + num_newDoFs);

			UG_LOG("AGAIN: in init(): numDoFs = " <<  u.size() << "\n");

			m_spInterfaceMapper->set_bScaleDoFs(bScaleDoFs);
			m_spInterfaceHandlerLocal->set_bScaleDoFs(bScaleDoFs);

			m_spInterfaceHandlerLocal->init_integral();

 		// lieber in jedem Schritt über 'modify_GlobalSol()' (mapper!) setzten!
			//set_interface_values(u, numDoFs, num_newDoFs);

			// not necessary anymore: only local evaluations within diffusion problem!
			//m_spCutElementHandler->template init_marker<TDomain>(dd, baseLevel, topLevel);
		}

	   /// checks if grid data is updated and returns 'levIndex'-pair for 'gridLevel' in 'm_Map'
	    int get_Index(const GridLevel& gridLevel)
	    {
	    	ConstSmartPtr<DoFDistribution> dd = m_spApproxSpace->dof_distribution(gridLevel);

 	 		const int levIndex = m_spCutElementHandler->get_Index(gridLevel, dd);

	    	return levIndex;
	    }
    
	    //ToDo: method needed?
	    void update_interface( const int topLevel, number deltaT);

		bool is_time_dependent() { return m_spInterfaceHandlerLocal->is_time_dependent();}

		/// helper functions for compute_error_on_circle()
		void interpolate_point(ConstSmartPtr<DoFDistribution> dd, const vector_type& u,
							   const MathVector<dim>& evalPos, MathVector<dim+1>& interpolation);
		void compute_error_on_circle(const vector_type& u, const int topLevel, number radius);

		void print_deltaP(const vector_type& u, const int topLevel);
		void print_pressure(const vector_type& u, const int topLevel);
		void print_pressure_nodal(const vector_type& u, const int topLevel);

	/// writing data to file; called via .lua
		void print_velocity(const vector_type& u, const int topLevel, number time, const char* filename);
		void print_velocity_many_particles(const vector_type& u, const int topLevel, number time, const char* filename);

        size_t get_numDoFs(const vector_type& u) {return u.size(); }

	private:
	///	current ApproxSpace
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

		SmartPtr<MovingInterfaceBase::DiffusionInterfaceProvider<dim> > m_spInterfaceProvider;
		SmartPtr<MovingInterfaceBase::CutElementHandlerImmersed<dim> > m_spCutElementHandler;
		SmartPtr<InterfaceHandlerLocalDiffusion<dim> > m_spInterfaceHandlerLocal;

		SmartPtr<DiffusionInterfaceMapper<TDomain, TAlgebra> > m_spInterfaceMapper;

};

} // end namespace MovingInterfaceDiffusion
} // end namespace ug


#include "diffusion_interface_impl.h"

#endif /* DIFFUSION_INTERFACE_DIFFUSION_H_ */




