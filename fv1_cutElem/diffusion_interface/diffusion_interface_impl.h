
/*
 * diffusion_interface_impl.h
 *
 *  Created on: 24.08.2017
 *      Author: suze
 */

#ifndef IMMERSED_INTERFACE_DIFFUSION_IMPL_H_
#define IMMERSED_INTERFACE_DIFFUSION_IMPL_H_

#include "lib_disc/spatial_disc/disc_util/geom_provider.h" 

namespace ug {
namespace ConvectionDiffusionPlugin {
        
///////////////////////////////////////////////////////////
// Implementation of the methods class
// 	 		'ImmersedInterfaceDiffusion'
///////////////////////////////////////////////////////////


template<typename TDomain, typename TAlgebra>
ImmersedInterfaceDiffusion<TDomain, TAlgebra>::ImmersedInterfaceDiffusion(
                            SmartPtr<IAssemble<TAlgebra> > ass,
                            SmartPtr<ConvectionDiffusionPlugin::ConvectionDiffusionFV1_cutElem<TDomain> > spMaster,
                            SmartPtr<DiffusionInterfaceProvider<dim> > interfaceProvider,
                            SmartPtr<CutElementHandler_TwoSided<dim> > cutElementHandler) :
    m_spInterfaceProvider(interfaceProvider),
    m_spCutElementHandler(cutElementHandler),
    m_spInterfaceHandlerLocal(new InterfaceHandlerLocalDiffusion<dim>(interfaceProvider, cutElementHandler)),
    m_spInterfaceMapper(new DiffusionInterfaceMapper<TDomain, TAlgebra> (m_spInterfaceHandlerLocal))
{
    if (interfaceProvider->num_particles() == 0)
        UG_THROW("ImmersedInterfaceDiffusion::Constructor(): no particles initializen in 'globalHandler\n");
    
// initialize singleton and set local handler
    typedef DimFV1CutGeometry<dim, dim, InterfaceHandlerLocalDiffusion<dim> > TFVGeom;
    TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1),1);
    geo.set_interface_handler(m_spInterfaceHandlerLocal);
    
// initialize mapper within domainDisc:
    SmartPtr<AssemblingTuner<TAlgebra> > assAdapt = ass->ass_tuner();
	assAdapt->set_mapping(m_spInterfaceMapper.get());

// needs to be enabled, in order to call 'spAssTuner->modify_LocalData()' during element disc assembling
    assAdapt->enable_modify_solution(true);
    
    
}

    
template<typename TDomain, typename TAlgebra>
void ImmersedInterfaceDiffusion<TDomain, TAlgebra>::
init(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, const int baseLevel,
              const int topLevel, bool bScaleDoFs)
{
 // get data
    m_spApproxSpace = spApproxSpace;
    ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
    
  // get the level Index ONLY for the toplevel in order to resize the data for the toplevel
  //   --> not implemented for multigrid already!
    const int levIndex = m_spCutElementHandler->get_Index(topLevel, dd);

 // initialize the orientation of the interface: = 1, i.e. inside the circle line is outside the domain
    m_spCutElementHandler->set_orientation(1);
    
// call of 'update_interface_data()' (via init() of CutElementHandler):
// sets the marker (INSIDE, OUTSIDE, CUT_BY_INTERFACE) to each element
//  and fills the lists of cut elements
    m_spCutElementHandler->template init<TDomain>(dd, baseLevel, topLevel);

// store the number of regular DoFs (u.size()) and additional DoFs on the immersed interface:
    size_t numDoFs = u.size();
    size_t num_newDoFs = m_spCutElementHandler->get_numCutElements(levIndex);
    m_spInterfaceMapper->set_numDoFs(numDoFs);
    m_spInterfaceMapper->set_numNewDoFs(num_newDoFs);

    
    if ( m_spInterfaceHandlerLocal->get_Nitsche() )
    {
        initialize_interface_Nitsche(u, dd);
        const size_t num_NitscheDoFs = m_spInterfaceHandlerLocal->get_num_NitscheDoFs();
            
        num_newDoFs = num_NitscheDoFs;
        m_spInterfaceMapper->set_numNewDoFs(num_newDoFs);
            
     }
 
// initialize the scaling of new DoFs: *2 for jumping values (strong discontinuity) at interface
//  => 1 DoFs is defined for the solution on each side of the interface
// values for new DoFs are set to 0.0 by the 'resize()'-method (see vector.h):
    if ( bScaleDoFs )   u.resize(numDoFs + 2*num_newDoFs);
    else                u.resize(numDoFs + num_newDoFs);
    
    m_spInterfaceMapper->set_bScaleDoFs(bScaleDoFs);
    m_spInterfaceHandlerLocal->set_bScaleDoFs(bScaleDoFs);
    
// initialize the value for the l2 error computaton
    m_spInterfaceHandlerLocal->L2Error_init();
          
}

template<typename TDomain, typename TAlgebra>
void ImmersedInterfaceDiffusion<TDomain, TAlgebra>::
update(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                const int baseLevel, const int topLevel, const number time)
{
// get data
    ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));
    int topLev = spApproxSpace->num_levels()-1;
    if ( topLev != topLevel )
        UG_THROW("update: parameter 'topLevel' = " << topLevel << " != "
                     << topLev << "current top leven! \n");
    
// update data: update BoolMarker and cut element lists
    m_spCutElementHandler->template init<TDomain>(dd, baseLevel, topLevel);
        
}
    
template<typename TDomain, typename TAlgebra>
double ImmersedInterfaceDiffusion<TDomain, TAlgebra>::
compute_solution_value(const MathVector<dim>& vrtPos)
{
	double kappa_1 = 1.0;
	double kappa_2 = 10.0;
	double sqR = 0.4*0.4;
	double dist_x = vrtPos[0] - 0.1;
	double dist_y = vrtPos[1] - 0.2;
	double sqDist = dist_x*dist_x+dist_y*dist_y;

	double value = -4*kappa_1*kappa_2*kappa_2*sqR*sqDist + 2*sqR*sqR*kappa_2*(2*kappa_1*kappa_2 - 1);

	double dist = sqrt(sqDist);
	if ( dist > 0.4 )
		value = -2*kappa_2*sqDist*sqDist;
	
	return value;
}

template<typename TDomain, typename TAlgebra>
void ImmersedInterfaceDiffusion<TDomain, TAlgebra>::
set_analytic_solution(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, SmartPtr<MultiGrid> mg, const int topLevel)
{
	ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));

	typedef MathVector<dim> position_type;
 	typedef Attachment<position_type> position_attachment_type;
 	typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

	position_attachment_type m_aPos = GetDefaultPositionAttachment<position_attachment_type>();
	position_accessor_type m_aaPos;

	if(!mg->has_attachment<Vertex>(m_aPos))
		mg->attach_to<Vertex>(m_aPos);
	m_aaPos.access(*mg, m_aPos);

	typedef typename domain_traits<dim>::grid_base_object grid_base_object;

	typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
	iter = dd->template begin<grid_base_object>();
	iterEnd = dd->template end<grid_base_object>();

	//	loop elements in order to compute the volume and set rhs:
	for( ; iter != iterEnd; iter++)
	{
 	//	get element
		grid_base_object* elem = *iter;
		std::vector<MathVector<dim> > vCornerCoords;
		CollectCornerCoordinates(vCornerCoords, *elem, m_aaPos);

		std::vector<DoFIndex> ind;
		dd->dof_indices(elem, 0, ind);

	//	loop vertices
		for(size_t i = 0; i < elem->num_vertices(); ++i)
		{
			double value = compute_solution_value(vCornerCoords[i]);
			DoFRef(u, ind[i]) = value;
		}
	}

}

template<typename TDomain, typename TAlgebra>
void ImmersedInterfaceDiffusion<TDomain, TAlgebra>::
adjust_for_error(vector_type& u, vector_type& uCopy, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
                 SmartPtr<MultiGrid> mg, const int topLevel)
{
//  (1) remove additional nodes ON the interface, not lying on the original grid
//      --> resize vector 'u' to original size, i.e. decrease size
 	size_t numDoFsCopy = uCopy.size();
 
	u.resize(numDoFsCopy);


// get data
	typedef MathVector<dim> position_type;
 	typedef Attachment<position_type> position_attachment_type;
 	typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

 	position_attachment_type m_aPos = GetDefaultPositionAttachment<position_attachment_type>();
	position_accessor_type m_aaPos;

	if(!mg->has_attachment<Vertex>(m_aPos)) mg->attach_to<Vertex>(m_aPos);
	m_aaPos.access(*mg, m_aPos);

    ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));

    typedef typename domain_traits<dim>::grid_base_object grid_base_object;
	typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
	iter = dd->template begin<grid_base_object>();
	iterEnd = dd->template end<grid_base_object>();

//  (2) set solution to zero in near-interface nodes
//      --> loop elements in order to set solution to zero on cut element nodes of vector 'u':
	for( ; iter != iterEnd; iter++)
	{
 	//	get element
		grid_base_object* elem = *iter;

		ElementModus elemModus = m_spInterfaceHandlerLocal->compute_element_modus(elem);

    // choose only cut elements for adjustment
        bool do_adjust = false;
		switch(elemModus)
		{
	 		case INSIDE_DOM:        break;
	 		case OUTSIDE_DOM:       break;

			case CUT_BY_INTERFACE:  do_adjust = true;
                                    break;
			default:
				throw(UGError("Error in InterfaceHandlerLocalDiffusion::update(): switch(m_elemModus)!"));
		}
        
        if ( !do_adjust ) continue;

    //	get local indices
		std::vector<DoFIndex> ind;
		dd->dof_indices(elem, 0, ind);
        
    // loop vertices and set solution to zero
		for(size_t i = 0; i < elem->num_vertices(); ++i)
            DoFRef(u, ind[i]) = 0.0;

	}

}

template<typename TDomain, typename TAlgebra>
const size_t ImmersedInterfaceDiffusion<TDomain, TAlgebra>::
initialize_interface_Nitsche(vector_type& u, ConstSmartPtr<DoFDistribution> dd)
{
	m_spInterfaceHandlerLocal->m_MapInserted_Nitsche.clear();

	typedef typename domain_traits<dim>::grid_base_object grid_base_object;

	typename DoFDistribution::traits<grid_base_object>::const_iterator iter, iterEnd;
	iter = dd->template begin<grid_base_object>();
	iterEnd = dd->template end<grid_base_object>();

	size_t num_cutElements = 0;

	//	loop elements in order to compute the volume and set rhs:
	for( ; iter != iterEnd; iter++)
	{
 	//	get element
		grid_base_object* elem = *iter;

 		ElementModus elemModus = m_spCutElementHandler->compute_element_modus(elem);

		if ( elemModus == CUT_BY_INTERFACE )
		{
			num_cutElements += 1;

		// 	get global indices
			LocalIndices ind;
			dd->indices(elem, ind, false);

		// fill 'm_MapInserted_Nitsche' with global indices:
			for(size_t i = 0; i < 3; ++i)
				m_spInterfaceHandlerLocal->get_or_insert_indexPair_Nitsche(ind.index(0, i));

		}
	}

	return num_cutElements;

}


template<typename TDomain, typename TAlgebra>
number ImmersedInterfaceDiffusion<TDomain, TAlgebra>::
MeanElementDiameter(TDomain& domain, int level)
{
	typedef typename domain_traits<TDomain::dim>::grid_base_object TElem;
	typedef typename geometry_traits<TElem>::iterator ListIter;

	ListIter iter = domain.grid()->template begin<TElem>(level);
	ListIter iterEnd = domain.grid()->template end<TElem>(level);

	number mean = 0.0;
	size_t numIter = 0;
	for (; iter != iterEnd; ++iter) {
		mean += ElementDiameterSq(*domain.grid(), domain.position_accessor(),
				*iter);
		numIter++;
	}

	mean = mean / numIter;

#ifdef UG_PARALLEL
// share value between all procs
	pcl::ProcessCommunicator com;
// ToDO: PCL_RO_MIN oder MAX oder doch mean noch berechnen m.H.v. /numProcs???
//UG_THROW("in MeanElementDiameter(): was macht 'allredue' im Parallelen???\n");
	mean = com.allreduce(mean, PCL_RO_MIN);
#endif

	return std::sqrt(mean);
}
    
template<typename TDomain, typename TAlgebra>
void ImmersedInterfaceDiffusion<TDomain, TAlgebra>::
initialize_threshold(TDomain& domain, const int baseLevel, const int topLevel)
{
    UG_THROW("initialize_threshold(): not tested for diffusion case! \n");

	if (baseLevel < 0)
		UG_THROW(
				"initialize_threshold(): no cast of baselevel from 'int' tp 'size_t' possible! \n");
	if (topLevel < 0)
		UG_THROW(
				"initialize_threshold(): no cast of toplevel from 'int' tp 'size_t' possible! \n");

// compute level-dependent value for threshold:
	for (size_t lev = baseLevel; lev <= (size_t) topLevel; ++lev)
    {
		//const number maxLength = MaxElementDiameter(domain, lev);
		const number meanLength = MeanElementDiameter(domain, lev);

//		set_threshold(lev, meanLength * meanLength);
        set_threshold(lev, 0.25 * meanLength*meanLength);
	}


}


} // end namespace ConvectionDiffusionPlugin
} // end namespace ug

#endif /* IMMERSED_INTERFACE_DIFFUSION_IMPL_H_ */
