
/*
 * moving_interface_diffusion_impl.h
 *
 *  Created on: 24.08.2017
 *      Author: suze
 */

#ifndef DIFFUSION_INTERFACE_IMPL_H_
#define DIFFUSION_INTERFACE_IMPL_H_


namespace ug {
namespace MovingInterfaceDiffusion {
        
///////////////////////////////////////////////////////////
// Implementation of the methods class
// 	 		'MovingInterfaceDiffusion'
///////////////////////////////////////////////////////////


template<typename TDomain, typename TAlgebra>
MovingInterfaceDiffusion<TDomain, TAlgebra>::MovingInterfaceDiffusion(
                            SmartPtr<IAssemble<TAlgebra> > ass,
                            SmartPtr<ConvectionDiffusionPlugin::ConvectionDiffusionFV1<TDomain> > spMaster,
                            SmartPtr<MovingInterfaceBase::DiffusionInterfaceProvider<dim> > interfaceProvider,
                            SmartPtr<MovingInterfaceBase::CutElementHandlerImmersed<dim> > cutElementHandler) :
    m_spInterfaceHandlerLocal(new InterfaceHandlerLocalDiffusion<dim>(interfaceProvider, cutElementHandler)),
    m_spInterfaceProvider(interfaceProvider),
    m_spCutElementHandler(cutElementHandler),
    m_spInterfaceMapper(new DiffusionInterfaceMapper<TDomain, TAlgebra> (m_spInterfaceHandlerLocal))
{
    if (interfaceProvider->num_particles() == 0)
        UG_THROW("MovingParticle::Constructor(): no particles initializen in 'globalHandler\n");
    
    // initialize singleton and set local handler
    typedef DimFV1CutGeometry<dim, dim, InterfaceHandlerLocalDiffusion<dim> > TFVGeom;
    TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1),1);
    geo.set_interface_handler(m_spInterfaceHandlerLocal);
    
    // initialize mapper within domainDisc:
    SmartPtr<AssemblingTuner<TAlgebra> > assAdapt = ass->ass_tuner();
	assAdapt->set_mapping(m_spInterfaceMapper.get());

    assAdapt->enable_modify_solution(true);
    // => assTuner->modify_LocSol() = mapper->modify_LocSol()
    // see: ass_tuner.h: 114
    
    
}

template<typename TDomain, typename TAlgebra>
double MovingInterfaceDiffusion<TDomain, TAlgebra>::
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
	{
		value = -2*kappa_2*sqDist*sqDist;
		UG_LOG("value = " << value << "\n");
	}
	return value;
}

template<typename TDomain, typename TAlgebra>
void MovingInterfaceDiffusion<TDomain, TAlgebra>::
set_analytic_solution(vector_type& u, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, SmartPtr<MultiGrid> mg, const int topLevel)
{
	ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));

	typedef MathVector<dim> position_type;
 	typedef Attachment<position_type> position_attachment_type;
 	typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

	//SmartPtr<MultiGrid> m_spMG = mg.operator->();
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
		//	get vertex
			Vertex* vrt = elem->vertex(i);
			double value = compute_solution_value(vCornerCoords[i]);
			DoFRef(u, ind[i]) = value;
		}
	}

}

template<typename TDomain, typename TAlgebra>
void MovingInterfaceDiffusion<TDomain, TAlgebra>::
adjust_for_error(vector_type& u, vector_type& uCopy, SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, SmartPtr<MultiGrid> mg, const int topLevel)
{


	size_t numDoFs = u.size();
	UG_LOG("domain disc: numDoFs = " << numDoFs << "\n");
	size_t numDoFsCopy = uCopy.size();
	UG_LOG("domain disc: numDoFsCopy = " << numDoFsCopy << "\n");

	u.resize(numDoFsCopy);

	numDoFs = u.size();
	UG_LOG("**domain disc: numDoFs = " << numDoFs << "\n");



	ConstSmartPtr<DoFDistribution> dd = spApproxSpace->dof_distribution(GridLevel(topLevel, GridLevel::LEVEL));

	typedef MathVector<dim> position_type;
 	typedef Attachment<position_type> position_attachment_type;
 	typedef Grid::VertexAttachmentAccessor<position_attachment_type> position_accessor_type;

	//SmartPtr<MultiGrid> m_spMG = mg.operator->();
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

		ElementModus elemModus = m_spInterfaceHandlerLocal->compute_element_modus(elem, 1);
		bool do_adjust = false;

		switch(elemModus)
		{
	 		case INSIDE_DOM:	   break;
	 		case OUTSIDE_DOM:	   break;

			case CUT_BY_INTERFACE: do_adjust = true; break;
			default:
				throw(UGError("Error in InterfaceHandlerLocalDiffusion::update(): switch(m_elemModus)!"));
		}


		std::vector<MathVector<dim> > vCornerCoords;
		CollectCornerCoordinates(vCornerCoords, *elem, m_aaPos);

		std::vector<DoFIndex> ind;
		dd->dof_indices(elem, 0, ind);

	//	loop vertices
		for(size_t i = 0; i < elem->num_vertices(); ++i)
		{
		//	get vertex
			Vertex* vrt = elem->vertex(i);

			if (do_adjust) DoFRef(u, ind[i]) = 0.0;
		}

	}

}

template<typename TDomain, typename TAlgebra>
const size_t MovingInterfaceDiffusion<TDomain, TAlgebra>::
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

 		ElementModus elemModus = m_spCutElementHandler->compute_element_modus(elem, 1);

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
const size_t MovingInterfaceDiffusion<TDomain, TAlgebra>::
initialize_interface(vector_type& u, ConstSmartPtr<DoFDistribution> dd)
{
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

 		ElementModus elemModus = m_spCutElementHandler->compute_element_modus(elem, 1);

		if ( elemModus == CUT_BY_INTERFACE )
			num_cutElements += 1;

	}

	return num_cutElements;

}

template<typename TDomain, typename TAlgebra>
number MovingInterfaceDiffusion<TDomain, TAlgebra>::
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

	UG_LOG("nach com.allreduce: mean = " << std::sqrt(mean) << "\n");
	return std::sqrt(mean);
}
    
template<typename TDomain, typename TAlgebra>
void MovingInterfaceDiffusion<TDomain, TAlgebra>::
initialize_threshold(TDomain& domain, const int baseLevel, const int topLevel)
{
	UG_LOG("----------------- START initialize_threshold() ---------------- \n");

	if (baseLevel < 0)
		UG_THROW(
				"initialize_threshold(): no cast of baselevel from 'int' tp 'size_t' possible! \n");
	if (topLevel < 0)
		UG_THROW(
				"initialize_threshold(): no cast of toplevel from 'int' tp 'size_t' possible! \n");

// compute level-dependent value for threshold:
	for (size_t lev = baseLevel; lev <= (size_t) topLevel; ++lev) {
		const number maxLength = MaxElementDiameter(domain, lev);
		const number meanLength = MeanElementDiameter(domain, lev);
		UG_LOG("maxLength = " << maxLength << "\n");
		UG_LOG("meanLength = " << meanLength << "\n");
		UG_LOG("threshold_max = " << maxLength*maxLength << "\n");
		UG_LOG("threshold_mean = " << meanLength*meanLength << "\n");

//		set_threshold(lev, meanLength * meanLength);
        set_threshold(lev, 0.25 * meanLength*meanLength);
	}

	UG_LOG("----------------- END initialize_threshold() ---------------- \n");

}

template<typename TDomain, typename TAlgebra>
void MovingInterfaceDiffusion<TDomain, TAlgebra>::
update_interface( const int topLevel, number deltaT)
{
	if ( deltaT ==  0.0 )
		UG_THROW("InterfaceProvider:update_prtCoords: deltaT = " << deltaT << " => no update necessary!\n");

// get level index
	const int levIndex = get_Index(GridLevel(topLevel, GridLevel::LEVEL));
	UG_LOG("update_prtCoords() for levIndex = " << levIndex << "\n");
	UG_LOG("update_prtCoords() for deltaT = " << deltaT << "\n");

// update center
	for (size_t p = 0; p < m_spInterfaceProvider->num_particles(); ++p)
	{
#ifdef UG_PARALLEL
        std::vector<grid_base_object*> ElemList = m_spCutElementHandler->m_vvvElemListCut[levIndex][p];
//		std::vector<grid_base_object*> ElemList = m_vvvElemListCut[levIndex][p];
 		UG_LOG("1_ update_prtCoords() ElemList.size(): " << ElemList.size() << "\n");
		if ( ElemList.size() == 0 ) {
 			UG_LOG("2_ update_prtCoords() ElemList.size(): " << ElemList.size() << " => skip assembling! \n");
			continue;
		}
#endif

	// get data:
  		MathVector<dim> centerNew = m_spInterfaceProvider->get_center(p);
  		number soution = m_spInterfaceProvider->get_solution(p, 0);
		UG_LOG(" solution = " << soution << "\n");
 		UG_LOG("deltaT = " << deltaT << "\n");

// ToDo: Hier das interface irgendwie updaten??


	} // end particle loop

}

} // end namespace MovingParticle
} // end namespace ug

#endif /* DIFFUSION_INTERFACE_IMPL_H_ */
