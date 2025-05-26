/*
 * Copyright (c) 2012-2014:  G-CSC, Goethe University Frankfurt
 * Authors: Daniel Gonzalez
 * Based on the modules by Dmitry Logashenko, Andreas Vogel
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


#include "outflow_fv1.h"

#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusionOutflowFV1<TDomain>::
ConvectionDiffusionOutflowFV1(SmartPtr< ConvectionDiffusionBase<TDomain> > spMaster)
: ConvectionDiffusionOutflowBase<TDomain>(spMaster)
{
    m_imVelocity.set_comp_lin_defect(false);

//	register imports
	this->register_import(m_imVelocity);


//	initialize the imports from the master discretization
	set_velocity(spMaster->velocity ());

	//	update assemble functions
	register_all_funcs(false);
};


template<typename TDomain>
void ConvectionDiffusionOutflowFV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	if(bNonRegularGrid)
		UG_THROW("ConvectionDiffusion: only regular grid implemented.");

//	check number
	if(vLfeID.size() != 1)
		UG_THROW("ConvectionDiffusion: Need exactly "<<1<<" functions");

    if(vLfeID[0].type() != LFEID::LAGRANGE || vLfeID[0].order() != 1)
        UG_THROW("ConvectionDiffusion: 'fv1' expects Lagrange P1 trial space "
					"for concentration.");
    

	//	update assemble functions
	register_all_funcs(false);
}

////////////////////////////////////////////////////////////////////////////////
//	assembling functions
////////////////////////////////////////////////////////////////////////////////

/**
 * Prepares the element loop for a given element type: computes the FV-geo, ...
 * Note that there are separate loops for every type of the grid elements.
 */
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionOutflowFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
//	register subsetIndex at Geometry
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	Only first order implementation
	if(!(TFVGeom::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
						" Geometry set.");

//	check if velocity  has been set
	if(!m_imVelocity.data_given())
		UG_THROW("OutflowFV1::prep_elem_loop:"
						" velocity  has not been set, but is required.\n");

//	extract indices of boundary
	this->extract_scheduled_data();

//	request the subset indices as boundary subset. This will force the
//	creation of boundary subsets when calling geo.update
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.add_boundary_subset(*subsetIter);
}

/**
 * Finalizes the element loop for a given element type.
 */
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionOutflowFV1<TDomain>::
fsh_elem_loop()
{
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	remove the bnd subsets
	typename std::vector<int>::const_iterator subsetIter;
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
		geo.remove_boundary_subset(*subsetIter);
}


/**
 * General initializations of a given grid element for the assembling.
 */
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionOutflowFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();
	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("OutflowFV1::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	find and set the local and the global positions of the IPs for imports
	typedef typename TFVGeom::BF BF;
	typename std::vector<int>::const_iterator subsetIter;
	
	m_vLocIP.clear(); m_vGloIP.clear();
	for(subsetIter = m_vBndSubSetIndex.begin();
			subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
		const int bndSubset = *subsetIter;
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);
		for(size_t i = 0; i < vBF.size(); ++i)
		{
			m_vLocIP.push_back(vBF[i].local_ip());
			m_vGloIP.push_back(vBF[i].global_ip());
		}
	}
	// REMARK: The loop above determines the ordering of the integration points:
	// The "outer ordering" corresponds to the ordering of the subsets in
	// m_vBndSubSetIndex, and "inside" of this ordering, the ip's are ordered
	// according to the order of the boundary faces in the FV geometry structure.

    m_imVelocity.set_local_ips(&m_vLocIP[0], m_vLocIP.size());
    m_imVelocity.set_global_ips(&m_vGloIP[0], m_vGloIP.size());
	
}


/// Assembling of the convective flux  in the Jacobian of the trasport eq.
template<typename TDomain>
template<typename BF>
void ConvectionDiffusionOutflowFV1<TDomain>::
convective_flux_Jac
(
	const size_t ip, // index of the integration point (for the velocity)
	const BF& bf, // boundary face to assemble
	LocalMatrix& J, // local Jacobian to update
	const LocalVector& u, // local solution
    const number& stdValue // concentration at ip
)
{
// The convection velocity according to the current approximation:
	number volumetric_flux = VecDot (m_imVelocity[ip], bf.normal ()) * 1;
	
// We assume that there should be no inflow through the outflow boundary:
	if (volumetric_flux < 0)
        volumetric_flux = 0;
	
//	Add flux to local Jacobian
	for(size_t sh = 0; sh < bf.num_sh(); ++sh)
	{
		number t = volumetric_flux * bf.shape(sh);
        J(_C_, bf.node_id(), _C_, sh) += t;
	}
}

/// Assembling of the convective flux  in the defect of the trasport eq.
template<typename TDomain>
template<typename BF>
void ConvectionDiffusionOutflowFV1<TDomain>::
convective_flux_defect
(
	const size_t ip, // index of the integration point (for the velocity)
	const BF& bf, // boundary face to assemble
	LocalVector& d, // local defect to update
	const LocalVector& u, // local solution
	const number& stdValue // concetration at ip
)
{
// The volumetric flux at boundary:
	number volumetric_flux = VecDot (m_imVelocity[ip], bf.normal ());
	
// We assume that there should be no inflow through the outflow boundary:
	if (volumetric_flux < 0)
        volumetric_flux = 0;
	
// Add the flux to the defect:

    d(_C_, bf.node_id()) += volumetric_flux * stdValue;
}



template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionOutflowFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	typedef typename TFVGeom::BF BF;

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ip = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

	// 	loop the boundary faces
		typename std::vector<BF>::const_iterator bf;
		for(bf = vBF.begin(); bf != vBF.end(); ++bf, ++ip)
		{
            
            number stdValue=0;
            convective_flux_Jac<BF> (ip, *bf, J, u, stdValue);

			

		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionOutflowFV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	Only first order implemented
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	typedef typename TFVGeom::BF BF;

// 	loop registered boundary segments
	typename std::vector<int>::const_iterator subsetIter;
	size_t ip = 0;
	for(subsetIter = m_vBndSubSetIndex.begin();
		subsetIter != m_vBndSubSetIndex.end(); ++subsetIter)
	{
	//	get subset index corresponding to boundary
		const int bndSubset = *subsetIter;
		
	//	get the list of the ip's:
		if(geo.num_bf(bndSubset) == 0) continue;
		const std::vector<BF>& vBF = geo.bf(bndSubset);

	// 	loop the boundary faces
		typename std::vector<BF>::const_iterator bf;
		for(bf = vBF.begin(); bf != vBF.end(); ++bf, ++ip)
		{

        // A. Compute concentration at ip
            number stdValue=0;
            for(size_t sh = 0; sh < bf->num_sh(); ++sh)
            {
                stdValue += u(_C_, sh) * bf->shape(sh);
            }

            convective_flux_defect<BF> (ip, *bf, d, u, stdValue);
            
		
		}
	}
}





////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void ConvectionDiffusionOutflowFV1<Domain1d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		UG_THROW("OutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_2
template<>
void ConvectionDiffusionOutflowFV1<Domain2d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
	}
	else
	{
		UG_THROW("OutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_3
template<>
void ConvectionDiffusionOutflowFV1<Domain3d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
	}
	else
	{
		UG_THROW("OutflowFV1: Hanging Nodes not implemented.")
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void
ConvectionDiffusionOutflowFV1<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(	id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 	id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( 	id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(	id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(	id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(	id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(	id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(	id, &T::template add_rhs_elem<TElem, TFVGeom>);
}


////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ConvectionDiffusionOutflowFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConvectionDiffusionOutflowFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionOutflowFV1<Domain3d>;
#endif

} // namespace ConvectionDiffusionPlugin
} // namespace ug
