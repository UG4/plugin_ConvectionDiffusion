/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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


#include "convection_diffusion_fv1_cutElem.h"

#include "lib_disc/spatial_disc/disc_util/fv1Cut_geom.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"

namespace ug{
namespace ConvectionDiffusionPlugin{


DebugID DID_CONV_DIFF_FV1_CUTELEM("CONV_DIFF_FV1_CUTELEM");

////////////////////////////////////////////////////////////////////////////////
//	helper function
////////////////////////////////////////////////////////////////////////////////

/// helper function to prepare data for 'add_def_A_elem_local()' and 'add_jac_A_elem_local()'
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
get_local_data(TFVGeom& geo, const LocalVector& u, LocalVector& locUOut, MathMatrix<dim, dim> diffusionOut, LocalVector& jumpOut, LocalVector& jump_gradOut, LocalVector& sourceOut)
{
    const bool bElementIsCut = geo.get_element_modus();
    
    LocalIndices ind = u.get_indices();
    jumpOut.resize(ind); jump_gradOut.resize(ind); sourceOut.resize(ind);

    const int orientation = geo.get_orientation();

    std::vector<double> imSource;
    if ( m_imSource.data_given() ) {
        for ( size_t i = 0; i < 3; ++i )
            imSource.push_back(m_imSource[i]);
    }
    
// A. set data for non-cut elements
    if ( !bElementIsCut )
    {
     // set diffusion
        diffusionOut *= geo.get_diffusion(geo.get_boolian_for_diffusion());
        
     // set boundary condition values source, jump, jump_grad
        jumpOut = 0.0; jump_gradOut = 0.0;
        sourceOut = geo.set_source(imSource, ind, 3, false);
        LocalVector source = geo.set_source(imSource, ind, 3, false);
        
        int s = 0;
    }
// B. set data for cut element
    else
    {
    // set diffusion
        diffusionOut *= geo.get_diffusion();

    // set boundary condition values source, jump, jump_grad
        int indexSize = 3;
        if ( geo.get_roid() == ROID_QUADRILATERAL )
            indexSize = 4;
            
        geo.set_local_sol(locUOut, indexSize, u, orientation);
            
        jumpOut      = geo.set_jump_values(ind, indexSize);
        jump_gradOut = geo.set_jump_grad_values(ind, indexSize);
        sourceOut    = geo.set_source(imSource, ind, indexSize, true);
    }
    
}
        
////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusionFV1_cutElem<TDomain>::
ConvectionDiffusionFV1_cutElem(const char* functions, const char* subsets)
 : ConvectionDiffusionBase<TDomain>(functions,subsets),
   m_spConvShape(new ConvectionShapesNoUpwind<dim>),
   m_bNonRegularGrid(false)
{
	register_all_funcs(m_bNonRegularGrid);
}

template<typename TDomain>
void ConvectionDiffusionFV1_cutElem<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
//	check number
	if(vLfeID.size() != 1)
		UG_THROW("ConvectionDiffusion: Wrong number of functions given. "
				"Need exactly "<<1);

	if(vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("ConvectionDiffusion FV Scheme only implemented for 1st order.");

    m_LFEID = vLfeID[0];

//	remember
	m_bNonRegularGrid = bNonRegularGrid;

//	update assemble functions
	register_all_funcs(m_bNonRegularGrid);
}

template<typename TDomain>
bool ConvectionDiffusionFV1_cutElem<TDomain>::
use_hanging() const
{
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// Assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ConvectionDiffusionFV1_cutElem<TDomain>::
prep_assemble_loop()
{
	if (m_sss.valid())
		m_sss->clear_markers();
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{

//	check, that upwind has been set
	if(m_spConvShape.invalid())
		UG_THROW("ConvectionDiffusionFV1_cutElem::prep_elem_loop:"
						" Upwind has not been set.");

//	set local positions
//	if(!TFVGeom::usesHangingNodes)
    if(!TFVGeom::usesHangingNodes && TFVGeom::staticLocalData)
	{
		static const int refDim = TElem::dim;
//		TFVGeom& geo = GeomProvider<TFVGeom>::get();
        TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID, 1);
        const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		m_imDiffusion.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imVelocity.template 		set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imFlux.template 			set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imVectorSource.template 	set_local_ips<refDim>(vSCVFip,numSCVFip, false);
		m_imReactionRate.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReaction.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReactionRateExpl.template set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imReactionExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imSourceExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMassScale.template 		set_local_ips<refDim>(vSCVip,numSCVip, false);
		m_imMass.template 			set_local_ips<refDim>(vSCVip,numSCVip, false);

		//	init upwind for element type
		if(!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
			UG_THROW("ConvectionDiffusionFV1_cutElem::prep_elem_loop:"
						" Cannot init upwind for element type.");
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
    //static TFVGeom& geo = GeomProvider<TFVGeom>::get();
    TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID,1);
    //TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1),1);
    
// fix: set orientation initially globally!
    geo.set_orientation(1);

    
	try{
		UG_DLOG(DID_CONV_DIFF_FV1_CUTELEM, 2, ">>OCT_DISC_DEBUG: " << "convection_diffusion_fv1.cpp: " << "prep_elem(): update(): "<< roid << std::endl);
//        geo.update(elem, roid, vCornerCoords, &(this->subset_handler()));
        geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}UG_CATCH_THROW("ConvectionDiffusionFV1_cutElem::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	set local positions
//	if(TFVGeom::usesHangingNodes)
    if(TFVGeom::usesHangingNodes || !TFVGeom::staticLocalData)
	{
		const int refDim = TElem::dim;
		const MathVector<refDim>* vSCVFip = geo.scvf_local_ips();
		const size_t numSCVFip = geo.num_scvf_ips();
		const MathVector<refDim>* vSCVip = geo.scv_local_ips();
		const size_t numSCVip = geo.num_scv_ips();
		m_imDiffusion.template 		set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imVelocity.template 		set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imFlux.template 			set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imSource.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imVectorSource.template 	set_local_ips<refDim>(vSCVFip,numSCVFip);
		m_imReactionRate.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReaction.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReactionRateExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imReactionExpl.template 	set_local_ips<refDim>(vSCVip,numSCVip);
		m_imSourceExpl.template		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imMassScale.template 		set_local_ips<refDim>(vSCVip,numSCVip);
		m_imMass.template 			set_local_ips<refDim>(vSCVip,numSCVip);
/*
		if(m_spConvShape.valid())
			if(!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
				UG_THROW("ConvectionDiffusionFV1_cutElem::prep_elem_loop:"
								" Cannot init upwind for element type.");
 */
	}

	//	set global positions
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();
	m_imDiffusion.			set_global_ips(vSCVFip, numSCVFip);
	m_imVelocity.			set_global_ips(vSCVFip, numSCVFip);
	m_imFlux.				set_global_ips(vSCVFip, numSCVFip);
	m_imSource.				set_global_ips(vSCVip, numSCVip);
	m_imVectorSource.		set_global_ips(vSCVFip, numSCVFip);
	m_imReactionRate.		set_global_ips(vSCVip, numSCVip);
	m_imReactionRateExpl.	set_global_ips(vSCVip, numSCVip);
	m_imReactionExpl.		set_global_ips(vSCVip, numSCVip);
	m_imSourceExpl.			set_global_ips(vSCVip, numSCVip);
	m_imReaction.			set_global_ips(vSCVip, numSCVip);
	m_imMassScale.			set_global_ips(vSCVip, numSCVip);
	m_imMass.				set_global_ips(vSCVip, numSCVip);
}

template <class TVector>
static TVector CalculateCenter(GridObject* o, const TVector* coords)
{
	TVector v;
	VecSet(v, 0);

	size_t numCoords = 0;
	switch(o->base_object_id()){
		case VERTEX: numCoords = 1; break;
		case EDGE: numCoords = static_cast<Edge*>(o)->num_vertices(); break;
		case FACE: numCoords = static_cast<Face*>(o)->num_vertices(); break;
		case VOLUME: numCoords = static_cast<Volume*>(o)->num_vertices(); break;
		default: UG_THROW("Unknown element type."); break;
	}

	for(size_t i = 0; i < numCoords; ++i)
		VecAdd(v, v, coords[i]);

	if(numCoords > 0)
		VecScale(v, v, 1. / (number)numCoords);

	return v;
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
    
    bool debug = false;
    bool boundary = false;
        
    // get finite volume geometry
    //	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID,1);
    const bool bElementIsCut = geo.get_element_modus();
    
    // First call with orientation = 1:
    int orientation = 1;
    geo.set_orientation(orientation);
    
    geo.init_integral();
    
    // normal assembling if not cut by interface:
    if ( !bElementIsCut )
    {
        LocalVector dummyU;
        LocalIndices ind = u.get_indices();
        dummyU.resize(ind);
        dummyU = 0;
        
        this->template add_jac_A_elem_local<TElem,TFVGeom> (geo, u, J, dummyU, elem, vCornerCoords);
        return;
    }
        
    // get data:
    geo.resize_local_data(u);
    LocalMatrix& locJ_tri  = geo.get_jacobian_tri();
    LocalMatrix& locJ_quad = geo.get_jacobian_quad();
    
    LocalVector& locU_tri  = geo.get_solution_tri();
    LocalVector& locU_quad = geo.get_solution_quad();
    
    // reset data:
    locJ_tri = 0;
    locJ_quad = 0;
    
    LocalIndices ind = u.get_indices();
    
    // call elem disc twice:
    
    if ( debug ) geo.print_InterfaceIDdata();
        
    ReferenceObjectID roidCheck = geo.get_roid();
    if ( roidCheck == ROID_TRIANGLE )
    {
        geo.set_local_sol(locU_tri, 3, u, orientation);
        LocalVector jump_tri = geo.set_jump_values(ind, 3);
        
        this->template add_jac_A_elem_local<TElem,TFVGeom> (geo, u, locJ_tri, locU_tri, elem, vCornerCoords);
        
        if ( boundary )
        {
            geo.reset_jacobian_on_interface(locJ_tri, 3);
            this->template add_jac_A_elem_boundary<TElem,TFVGeom> (geo, locJ_tri, locU_tri, elem, vCornerCoords);
        }
        
        geo.set_jacobian_tri(locJ_tri);
        geo.set_DoF_tag_tri(false);
    }
    if ( roidCheck == ROID_QUADRILATERAL )
    {
        geo.set_local_sol(locU_quad, 4, u, orientation);
        LocalVector jump_quad = geo.set_jump_values(ind, 4);
        
        this->template add_jac_A_elem_local<TElem,TFVGeom> (geo, u, locJ_quad, locU_quad, elem, vCornerCoords);
        
        if ( boundary )
        {
            geo.reset_jacobian_on_interface(locJ_quad, 4);
            this->template add_jac_A_elem_boundary<TElem,TFVGeom> (geo, locJ_quad, locU_quad, elem, vCornerCoords);
        }
        
        geo.set_jacobian_quad(locJ_quad);
        geo.set_DoF_tag_quad(false);
    }
        
        
    // Second call with orientation = 1:
    bool shiftTag = geo.get_bScaleDoFs(); // shiftTag = true in case of double DoFs on interface!
    
    orientation *= -1;
    geo.set_orientation(orientation);
    try{
        geo.update(elem, vCornerCoords, &(this->subset_handler()));
    }UG_CATCH_THROW("ConvectionDiffusionFV1_cutElem::update:"
                    " Cannot update Finite Volume Geometry.");
    
    if ( debug ) geo.print_InterfaceIDdata();
    
    roidCheck = geo.get_roid();
    if ( roidCheck == ROID_TRIANGLE )
    {
        geo.set_local_sol(locU_tri, 3, u, orientation);
        LocalVector jump_tri = geo.set_jump_values(ind, 3);
        
        this->template add_jac_A_elem_local<TElem,TFVGeom> (geo, u, locJ_tri, locU_tri, elem, vCornerCoords);
        
        if ( boundary )
        {
            geo.reset_jacobian_on_interface(locJ_tri, 3);
            this->template add_jac_A_elem_boundary<TElem,TFVGeom> (geo, locJ_tri, locU_tri, elem, vCornerCoords);
        }
        
        geo.set_jacobian_tri(locJ_tri);
        geo.set_DoF_tag_tri(shiftTag);
        
    }
    if ( roidCheck == ROID_QUADRILATERAL )
    {
        geo.set_local_sol(locU_quad, 4, u, orientation);
        LocalVector jump_quad = geo.set_jump_values(ind, 4);
        
        this->template add_jac_A_elem_local<TElem,TFVGeom> (geo, u, locJ_quad, locU_quad, elem, vCornerCoords);
        
        if ( boundary )
        {
            geo.reset_jacobian_on_interface(locJ_quad, 4);
            this->template add_jac_A_elem_boundary<TElem,TFVGeom> (geo, locJ_quad, locU_quad, elem, vCornerCoords);
        }
        
        geo.set_jacobian_quad(locJ_quad);
        geo.set_DoF_tag_quad(shiftTag);
    }
    
        
}
    
    
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
add_jac_A_elem_local(TFVGeom& geo, const LocalVector& u, LocalMatrix& J, LocalVector& locU, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
    MathMatrix<dim, dim> diffusion;
    MatSet(diffusion, 0);
    MatDiagSet(diffusion, 1.0);
    
    LocalVector source, jump;
    
    get_local_data<TElem,TFVGeom>(geo, u, locU, diffusion, jump, jump, jump);
    
    
//	Diff. Tensor times Gradient
	MathVector<dim> Dgrad;

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

//	Diffusion and Velocity Term
	if(m_imDiffusion.data_given() || m_imVelocity.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		////////////////////////////////////////////////////
		// Diffusive Term
		////////////////////////////////////////////////////
			if(m_imDiffusion.data_given())
			{
			//	DID_CONV_DIFF_FV1_CUTELEM
				number D_diff_flux_sum = 0.0;

			// 	loop shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
				// 	Compute Diffusion Tensor times Gradient
					MatVecMult(Dgrad, diffusion, scvf.global_grad(sh));

				//	Compute flux at IP
					const number D_diff_flux = VecDot(Dgrad, scvf.normal());
					UG_DLOG(DID_CONV_DIFF_FV1_CUTELEM, 2, ">>OCT_DISC_DEBUG: " << "convection_diffusion_fv1.cpp: " << "add_jac_A_elem(): " << "sh # "  << sh << " ; normalSize scvf # " << ip << ": " << VecLength(scvf.normal()) << "; \t from "<< scvf.from() << "; to " << scvf.to() << "; D_diff_flux: " << D_diff_flux << "; scvf.global_grad(sh): " << scvf.global_grad(sh) << std::endl);

				// 	Add flux term to local matrix // HIER MATRIXINDIZES!!!
					UG_ASSERT((scvf.from() < J.num_row_dof(_C_)) && (scvf.to() < J.num_col_dof(_C_)),
							  "Bad local dof-index on element with object-id " << elem->base_object_id()
							  << " with center: " << CalculateCenter(elem, vCornerCoords));

					J(_C_, scvf.from(), _C_, sh) -= D_diff_flux;
					J(_C_, scvf.to()  , _C_, sh) += D_diff_flux;

				//	DID_CONV_DIFF_FV1_CUTELEM
					D_diff_flux_sum += D_diff_flux;
				}

				UG_DLOG(DID_CONV_DIFF_FV1_CUTELEM, 2, "D_diff_flux_sum = " << D_diff_flux_sum << std::endl << std::endl);
			}
            
            
        /////////////////////////////////////////////////////////////////////////////
        // Additional diffusive Term due to jump in solution at the interface
        // u^+ - u^- = jump
        /////////////////////////////////////////////////////////////////////////////
            if ( 0 )
            {
                // 	loop shape functions
                for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                {
                    // 	Compute Diffusion Tensor times Gradient
                    MatVecMult(Dgrad, diffusion, scvf.global_grad(sh));
                    
                    //	Compute flux at IP
                    const number D_diff_flux = jump(_C_,sh) * VecDot(Dgrad, scvf.normal());
                    
                    J(_C_, scvf.from(), _C_, sh) -= D_diff_flux;
                    J(_C_, scvf.to()  , _C_, sh) += D_diff_flux;
                }
                
            }

		////////////////////////////////////////////////////
		// Convective Term
		////////////////////////////////////////////////////
			if(m_imVelocity.data_given())
			{
			//	Add Flux contribution
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					const number D_conv_flux = convShape(ip, sh);

				//	Add flux term to local matrix
					J(_C_, scvf.from(), _C_, sh) += D_conv_flux;
					J(_C_, scvf.to(),   _C_, sh) -= D_conv_flux;
				}
			}

			// no explicit dependency on flux import
		}
	}

	//UG_LOG("Local Matrix is: \n"<<J<<"\n");

////////////////////////////////////////////////////
// Reaction Term (using lumping)
////////////////////////////////////////////////////

	if(m_imReactionRate.data_given())
	{
	// 	loop Sub Control Volume (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);
			
		// 	get associated node
			const int co = scv.node_id();
			
		// 	Add to local matrix
			J(_C_, co, _C_, co) += m_imReactionRate[ip] * scv.volume();
		}
	}
	
//	reaction term does not explicitly depend on the associated unknown function

////////////////////////////////
// Singular sources and sinks
////////////////////////////////

	if (m_sss.valid()) {
		const typename TDomain::position_accessor_type& aaPos = this->domain()->position_accessor();
		const typename TDomain::grid_type& grid = *this->domain()->grid();
		const number time = this->time();
		MathVector<1> out;
		for(size_t i = 0; i < geo.num_scv(); i++) {
			const typename TFVGeom::SCV& scv = geo.scv(i);
			const int co = scv.node_id();
			const number len = m_sss->get_contrib_of_scv((TElem*)elem, (Grid&)grid, aaPos, geo, co, time, out);
			if (len == 0.0) continue;
			out[0] *= len;
			if (out[0] < 0.0)
			// sink
				J(_C_, co, _C_, co) -= out[0];
		}
	}
}
    
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
add_jac_A_elem_boundary(TFVGeom& geo, LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
    double diffusion = 0.0;
    
    if ( !geo.get_element_modus() )
        diffusion = geo.get_diffusion(geo.get_boolian_for_diffusion());
    else
        diffusion = geo.get_diffusion();
        
    std::vector<typename TFVGeom::BF>& vBF = geo.get_boundary_faces();
        
    if ( vBF.size() > 2 )
        UG_THROW("add_def_A_elem(): vBF.size() is greater than 2: " << vBF.size() << "\n");
        
    //	loop integration points
    for(size_t ip = 0; ip < vBF.size(); ++ip)
    {
        typename TFVGeom::BF bf = vBF[ip];
        
        //	loop trial space
        for(size_t sh = 0; sh < bf.num_sh(); ++sh)
        {
            UG_LOG("bf.node_id(): " << bf.node_id() << "\n");
            UG_LOG("bf.global_grad(sh): " << bf.global_grad(sh) << "\n");
            UG_LOG("normal(): " << bf.normal() << "\n");
                
            //	add to local matrix
            J(_C_, bf.node_id(), _C_, sh) += VecDot(bf.global_grad(sh), bf.normal());
            J(_C_, bf.node_id(), _C_, sh) *= diffusion;
        }
    }
        
}
    
    
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
    //	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID,1);

	if(!m_imMassScale.data_given()) return;

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local matrix
		J(_C_, co, _C_, co) += scv.volume() * m_imMassScale[ip];
	}

//	m_imMass part does not explicitly depend on associated unknown function
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
    bool output = false;
    bool output_integral = false;
        
    bool debug = false;
    bool boundary = false;
    bool add = true;
    
    // get finite volume geometry
    //	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID,1);
    const bool bElementIsCut = geo.get_element_modus();
    
    // First call with orientation = 1:
    int orientation = 1;
    geo.set_orientation(orientation);
    
    // necessary for call of 'get_solution_tri' an 'get_solution_quad':
    geo.resize_local_data(u);
    std::vector<double> imSource;
    if ( m_imSource.data_given() ) {
        for ( size_t i = 0; i < 3; ++i )
            imSource.push_back(m_imSource[i]);
    }
    
    // normal assembling if not cut by interface:
    if ( !bElementIsCut )
    {
        LocalVector dummyU;
        LocalIndices ind = d.get_indices();
        dummyU.resize(ind);
        dummyU = 0;
        LocalVector source = geo.set_source(imSource, ind, 3, false);
        
        if ( output )
        {
            for ( size_t i = 0; i < 3; ++i)
                UG_LOG("*** corner" << vCornerCoords[i][0] << " and " << vCornerCoords[i][1] << "\n" );
            UG_LOG("\n" );
        }
        
        this->template add_def_A_elem_local<TElem,TFVGeom> (geo, u, d, dummyU, elem, vCornerCoords);
        
        number intValElem = this->template add_l2error_A_elem<TElem,TFVGeom> (geo, ROID_TRIANGLE, d, u, elem);
        geo.add_to_integral(intValElem);
        if ( output_integral )
            UG_LOG("------------------> usual: integral = " << sqrt(geo.get_integral()) << "\n");
        
        return;
    }
        
    // get data:
    LocalVector& locD_tri  = geo.get_defect_tri();
    LocalVector& locD_quad = geo.get_defect_quad();
    
    LocalVector& locU_tri  = geo.get_solution_tri();
    LocalVector& locU_quad = geo.get_solution_quad();
    
    // reset data:
    locD_tri = 0;
    locD_quad = 0;
    
    // call elem disc twice:
    
    if ( debug ) geo.print_InterfaceIDdata();
    
    LocalIndices ind = d.get_indices();
    
    ReferenceObjectID roidCheck = geo.get_roid();
    if ( roidCheck == ROID_TRIANGLE )
    {
        geo.set_local_sol(locU_tri, 3, u, orientation);
        
        LocalVector jump_tri = geo.set_jump_values(ind, 3);
        LocalVector jump_tri_grad = geo.set_jump_grad_values(ind, 3);
        LocalVector source_tri = geo.set_source(imSource, ind, 3, true);
        
        if ( output ) UG_LOG(" tri 1: orientaten: " << orientation << "\n");
        
        this->template add_def_A_elem_local<TElem,TFVGeom> (geo, u, locD_tri, locU_tri, elem, vCornerCoords);
            
        if ( boundary )
        {
            geo.reset_defect_on_interface(locD_tri, 3);
            this->template add_def_A_elem_boundary<TElem,TFVGeom> (geo, locD_tri, locU_tri, elem, vCornerCoords);
        }
        
        if ( output )
        {
            for ( size_t i = 0; i < 3; ++i)
            UG_LOG("corner" << vCornerCoords[i][0] << " and " << vCornerCoords[i][1] << "\n" );
            UG_LOG("\n" );
        }
        
        number intValElem = this->template add_l2error_A_elem<TElem,TFVGeom> (geo, ROID_TRIANGLE, locD_tri, locU_tri, elem);
        if ( add ) geo.add_to_integral(intValElem);
        
        if ( output_integral ) UG_LOG("------------------> tri1: integral = " << sqrt(geo.get_integral()) << "\n");
        
        geo.set_defect_tri(locD_tri);
        geo.set_DoF_tag_tri(false);
            
    }
    if ( roidCheck == ROID_QUADRILATERAL )
    {
        geo.set_local_sol(locU_quad, 4, u, orientation);
        
        LocalVector jump_quad = geo.set_jump_values(ind, 4);
        LocalVector jump_quad_grad = geo.set_jump_grad_values(ind, 4);
        LocalVector source_quad = geo.set_source(imSource, ind, 4, true);
        
        if ( output ) UG_LOG(" quad 1: orientaten: " << orientation << "\n");
        this->template add_def_A_elem_local<TElem,TFVGeom> (geo, u, locD_quad, locU_quad, elem, vCornerCoords);
            
        if ( boundary )
        {
            geo.reset_defect_on_interface(locD_quad, 4);
            this->template add_def_A_elem_boundary<TElem,TFVGeom> (geo, locD_quad, locU_quad, elem, vCornerCoords);
        }
        
        number intValElem = this->template add_l2error_A_elem<TElem,TFVGeom> (geo, ROID_QUADRILATERAL, locD_quad, locU_quad, elem);
        if ( add ) geo.add_to_integral(intValElem);
        
        if ( output_integral ) UG_LOG("------------------> quad1: integral = " << sqrt(geo.get_integral()) << "\n");
        
        
        geo.set_defect_quad(locD_quad);
        geo.set_DoF_tag_quad(false);
        
    }
    
    
    // Second call with orientation = -1:
    bool shiftTag = geo.get_bScaleDoFs();	// shiftTag = true in case of double DoFs on interface!
    
    orientation *= -1;
    if ( output ) UG_LOG(" ____2: orientaten: " << orientation << "\n");
    
    geo.set_orientation(orientation);
    try{
        geo.update(elem, vCornerCoords, &(this->subset_handler()));
    }UG_CATCH_THROW("ConvectionDiffusionFV1_cutElem::update:"
                        " Cannot update Finite Volume Geometry.");
        
    if ( debug ) geo.print_InterfaceIDdata();
        
    roidCheck = geo.get_roid();
    if ( roidCheck == ROID_TRIANGLE )
    {
        geo.set_local_sol(locU_tri, 3, u, orientation);
        LocalVector jump_tri = geo.set_jump_values(ind, 3);
        LocalVector jump_tri_grad = geo.set_jump_grad_values(ind, 3);
        LocalVector source_tri = geo.set_source(imSource, ind, 3, true);
        
        if ( output ) UG_LOG(" tri 2: orientaten: " << orientation << "\n");
        
        this->template add_def_A_elem_local<TElem,TFVGeom> (geo, u, locD_tri, locU_tri, elem, vCornerCoords);
            
        if ( boundary )
        {
            geo.reset_defect_on_interface(locD_tri, 3);
            this->template add_def_A_elem_boundary<TElem,TFVGeom> (geo, locD_tri, locU_tri, elem, vCornerCoords);
        }
        
        number intValElem = this->template add_l2error_A_elem<TElem,TFVGeom> (geo, ROID_TRIANGLE, locD_tri, locU_tri, elem);
        if ( add ) geo.add_to_integral(intValElem);
        
        if ( output ) UG_LOG("------------------> tri2: integral = " << sqrt(geo.get_integral()) << "\n");
        
        
        geo.set_defect_tri(locD_tri);
        geo.set_DoF_tag_tri(shiftTag);
    }
    if ( roidCheck == ROID_QUADRILATERAL )
    {
        geo.set_local_sol(locU_quad, 4, u, orientation);
        LocalVector jump_quad = geo.set_jump_values(ind, 4);
        LocalVector jump_quad_grad = geo.set_jump_grad_values(ind, 4);
        LocalVector source_quad = geo.set_source(imSource, ind, 4, true);
        
        if ( output ) UG_LOG(" quad 2: orientaten: " << orientation << "\n");
        
        this->template add_def_A_elem_local<TElem,TFVGeom> (geo, u, locD_quad, locU_quad, elem, vCornerCoords);
            
        if ( boundary )
        {
            geo.reset_defect_on_interface(locD_quad, 4);
            this->template add_def_A_elem_boundary<TElem,TFVGeom> (geo, locD_quad, locU_quad, elem, vCornerCoords);
        }
        
        number intValElem = this->template add_l2error_A_elem<TElem,TFVGeom> (geo, ROID_QUADRILATERAL, locD_quad, locU_quad, elem);
        if ( add ) geo.add_to_integral(intValElem);
        
        if ( output_integral ) UG_LOG("------------------> quad2: integral = " << sqrt(geo.get_integral()) << "\n");
        
        geo.set_defect_quad(locD_quad);
        geo.set_DoF_tag_quad(shiftTag);
        
    }
}

    
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
add_def_A_elem_local(TFVGeom& geo, const LocalVector& u, LocalVector& d, LocalVector& locU, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
    
    MathMatrix<dim, dim> diffusion;
    MatSet(diffusion, 0);
    MatDiagSet(diffusion, 1.0);
    
    LocalVector jump, jump_grad, source;
    
    get_local_data<TElem,TFVGeom>(geo, u, locU, diffusion, jump, jump_grad, source);

    
//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

	if(m_imDiffusion.data_given() || m_imVelocity.data_given() || m_imFlux.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		/////////////////////////////////////////////////////
		// Diffusive Term
		/////////////////////////////////////////////////////
			if(m_imDiffusion.data_given())
			{
			//	to compute D \nabla c
				MathVector<dim> Dgrad_c, grad_c;

			// 	compute gradient and shape at ip
				VecSet(grad_c, 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(grad_c, locU(_C_,sh), scvf.global_grad(sh));

			//	scale by diffusion tensor
				MatVecMult(Dgrad_c, diffusion, grad_c);

			// 	Compute flux
				const number diff_flux = VecDot(Dgrad_c, scvf.normal());

			// 	Add to local defect
				d(_C_, scvf.from()) -= diff_flux;
				d(_C_, scvf.to()  ) += diff_flux;
			}

        /////////////////////////////////////////////////////////////////////////////
        // Additional diffusive Term due to jump in solution at the interface
        // u^+ - u^- = jump
        /////////////////////////////////////////////////////////////////////////////
            if ( 1 )
            {
                // scale diffusion by jump in solution:
                //	const double jump = 2.0;
                //	diffusion *= jump;
                
                //	to compute D \nabla c=Id_interface
                MathVector<dim> Dgrad, grad;
                
                // 	compute gradient and shape at ip
                VecSet(grad, 0.0);
                for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
                    VecScaleAppend(grad, jump(_C_,sh), scvf.global_grad(sh));
                
                //	scale by diffusion tensor
                MatVecMult(Dgrad, diffusion, grad);
                
                // 	Compute flux
                const number diff_flux = VecDot(Dgrad, scvf.normal());
                
                // 	Add to local defect
                d(_C_, scvf.from()) -= diff_flux;
                d(_C_, scvf.to()  ) += diff_flux;
                
            }
        }

	}

    /////////////////////////////////////////////////////
    // add rhs during same method!
    // --> in elem_disc_assemble_util, the method 'add_rhs_elem()' adds the local vector otherwise! NOT functional!!
    /////////////////////////////////////////////////////
    
    if ( 1 )
    { //m_imSource.data_given() ) {
        for ( size_t ip = 0; ip < geo.num_scv(); ++ip )
        {
            // get current SCV
            const typename TFVGeom::SCV& scv = geo.scv( ip );
            
            // get associated node
            int co = scv.node_id();
            d(_C_, co) -= source(_C_, co) * scv.volume();
            
            // Add to local rhs
            /*			if ( co > 2 )
             {
             d(_C_, co) -= m_imSource[2] * scv.volume();
             UG_LOG("m_imSource[2] * scv.volume(): " << m_imSource[2]  << "\n");
             UG_LOG("source(_C_, co) * scv.volume(): " << source(_C_, co)  << "\n");
             }
             else
             {
             d(_C_, co) -= m_imSource[co] * scv.volume();
             UG_LOG("m_imSource[co] * scv.volume(): " << m_imSource[co]  << "\n");
             UG_LOG("source(_C_, co) * scv.volume(): " << source(_C_, co)  << "\n");
             }
             
             */
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////
    // Additional source Term due to jump in gradient at the interface
    // (\nabla u^+ - \nabla u^-)\cdot n = h(x) * |n|
    /////////////////////////////////////////////////////////////////////////////
    if ( 1 )
    {
        std::vector<typename TFVGeom::BF>& vBF = geo.get_boundary_faces();
        MathVector<dim> Dgrad;
        VecSet(Dgrad, 0.0);
        
        if ( vBF.size() > 2 )
            UG_THROW("add_def_A_elem(): vBF.size() is greater than 2: " << vBF.size() << "\n");
        //	loop integration points
        for(size_t ip = 0; ip < vBF.size(); ++ip)
        {
            typename TFVGeom::BF bf = vBF[ip];
            // Add to local rhs
            d(_C_, bf.node_id()) += jump_grad(_C_,bf.node_id()) * bf.Vol;
        }
    }

////////////////////////////////
// Singular sources and sinks
////////////////////////////////

	if (m_sss.valid()) {
		const typename TDomain::position_accessor_type& aaPos = this->domain()->position_accessor();
		const typename TDomain::grid_type& grid = *this->domain()->grid();
		const number time = this->time();
		MathVector<1> out;
		for(size_t i = 0; i < geo.num_scv(); i++) {
			const typename TFVGeom::SCV& scv = geo.scv(i);
			const int co = scv.node_id();
			const number len = m_sss->get_contrib_of_scv((TElem*)elem, (Grid&)grid, aaPos, geo, co, time, out);
			if (len == 0.0) continue;
			out[0] *= len;
			if (out[0] > 0.0)
			// source
				d(_C_, co) -= out[0];
			else
			// sink
				d(_C_, co) -= out[0] * locU(_C_, co);
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
add_def_A_elem_boundary(TFVGeom& geo, LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
    double diffusion = 0.0;
    
    if ( !geo.get_element_modus() )
        diffusion = geo.get_diffusion(geo.get_boolian_for_diffusion());
    else
        diffusion = geo.get_diffusion();
    
    std::vector<typename TFVGeom::BF>& vBF = geo.get_boundary_faces();
    MathVector<dim> Dgrad;
    VecSet(Dgrad, 0.0);
    
    UG_LOG("---------- vBF.size(): " << vBF.size() << "\n");
    
    if ( vBF.size() > 2 )
        UG_THROW("add_def_A_elem(): vBF.size() is greater than 2: " << vBF.size() << "\n");
    
    // set solution
    /*	for(size_t sh = 0; sh < geo.num_sh(); ++sh)
        {
        u(_C_, sh) = 1.0;
        u(_C_, sh) -= corners[sh][0]*(2*corners[sh][0] - 1.0);
        u(_C_, sh) -= corners[sh][1]*(2*corners[sh][1] - 1.0);
        }
        */
    //		u(_C_, sh) = 1 - 2*(corners[sh][0]*corners[sh][0] + corners[sh][1]*corners[sh][1]) - (corners[sh][0]+corners[sh][1]);
    //		u(_C_, sh) = corners[sh][1]*(2*corners[sh][1] - 1.0);
    //		u(_C_, sh) = (1.0 - corners[sh][0] - corners[sh][1]) * (1.0 - 2*corners[sh][0] - 2*corners[sh][1]);
    //		u(_C_, sh) = corners[sh][1];
    //		for(size_t sh = 0; sh < geo.num_sh(); ++sh)
    //			u(_C_, sh) = corners[sh][0]*(2*corners[sh][0] - 1.0); //corners[sh][0];
    //	u(_C_, sh) = corners[sh][0]*(2* corners[sh][0]-1);
    
    ////////////////////////////////////////////////////////////////////////////////
    //	NO loop integration points!
    //	/* for(size_t ip = 0; ip < vBF.size(); ++ip) */
    // 	Reason: the length of the normal is already the length of the total face (NOT the scvf!)
    ////////////////////////////////////////////////////////////////////////////////
    
    //	loop integration points
    for(size_t ip = 0; ip < vBF.size(); ++ip)
    {
        typename TFVGeom::BF bf = vBF[ip];
        VecSet(Dgrad, 0.0);
        
        //	loop trial space
        for(size_t sh = 0; sh < bf.num_sh(); ++sh)
        {
            //	Diffusion
            UG_LOG("bf.num_sh(): " << bf.num_sh() << "\n");
            UG_LOG("Dgrad: " << Dgrad << "\n");
            UG_LOG("bf.normal(): " << bf.normal() << "\n");
            UG_LOG("bf.global_grad(ip, sh): " << bf.global_grad(sh) << "\n");
            
            UG_LOG("u(_C_, sh): " << u(_C_, sh) << "\n");
            
            VecScaleAppend(Dgrad, diffusion * u(_C_, sh), bf.global_grad(sh));
        }
        
        //	add to local vector
        d(_C_, bf.node_id()) += VecDot(Dgrad, bf.normal());
        
    }
        
    UG_LOG("---------- end ----------- \n\n");
        
}

    
////////////////////////////////////////////////////////////////////////////////
///
///     methods for cut element error computation via 'add_l2error_A_elem()'
///
////////////////////////////////////////////////////////////////////////////////

template <int dim>
number get_exact_sol_test(MathVector<dim> position)
{
    return sin(2*M_PI*position[0]) + sin(2*M_PI*position[1]);
}
    
template <int dim>
number get_exact_sol_Gangl(MathVector<dim> position)
{
    double kappa_2 = 10.0;
    double dist_x = position[0] - 0.0;
    double dist_y = position[1] - 0.0;
    double sqR = 0.4*0.4;
    
    double sqDist = dist_x*dist_x + dist_y*dist_y;
    double dist = sqrt(sqDist);
    
    double returnValue = -4*kappa_2*kappa_2*sqR*sqDist + 2*sqR*sqR*kappa_2*(2*kappa_2 - 1);
    
    if ( dist >= 0.4 )
        returnValue = -2*kappa_2*sqDist*sqDist;
        
    return returnValue;
}

template <int dim>
MathVector<dim> get_exact_grad_Gangl(MathVector<dim> position)
{
    double kappa_2 = 10.0;
    double dist_x = position[0] - 0.0;
    double dist_y = position[1] - 0.0;
    double sqR = 0.4*0.4;
    
    double sqDist = dist_x*dist_x + dist_y*dist_y;
    double dist = sqrt(sqDist);
    
    double factor = -8*kappa_2*kappa_2*sqR;
    
    if ( dist >= 0.4 )
        factor = -8*kappa_2*sqDist;
        
    MathVector<dim> returnVector;
    returnVector[0] = factor*dist_x;
    returnVector[1] = factor*dist_y;
    
    return returnVector;
}

template <int dim>
MathVector<dim> get_exact_grad_FedkiwEx5(MathVector<dim> position)
{
    double center_x = 0.0;
    double center_y = 0.0;
    double radius = 0.5;
    
    double dist_x = position[0] - center_x;
    double dist_y = position[1] - center_y;
    
    double sqDist = dist_x*dist_x + dist_y*dist_y;
    double dist = sqrt(sqDist);
    
    double absValue = position[0]*position[0] + position[1]* position[1];
    
    MathVector<dim> returnVector;
    returnVector[0] = 0.0;
    returnVector[1] = 0.0;
    
    double factor = 1.0/absValue;
    if ( dist >= radius )
    {
        returnVector[0] = factor*position[0];
        returnVector[1] = factor*position[1];
    }
    
    return returnVector;
    
}
    
template <int dim>
number get_exact_sol_FedkiwEx6(MathVector<dim> position)
{
    double center_x = 0.0;
    double center_y = 0.0;
    double radius = 0.5;
    
    double dist_x = position[0] - center_x;
    double dist_y = position[1] - center_y;
    
    double sqDist = dist_x*dist_x + dist_y*dist_y;
    double dist = sqrt(sqDist);
    
    double returnValue = 0.0;
    
    if ( dist <= radius )
        returnValue = exp(position[0])*cos(position[1]);
    
    return returnValue;
}

template <int dim>
number get_exact_sol_FedkiwEx5(MathVector<dim> position)
{
    double center_x = 0.0;
    double center_y = 0.0;
    
    double dist_x = position[0]-center_x;
    double dist_y = position[1]-center_y;
    
    double sqDist = dist_x*dist_x + dist_y*dist_y;
    double dist = sqrt(sqDist);
    
    double returnValue = 1.0;
        
    if ( dist > 0.5 )
        returnValue = 1.0 + log(2*dist);
    
    return returnValue;
}
    
template <int dim>
number get_exact_sol_FedkiwEx3(MathVector<dim> position)
{
    double center_x = 0.5;
    double center_y = 0.5;
    double radius = 0.25;
        
    double dist_x = position[0] - center_x;
    double dist_y = position[1] - center_y;
    
    double sqDist = dist_x*dist_x + dist_y*dist_y;
    double dist = sqrt(sqDist);
    
    double absValue = position[0]*position[0] + position[1]* position[1];
    
    double returnValue = 0.0;
    
    if ( dist <= radius )
        returnValue = exp(-absValue);
    
    return returnValue;
}
    
template <int dim>
MathVector<dim> get_exact_grad(MathVector<dim> position)
{
        
}
    
//////////////////////////////////////////////////////////////////////
// code see ugbase/lib_disc/function_spaces/integrate.h: evaluate() for
// --> L2ErrorIntegrand (for value)
// --> H1ErrorIntegrand (for gradient): lines 1873-1910
//
//	called bei Integrate() via method 'integrand.values':
//	integrand.values(&(vValue[0]), &(vGlobIP[0]),
//		                 pElem, &vCorner[0], rQuadRule.points(),
//						 &(vJT[0]),
//						 numIP);
//
//////////////////////////////////////////////////////////////////////
template<typename TDomain>
template<typename TElem, typename TFVGeom>
number ConvectionDiffusionFV1_cutElem<TDomain>::
add_l2error_A_elem(TFVGeom& geo, ReferenceObjectID roid, LocalVector& d, const LocalVector& u, GridObject* elem)
{
    bool output = false;
        
   // number integral = 0;
    
    std::vector<MathVector<dim> > vCorner;
    std::vector<MathVector<dim> > vGlobIP;
    std::vector<MathVector<dim> > vLocIP;
    std::vector<MathMatrix<dim, dim> > vJT;
    std::vector<number> vValue;
    std::vector<number> vValueGrad;
    
    QuadType type = GetQuadratureType("best");
    
    const QuadratureRule<dim>& rQuadRule
    = QuadratureRuleProvider<dim>::get(roid, 1, type);
    
    //	get reference element mapping by reference object id
    DimReferenceMapping<dim, dim>& mapping
    = ReferenceMappingProvider::get<dim, dim>(roid);
    
    //	number of integration points
    const size_t numIP = rQuadRule.size();
    
    //	get all corner coordinates
    //	CollectCornerCoordinates(vCorner, *pElem, aaPos, true);
    
    const DimReferenceElement<dim>& rRefElem
    = ReferenceElementProvider::get<dim>(roid);
    
    vCorner.clear();
    // 	remember global position of nodes
    for(size_t i = 0; i < rRefElem.num(0); ++i)
        vCorner.push_back(geo.get_corner(i));
    
    if ( output )
    {
        for ( size_t i = 0; i < rRefElem.num(0); ++i)
            UG_LOG("vCorner" << vCorner[i][0] << " and " << vCorner[i][1] << "\n" );
        UG_LOG("\n" );
    }
    
    //	update the reference mapping for the corners
    mapping.update(vCorner);
    
    //	compute global integration points
    vGlobIP.resize(numIP);
    mapping.local_to_global(&(vGlobIP[0]), rQuadRule.points(), numIP);
    
    if ( output ) UG_LOG("vGlobIP" << vGlobIP[0][0] << " and " << vGlobIP[0][1] << "\n" );
    if ( output ) UG_LOG("\n" );
    
    //	compute local integration points
    vLocIP.resize(numIP);
    for(size_t ip = 0; ip < numIP; ++ip)
        vLocIP[ip] = rQuadRule.points()[ip];
    
    if ( output ) UG_LOG("vLocIP" << vLocIP[0][0] << " and " << vLocIP[0][1] << "\n" );
    if ( output ) UG_LOG("\n" );
    
    
    //	compute transformation matrices
    vJT.resize(numIP);
    mapping.jacobian_transposed(&(vJT[0]), rQuadRule.points(), numIP);
    
    const size_t num_sh = geo.num_scvf();
    
    if ( num_sh != rRefElem.num(0) )
        UG_THROW("wrong number of corners: sh = " << num_sh << ", but rRefElem.num(0) = " << rRefElem.num(0) << "\n");
    
    //	compute integrand values at integration points
    vValue.resize(numIP);
    vValueGrad.resize(numIP);
    
    try
    {
        //	loop all integration points
        for(size_t ip = 0; ip < numIP; ++ip)
        {
            //	compute exact solution at integration point
//            number exactSolIP = get_exact_sol_FedkiwEx5<dim>(vGlobIP[ip]);
            number exactSolIP = get_exact_sol_Gangl<dim>(vGlobIP[ip]);
            
            //	compute exact gradient at integration point
            MathVector<dim> exactGradIP = get_exact_grad_FedkiwEx5<dim>(vGlobIP[ip]);
            
            // 	compute approximated solution at integration point
            number approxSolIP = 0.0;
            MathVector<dim> locTmp; VecSet(locTmp, 0.0);
            
            const typename TFVGeom::SCV& scv = geo.scv(ip);
            
            for(size_t sh = 0; sh < num_sh; ++sh)
            {
                //	add shape fct at ip * value at shape
                approxSolIP += u(_C_,sh) * geo.get_shape(sh, vLocIP[ip], roid);
                
                //	add gradient at ip
                VecScaleAppend(locTmp, u(_C_,sh), scv.local_grad(sh));
            }
            
            //	get squared of difference
            vValue[ip] = (exactSolIP - approxSolIP);
            vValue[ip] *= vValue[ip];
            
            //	compute global gradient
            MathVector<dim> approxGradIP;
            MathMatrix<dim, dim> JTInv;
            Inverse(JTInv, vJT[ip]);
            MatVecMult(approxGradIP, JTInv, locTmp);
            
            // get error of gradient
            vValueGrad[ip] = VecDistanceSq(approxGradIP, exactGradIP);
            
            
        }
        /*		integrand.values(&(vValue[0]), &(vGlobIP[0]),
            pElem, &vCorner[0], rQuadRule.points(),
            &(vJT[0]),
            numIP);
            */
    }
    UG_CATCH_THROW("Unable to compute values of integrand at integration point.");
        
    //	reset contribution of this element
    number intValElem = 0;
    
    //	loop integration points
    for(size_t ip = 0; ip < numIP; ++ip)
    {
        //	get quadrature weight
        const number weightIP = rQuadRule.weight(ip);
        
        //	get determinate of mapping
        const number det = SqrtGramDeterminant(vJT[ip]);
        
        //	add contribution of integration point
        intValElem += vValue[ip] * weightIP * det;
        //		intValElem += vValueGrad[ip] * weightIP * det;
            
    }
    
    //	add to global sum
    
    if ( output ) UG_LOG("added: " <<  intValElem << "\n\n");
    
    return intValElem;
}
    
    
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
add_def_A_expl_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
//	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID,1);

//	reaction rate
	if(m_imReactionRateExpl.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_C_, co) += u(_C_, co) * m_imReactionRateExpl[ip] * scv.volume();
		}
	}

//	reaction
	if(m_imReactionExpl.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_C_, co) += m_imReactionExpl[ip] * scv.volume();
		}
	}

	if(m_imSourceExpl.data_given())
	{
		// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
			// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

			// 	get associated node
			const int co = scv.node_id();

			// 	Add to local rhs
			d(_C_, co) -= m_imSourceExpl[ip] * scv.volume();
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
//	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID,1);

	if(!m_imMassScale.data_given() && !m_imMass.data_given()) return;

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	//	mass value
		number val = 0.0;

	//	multiply by scaling
		if(m_imMassScale.data_given())
			val += m_imMassScale[ip] * u(_C_, co);

	//	add mass
		if(m_imMass.data_given())
			val += m_imMass[ip];

	// 	Add to local defect
		d(_C_, co) += val * scv.volume();
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
    /////////////////////////////////////////////////////
    // add rhs ALLREADY(!) during 'add_def_A_elem_local()' method!
    // --> in elem_disc_assemble_util, the method 'add_rhs_elem()' adds the local vector otherwise! NOT functional!!
    /////////////////////////////////////////////////////
    
    return;
    
    // get finite volume geometry
    //	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
    static const TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID,1);
    
	// loop Sub Control Volumes (SCV)
	if ( m_imSource.data_given() ) {
		for ( size_t ip = 0; ip < geo.num_scv(); ++ip ) {
			// get current SCV
			const typename TFVGeom::SCV& scv = geo.scv( ip );

			// get associated node
			const int co = scv.node_id();

			// Add to local rhs
			d(_C_, co) += m_imSource[ip] * scv.volume();
			//UG_LOG("d(_C_, co) = " << d(_C_, co) << "; \t ip " << ip << "; \t co " << co << "; \t scv_vol " << scv.volume() << "; \t m_imSource[ip] " << m_imSource[ip] << std::endl);
		}
	}

	// loop Sub Control Volumes (SCVF)
	if ( m_imVectorSource.data_given() ) {
		for ( size_t ip = 0; ip < geo.num_scvf(); ++ip ) {
			// get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf( ip );

			// Add to local rhs
			d(_C_, scvf.from()) -= VecDot(m_imVectorSource[ip], scvf.normal() );
			d(_C_, scvf.to()  ) += VecDot(m_imVectorSource[ip], scvf.normal() );
		}
	}
}



////////////////////////////////////////////////////////////////////////////////
//	upwind
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ConvectionDiffusionFV1_cutElem<TDomain>::
set_upwind(SmartPtr<IConvectionShapes<dim> > shapes) {m_spConvShape = shapes;}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
const typename ConvectionDiffusionFV1_cutElem<TDomain>::conv_shape_type&
ConvectionDiffusionFV1_cutElem<TDomain>::
get_updated_conv_shapes(const FVGeometryBase& geo)
{
//	compute upwind shapes for transport equation
//	\todo: we should move this computation into the preparation part of the
//			disc, to only compute the shapes once, reusing them several times.
	if(m_imVelocity.data_given())
	{
	//	get diffusion at ips
		const MathMatrix<dim, dim>* vDiffusion = NULL;
		if(m_imDiffusion.data_given()) vDiffusion = m_imDiffusion.values();

	//	update convection shapes
		if(!m_spConvShape->update(&geo, m_imVelocity.values(), vDiffusion, true))
		{
			UG_LOG("ERROR in 'ConvectionDiffusionFV1_cutElem::get_updated_conv_shapes': "
					"Cannot compute convection shapes.\n");
		}
	}

//	return a const (!!) reference to the upwind
	return *const_cast<const IConvectionShapes<dim>*>(m_spConvShape.get());
}

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void ConvectionDiffusionFV1_cutElem<Domain1d>::
register_all_funcs(bool bHang)
{
    register_func<RegularEdge, DimFV1CutGeometry<dim, dim, InterfaceHandlerLocalDiffusion<dim> > >();
    
/*
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		register_func<RegularEdge, HFV1Geometry<RegularEdge, dim> >();
	}
 */
}
#endif

#ifdef UG_DIM_2
template<>
void ConvectionDiffusionFV1_cutElem<Domain2d>::
register_all_funcs(bool bHang)
{
    register_func<Triangle, DimFV1CutGeometry<dim, dim, InterfaceHandlerLocalDiffusion<dim> > >();
    
/*
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
	}
	else
	{
		register_func<RegularEdge, HFV1Geometry<RegularEdge, dim> >();
		register_func<Triangle, HFV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, HFV1Geometry<Quadrilateral, dim> >();
	}
 */
}
#endif

#ifdef UG_DIM_3
template<>
void ConvectionDiffusionFV1_cutElem<Domain3d>::
register_all_funcs(bool bHang)
{
    register_func<Tetrahedron, DimFV1CutGeometry<dim, dim, InterfaceHandlerLocalDiffusion<dim> > >();
    
 /*

//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
		register_func<Octahedron, FV1Geometry<Octahedron, dim> >();
	}
	else
	{
		register_func<RegularEdge, HFV1Geometry<RegularEdge, dim> >();
		register_func<Triangle, HFV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, HFV1Geometry<Quadrilateral, dim> >();
		register_func<Tetrahedron, HFV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, HFV1Geometry<Prism, dim> >();
		register_func<Pyramid, HFV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, HFV1Geometry<Hexahedron, dim> >();
		register_func<Octahedron, HFV1Geometry<Octahedron, dim> >();
	}
*/
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1_cutElem<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	static const int refDim = reference_element_traits<TElem>::dim;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFVGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TFVGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template fsh_elem_loop<TElem, TFVGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFVGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFVGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFVGeom>);
	this->set_add_def_A_expl_elem_fct(id, &T::template add_def_A_expl_elem<TElem, TFVGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(  id, &T::template add_rhs_elem<TElem, TFVGeom>);

}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ConvectionDiffusionFV1_cutElem<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConvectionDiffusionFV1_cutElem<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionFV1_cutElem<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug

