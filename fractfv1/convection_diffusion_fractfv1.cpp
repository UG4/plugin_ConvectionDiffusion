/*
 * Copyright (c) 2013-2018:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
 * Based on the modules by Andreas Vogel
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

#include "convection_diffusion_fractfv1.h"

#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusionFractFV1<TDomain>::
ConvectionDiffusionFractFV1 (const char* functions, const char* subsets)
:	ConvectionDiffusionBase<TDomain> (functions, subsets),
	m_spConvShape (new ConvectionShapesNoUpwind<dim>)
{
//	initialize the fracture-specific input-export parameters that are not in the base
	m_imAperture.set_comp_lin_defect (false);

//	register the fracture-specific input-export parameters that are not in the base
	this->register_import (m_imAperture);
	this->register_import (m_imOrthoDiffusion);
	this->register_import (m_imOrthoVelocity);
	this->register_import (m_imOrthoFlux);
	this->register_import (m_imOrthoVectorSource);

//	register the element assembling functions
	register_all_funcs ();
}

////////////////////////////////////////////////////////////////////////////////
//	Local discretization interface
////////////////////////////////////////////////////////////////////////////////

/// checks the grid and the shape functions
template<typename TDomain>
void ConvectionDiffusionFractFV1<TDomain>::prepare_setting
(
	const std::vector<LFEID> & vLfeID,
	bool bNonRegular
)
{
//	check the grid
	if (bNonRegular)
		UG_THROW ("ERROR in ConvectionDiffusionFractFV1::prepare_setting:"
				" This discretization does not support hanging nodes.\n");

//	check number of the components
	if (vLfeID.size () != 1)
		UG_THROW ("ConvectionDiffusionFractFV1::prepare_setting: Wrong number of functions given. Need exactly 1.");

//	check whether these are the LagrangeP1 elements
	if (vLfeID[0] != LFEID(LFEID::LAGRANGE, dim, 1))
		UG_THROW ("ConvectionDiffusionFractFV1::prepare_setting: ConvectionDiffusion FV Scheme only implemented for 1st order.");
}

////////////////////////////////////////////////////////////////////////////////
// Assembling functions
////////////////////////////////////////////////////////////////////////////////

///	computes and returns the upwind shapes for the given velocity
template<typename TDomain>
const typename ConvectionDiffusionFractFV1<TDomain>::conv_shape_type&
ConvectionDiffusionFractFV1<TDomain>::get_updated_conv_shapes
(
	bool computeDeriv ///< whether to compute the derivatives of the shapes
)
{
	if(m_imVelocity.data_given())
	{
	//	get diffusion at ips
		const MathMatrix<dim, dim>* vDiffusion = NULL;
		if (m_imDiffusion.data_given ()) vDiffusion = m_imDiffusion.values ();

	//	update convection shapes
		if (!m_spConvShape->update (m_pFractGeo, m_imVelocity.values (), vDiffusion, computeDeriv))
		{
			UG_LOG("ERROR in 'ConvectionDiffusionFractFV1::get_updated_conv_shapes': "
					"Cannot compute convection shapes.\n");
		}
	}

//	return a const (!!) reference to the upwind
	return *const_cast<const IConvectionShapes<dim>*>(m_spConvShape.get());
}

/// prepares the loop over the elements: checks whether the parameters are set, ...
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::prep_elem_loop
(
	ReferenceObjectID roid, ///< only elements with this roid are looped over
	int si ///< and only in this subdomain
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	
//	check the imports
	if (!m_imAperture.data_given())
		UG_THROW ("ConvectionDiffusionFractFV1::prep_elem_loop: Missing Import 'fracture width (aperture)'.");
		
//	check whether we are in a degenerated fracture
	if (! m_spFractManager.valid())
		UG_THROW ("ConvectionDiffusionFractFV1::prep_elem_loop: No fracture manager specified.");
	if (! m_spFractManager->is_closed ())
		UG_THROW ("ConvectionDiffusionFractFV1::prep_elem_loop: Fracture manager not closed.");
	if (! m_spFractManager->contains (si))
		UG_THROW ("ConvectionDiffusionFractFV1::prep_elem_loop:"
			" The fract. discretization is used for subset " << si << " which is not a degenerated fracture.");
	
//	check, that upwind has been set
	if(m_spConvShape.invalid())
		UG_THROW("ConvectionDiffusionFractFV1::prep_elem_loop: Upwind has not been set.");

//	initialize the pointer to the FV geometry for fracture elements
	m_pFractGeo = &GeomProvider<TFractFVGeom>::get (LFEID (LFEID::LAGRANGE, low_dim, 1), 1);
	
//	set up local ip coordinates for corner import parameters
//	REMARK: Note that for the fracture elements, values of the corner import
//	parameters are indexed not by scv (as for the normal elements) but by
//	the corner indices in the reference element
}

template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::fsh_elem_loop ()
{}

template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::prep_elem
(
	const LocalVector & u, ///< local solution at the dofs associated with elem
	GridObject * elem, ///< element to prepare
	ReferenceObjectID roid,  // id of reference element used for assembling
	const position_type vCornerCoords [] ///< coordinates of the corners of the element
)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	
	static const int refDim = TElem::dim; // dimensionality of the element, not the side!
	TElem * pElem = static_cast<TElem*> (elem);
	ref_elem_type& rRefElem = Provider<ref_elem_type>::get ();

//	get the non-degenerated sides of the fracture element
	try
	{
		m_spFractManager->get_layer_sides
			(pElem,
				m_numFractCo, m_innerFractSide, m_innerFractSideIdx, m_innerSideCo,
				m_outerFractSide, m_outerFractSideIdx, m_outerSideCo,
				m_assCo);
	}
	UG_CATCH_THROW("ConvectionDiffusionFractFV1::prep_elem: Cannot find orientation of a fracture element.");
	
//	compute the FV geometry of the inner side
	position_type vSideCornerCoords [maxFractSideCorners]; // global side corner coords
	try
	{
		for (size_t co = 0; co < m_numFractCo; co++)
			vSideCornerCoords [co] = vCornerCoords [m_innerSideCo [co]];
		m_pFractGeo->update (m_innerFractSide, vSideCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("ConvectionDiffusionFractFV1::prep_elem: Cannot update the Finite Volume Geometry for a fracture element.");
	size_t numSCVFip = m_pFractGeo->num_scvf_ips ();
	size_t numSCVip = m_pFractGeo->num_scv_ips ();
	
//	convert local coordinates of the side into the local coordinates of the element (for the input parameters)
	position_type vSideLocCornerCoords [maxFractSideCorners]; // local side corner coords
	position_type vSideLocSCVFipCoords [TFractFVGeom::maxNumSCVF]; // local scvf ip's in a fracture element
	position_type vSideLocSCVipCoords [TFractFVGeom::maxNumSCV]; // local scv ip's in a fracture element
	position_type elemLocCOE; // local coordinates of the mass center of a fracture element
	try
	{
		DimReferenceMapping<low_dim, dim>& rMapping
			= ReferenceMappingProvider::get<low_dim, dim> (m_innerFractSide->reference_object_id ());
		for (size_t co = 0; co < m_numFractCo; co++)
			vSideLocCornerCoords [co] = rRefElem.corner (m_innerSideCo [co]);
		rMapping.update (vSideLocCornerCoords);
		
		rMapping.local_to_global (elemLocCOE, *(m_pFractGeo->coe_local ()));
		rMapping.local_to_global (vSideLocSCVFipCoords, m_pFractGeo->scvf_local_ips (), numSCVFip);
		rMapping.local_to_global (vSideLocSCVipCoords, m_pFractGeo->scv_local_ips (), numSCVip);
	}
	UG_CATCH_THROW("ConvectionDiffusionFractFV1::prep_elem: Cannot transform local side coordinates to local element coordinates in a fracture element.");
	

//	set local positions
	
	m_imDiffusion.template 		set_local_ips<refDim> (vSideLocSCVFipCoords, numSCVFip);
	m_imVelocity.template 		set_local_ips<refDim> (vSideLocSCVFipCoords, numSCVFip);
	m_imFlux.template 			set_local_ips<refDim> (vSideLocSCVFipCoords, numSCVFip);
	m_imVectorSource.template 	set_local_ips<refDim> (vSideLocSCVFipCoords, numSCVFip);
	
	m_imSource.template 		set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	m_imReactionRate.template 	set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	m_imReaction.template 		set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	m_imReactionRateExpl.template 	set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	m_imReactionExpl.template 	set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	m_imSourceExpl.template		set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	m_imMassScale.template 		set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	m_imMass.template 			set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);

	m_imOrthoDiffusion.template set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	m_imOrthoVelocity.template set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	m_imOrthoFlux.template set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	m_imOrthoVectorSource.template set_local_ips<refDim> (vSideLocSCVipCoords, numSCVip);
	
	m_imAperture.template 		set_local_ips<dim> (&elemLocCOE, 1);
	
//	set global positions

	const position_type* vSCVFip = m_pFractGeo->scvf_global_ips ();
	m_imDiffusion.			set_global_ips (vSCVFip, numSCVFip);
	m_imVelocity.			set_global_ips (vSCVFip, numSCVFip);
	m_imFlux.				set_global_ips (vSCVFip, numSCVFip);
	m_imVectorSource.		set_global_ips (vSCVFip, numSCVFip);
	
	const position_type* vSCVip = m_pFractGeo->scv_global_ips ();
	m_imSource.				set_global_ips (vSCVip, numSCVip);
	m_imReactionRate.		set_global_ips (vSCVip, numSCVip);
	m_imReactionRateExpl.	set_global_ips (vSCVip, numSCVip);
	m_imReactionExpl.		set_global_ips (vSCVip, numSCVip);
	m_imSourceExpl.			set_global_ips (vSCVip, numSCVip);
	m_imReaction.			set_global_ips (vSCVip, numSCVip);
	m_imMassScale.			set_global_ips (vSCVip, numSCVip);
	m_imMass.				set_global_ips (vSCVip, numSCVip);

	const position_type* coe_global = m_pFractGeo->coe_global ();
	m_imAperture.set_global_ips (coe_global, 1);
	
//	init upwind for element type
	if(m_spConvShape.valid())
		if(!m_spConvShape->template set_geometry_type<TFractFVGeom> (*m_pFractGeo))
			UG_THROW("ConvectionDiffusionFractFV1::prep_elem:"
							" Cannot init upwind for element type.");
}

/// computes the local stiffness matrix
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::add_jac_A_elem
(
	LocalMatrix & J, ///< the local matrix to update
	const LocalVector & u, ///< current approximation of the solution
	GridObject * elem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	TElem * pElem = static_cast<TElem*> (elem);

	this->template fract_add_jac_A_elem<TElem> (J, u, pElem, vCornerCoords);
	this->template fract_bulk_add_jac_A_elem<TElem> (J, u, pElem, vCornerCoords);
}

/// computes the local stiffness matrix on a fracture element
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::fract_add_jac_A_elem
(
	LocalMatrix & J, ///< the local matrix to update
	const LocalVector & u, ///< current approximation of the solution
	TElem * pElem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	number half_fr_width = m_imAperture[0] / 2;
	
	const size_t numSh = m_pFractGeo->num_sh ();
	const size_t numScvf = m_pFractGeo->num_scvf ();
	
//	Diffusion and Velocity Term
	if (m_imDiffusion.data_given () || m_imVelocity.data_given ())
	{
	//	get conv shapes
		const IConvectionShapes<dim>& convShape = get_updated_conv_shapes (false);
	
	// 	loop Sub Control Volume Faces (SCVF)
		for (size_t ip = 0; ip < numScvf; ++ip)
		{
		// 	get current SCVF
			const typename TFractFVGeom::SCVF& scvf = m_pFractGeo->scvf (ip);

		////////////////////////////////////////////////////
		// Diffusive Term
		////////////////////////////////////////////////////
			if (m_imDiffusion.data_given ())
			{
			//	Diff. Tensor times Gradient
				MathVector<dim> Dgrad;
				
			// 	loop shape functions
				for (size_t sh = 0; sh < numSh; ++sh)
				{
				// 	Compute Diffusion Tensor times Gradient
					MatVecMult (Dgrad, m_imDiffusion[ip], scvf.global_grad (sh));

				//	Compute flux at IP
					const number D_diff_flux = VecDot (Dgrad, scvf.normal ()) * half_fr_width;

					J(_C_, m_innerSideCo[scvf.from()], _C_, m_innerSideCo[sh]) -= D_diff_flux;
					J(_C_, m_innerSideCo[scvf.to()  ], _C_, m_innerSideCo[sh]) += D_diff_flux;
				}
			}

		////////////////////////////////////////////////////
		// Convective Term
		////////////////////////////////////////////////////
			if (m_imVelocity.data_given ())
			{
			//	Add Flux contribution
				for (size_t sh = 0; sh < numSh; ++sh)
				{
					const number D_conv_flux = convShape (ip, sh) * half_fr_width;

				//	Add flux term to local matrix
					J(_C_, m_innerSideCo[scvf.from()], _C_, m_innerSideCo[sh]) += D_conv_flux;
					J(_C_, m_innerSideCo[scvf.to()  ], _C_, m_innerSideCo[sh]) -= D_conv_flux;
				}
			}

			// no explicit dependency on flux import
		}
	}

//	Reaction rate
	if(m_imReactionRate.data_given())
	{
	// 	loop Sub Control Volume (SCV)
		for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);
			
		// 	get associated node
			const int co = m_innerSideCo [scv.node_id ()];
			
		// 	Add to local matrix
			J(_C_, co, _C_, co) += m_imReactionRate[ip] * scv.volume () * half_fr_width;
		}
	}
	
//	Reaction term does not explicitly depend on the associated unknown function
}

/// computes the local stiffness matrix of the fracture-bulk interaction terms on a fracture element
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::fract_bulk_add_jac_A_elem
(
	LocalMatrix & J, ///< the local matrix to update
	const LocalVector & u, ///< current approximation of the solution
	TElem * pElem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	number half_fr_width = m_imAperture[0] / 2;
	
//	loop over the corners of the inner side
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ip++)
	{
	// 	Get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv(ip);
		number s = scv.volume ();

	// 	Get associated node of the element (not side!)
		const int co = m_innerSideCo [scv.node_id()];
	
	//	Assemble Diffusion and Convection
		number D_flux = 0, D_flux_fr = 0;
		
		if (m_imOrthoDiffusion.data_given ())
		{
			D_flux = m_imOrthoDiffusion[ip] / half_fr_width;
			D_flux_fr = - D_flux;
		}
        
        if (m_imOrthoVelocity.data_given ())
        {
			/* We use the full upwind here: */
			if (m_imOrthoVelocity[ip] >= 0)
				D_flux_fr -= m_imOrthoVelocity[ip];
			else
				D_flux -= m_imOrthoVelocity[ip];
        }
		
		J(_C_, m_assCo [co], _C_, m_assCo [co]) += D_flux * s;
		J(_C_, m_assCo [co], _C_, co          ) += D_flux_fr * s;
		
		J(_C_, co, _C_, m_assCo [co]) -= D_flux * s;
		J(_C_, co, _C_, co          ) -= D_flux_fr * s;
	}
}


/// computes the mass matrix of a time-dependent problem
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::add_jac_M_elem
(
	LocalMatrix & J, ///< the local matrix to update
	const LocalVector & u, ///< current approximation of the solution
	GridObject * elem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	if (!m_imMassScale.data_given ()) return;

	number half_fr_width = m_imAperture[0] / 2;
	
// 	loop Sub Control Volumes (SCV)
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);

	// 	get associated node
		const int co = m_innerSideCo [scv.node_id ()];

	// 	Add to local matrix
		J(_C_, co, _C_, co) += scv.volume() * m_imMassScale[ip] * half_fr_width;
	}

//	m_imMass part does not explicitly depend on associated unknown function
}
	
/// computes the stiffness part of the local defect
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::add_def_A_elem
(
	LocalVector & d, ///< the local defect to update
	const LocalVector & u, ///< current approximation of the solution
	GridObject * elem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	TElem * pElem = static_cast<TElem*> (elem);

	this->template fract_add_def_A_elem<TElem> (d, u, pElem, vCornerCoords);
	this->template fract_bulk_add_def_A_elem<TElem> (d, u, pElem, vCornerCoords);
}

/// computes the stiffness part of the local defect on a fracture element
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::fract_add_def_A_elem
(
	LocalVector & d, ///< the local defect to update
	const LocalVector & u, ///< current approximation of the solution
	TElem * pElem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	number half_fr_width = m_imAperture[0] / 2;

//	Diff. Tensor times Gradient
	MathVector<dim> Dgrad;

    const size_t numSh = m_pFractGeo->num_sh ();
	const size_t numScvf = m_pFractGeo->num_scvf ();

//	Diffusion and Velocity Term
	if (m_imDiffusion.data_given () || m_imVelocity.data_given () || m_imFlux.data_given ())
	{
	//	get conv shapes
		const IConvectionShapes<dim>& convShape = get_updated_conv_shapes (false);

	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < numScvf; ++ip)
		{
		// 	get current SCVF
			const typename TFractFVGeom::SCVF& scvf = m_pFractGeo->scvf (ip);

		////////////////////////////////////////////////////
		// Diffusive Term
		////////////////////////////////////////////////////
			if(m_imDiffusion.data_given ())
			{
			//	to compute D \nabla c
				MathVector<dim> Dgrad_c, grad_c;

			// 	compute gradient and shape at ip
				VecSet (grad_c, 0.0);
				for(size_t sh = 0; sh < numSh; ++sh)
					VecScaleAppend (grad_c,
						u (_C_, m_innerSideCo[sh]), scvf.global_grad (sh));

			//	scale by diffusion tensor
				MatVecMult (Dgrad_c, m_imDiffusion[ip], grad_c);

			// 	Compute flux
				const number diff_flux = VecDot (Dgrad_c, scvf.normal ()) * half_fr_width;

			// 	Add to local defect
				d(_C_, m_innerSideCo [scvf.from()]) -= diff_flux;
				d(_C_, m_innerSideCo [scvf.to()  ]) += diff_flux;
			}

		////////////////////////////////////////////////////
		// Convective Term
		////////////////////////////////////////////////////
			if(m_imVelocity.data_given())
			{
			//	sum up convective flux using convection shapes
				number conv_flux = 0.0;
				for(size_t sh = 0; sh < numSh; ++sh)
					conv_flux += u (_C_, m_innerSideCo[sh]) * convShape (ip, sh);
				conv_flux *= half_fr_width;

			//  add to local defect
				d(_C_, m_innerSideCo [scvf.from()]) += conv_flux;
				d(_C_, m_innerSideCo [scvf.to()  ]) -= conv_flux;
			}

		/////////////////////////////////////////////////////
		// Flux Term
		/////////////////////////////////////////////////////
			if(m_imFlux.data_given())
			{
			//	sum up flux
				const number flux = VecDot (m_imFlux[ip], scvf.normal ()) * half_fr_width;

			//  add to local defect
				d(_C_, m_innerSideCo [scvf.from()]) += flux;
				d(_C_, m_innerSideCo [scvf.to()  ]) -= flux;
			}
		}
	}

//	Reaction rate
	if (m_imReactionRate.data_given ())
	{
	// 	loop Sub Control Volumes (SCV)
		for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);

		// 	get associated node
			const int co = m_innerSideCo [scv.node_id ()];

		// 	Add to local defect
			d(_C_, co) += u(_C_, co) * m_imReactionRate[ip] * scv.volume () * half_fr_width;
		}
	}

//	Reaction term
	if (m_imReaction.data_given ())
	{
	// 	loop Sub Control Volumes (SCV)
		for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);

		// 	get associated node
			const int co = m_innerSideCo [scv.node_id ()];

		// 	Add to local defect
			d(_C_, co) += m_imReaction[ip] * scv.volume () * half_fr_width;
		}
	}
}

/// computes the stiffness fracture-bulk interaction terms of the local defect on a fracture element
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::fract_bulk_add_def_A_elem
(
	LocalVector & d, ///< the local defect to update
	const LocalVector & u, ///< current approximation of the solution
	TElem * pElem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	number orthC_f, orthC_m;
	
	number half_fr_width = m_imAperture[0] / 2;
	
//	loop over the corners of the inner side
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ip++)
	{
	// 	Get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv(ip);

	// 	Get associated node of the element (not side!)
		const int co = m_innerSideCo [scv.node_id()];
		
	//	Get the corner values
		orthC_f = u(_C_, co);
		orthC_m = u(_C_, m_assCo[co]);
	
	//	Assemble Diffusion, Convection and the Flux
        number flux = 0;
		
        if (m_imOrthoDiffusion.data_given ())
        	flux = m_imOrthoDiffusion[ip] * (orthC_m - orthC_f) / half_fr_width;
        
        if (m_imOrthoVelocity.data_given ())
        {
			number orthVelocity = m_imOrthoVelocity[ip];
			/* We use the full upwind here: */
			flux -= orthVelocity * ((orthVelocity >= 0)? orthC_f : orthC_m);
		}
        
		if (m_imOrthoFlux.data_given ())
			flux -= m_imOrthoFlux[ip];
        
        flux *= scv.volume ();
        d(_C_, m_assCo[co]) += flux;
        d(_C_, co) -= flux;
	}
}

/// computes the stiffness part of the local defect
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::add_def_A_expl_elem
(
	LocalVector & d, ///< the local defect to update
	const LocalVector & u, ///< current approximation of the solution
	GridObject * elem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	TElem * pElem = static_cast<TElem*> (elem);

	this->template fract_add_def_A_expl_elem<TElem> (d, u, pElem, vCornerCoords);
	this->template fract_bulk_add_def_A_expl_elem<TElem> (d, u, pElem, vCornerCoords);
}

/// computes the stiffness part of the local defect on a fracture element
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::fract_add_def_A_expl_elem
(
	LocalVector & d, ///< the local defect to update
	const LocalVector & u, ///< current approximation of the solution
	TElem * pElem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	number half_fr_width = m_imAperture[0] / 2;

//	Reaction rate
	if (m_imReactionRateExpl.data_given ())
	{
	// 	loop Sub Control Volumes (SCV)
		for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);

		// 	get associated node
			const int co = m_innerSideCo [scv.node_id ()];

		// 	Add to local defect
			d(_C_, co) += u(_C_, co) * m_imReactionRateExpl[ip] * scv.volume () * half_fr_width;
		}
	}

//	reaction
	if (m_imReactionExpl.data_given ())
	{
	// 	loop Sub Control Volumes (SCV)
		for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);

		// 	get associated node
			const int co = m_innerSideCo [scv.node_id ()];

		// 	Add to local defect
			d(_C_, co) += m_imReactionExpl[ip] * scv.volume () * half_fr_width;
		}
	}

//	source
	if (m_imSourceExpl.data_given ())
	{
		// 	loop Sub Control Volumes (SCV)
		for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
		{
			// 	get current SCV
			const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);

			// 	get associated node
			const int co = m_innerSideCo [scv.node_id ()];

			// 	Add to local rhs
			d(_C_, co) -= m_imSourceExpl[ip] * scv.volume () * half_fr_width;
		}
	}
}

/// computes the stiffness fracture-bulk interaction terms of the local defect on a fracture element
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::fract_bulk_add_def_A_expl_elem
(
	LocalVector & d, ///< the local defect to update
	const LocalVector & u, ///< current approximation of the solution
	TElem * pElem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
}

/// computes the mass part of the defect of a time-dependent problem
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::add_def_M_elem
(
	LocalVector & d, ///< the local defect to update
	const LocalVector & u, ///< current approximation of the solution
	GridObject * elem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	if (!m_imMassScale.data_given () && !m_imMass.data_given ()) return;

	number half_fr_width = m_imAperture[0] / 2;
	
// 	loop Sub Control Volumes (SCV)
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);

	// 	get associated node
		const int co = m_innerSideCo [scv.node_id ()];

	//	mass value
		number val = 0;

	//	multiply by scaling
		if (m_imMassScale.data_given ())
			val += m_imMassScale[ip] * u(_C_, co);

	//	add mass
		if (m_imMass.data_given ())
			val += m_imMass[ip];

	// 	Add to local defect
		d(_C_, co) += val * scv.volume () * half_fr_width;
	}
}

/// computes the right-hand side due to the sources
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::add_rhs_elem
(
	LocalVector & d, ///< the right-hand side to update
	GridObject * elem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	TElem * pElem = static_cast<TElem*> (elem);

	this->template fract_add_rhs_elem<TElem> (d, pElem, vCornerCoords);
	this->template fract_bulk_add_rhs_elem<TElem> (d, pElem, vCornerCoords);
}

/// computes the right-hand side due to the sources on a fracture element
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::fract_add_rhs_elem
(
	LocalVector & d, ///< the local defect to update
	TElem * pElem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	number half_fr_width = m_imAperture[0] / 2;

	if (m_imSource.data_given ())
		// loop Sub Control Volumes (SCV)
		for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
		{
			// get current SCV
			const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);

			// get associated node
			const int co = m_innerSideCo [scv.node_id ()];

			// Add to local rhs
			d(_C_, co) += m_imSource[ip] * scv.volume () * half_fr_width;
		}

	if (m_imVectorSource.data_given ())
		// loop Sub Control Volume Fraces (SCVF)
		for(size_t ip = 0; ip < m_pFractGeo->num_scvf (); ++ip)
		{
			// get current SCVF
			const typename TFractFVGeom::SCVF& scvf = m_pFractGeo->scvf (ip);
			
			// compute the flux
			number flux = VecDot (m_imVectorSource[ip], scvf.normal ()) * half_fr_width;

			// Add to local rhs
			d(_C_, m_innerSideCo [scvf.from()]) -= flux;
			d(_C_, m_innerSideCo [scvf.to()  ]) += flux;
		}
}

/// computes the fracture-bulk interaction terms of the local rhs on a fracture element
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::fract_bulk_add_rhs_elem
(
	LocalVector & d, ///< the local defect to update
	TElem * pElem, ///< grid element
	const position_type vCornerCoords [] ///< corner coordinates of the grid element
)
{
	if (m_imOrthoVectorSource.data_given ())
		// loop Sub Control Volumes (SCV)
		for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ip++)
		{
		// 	Get current SCV
			const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv(ip);

		// 	Get associated node of the element (not side!)
			const int co = m_innerSideCo [scv.node_id()];
		
			number flux = m_imOrthoVectorSource[ip] * scv.volume ();
			d(_C_, m_assCo[co]) += flux;
			d(_C_, co) -= flux;
		}
}	

///	computes the linearized defect w.r.t to the fracture velocity
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_velocity
(
	const LocalVector& u,
	std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes (true);

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < m_pFractGeo->num_scvf (); ++ip)
	{
	// get current SCVF
		const typename TFractFVGeom::SCVF& scvf = m_pFractGeo->scvf (ip);

	//	sum up contributions of convection shapes
		MathVector<dim> linDefect;
		VecSet(linDefect, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend (linDefect, u (_C_, m_innerSideCo[sh]), convShape.D_vel (ip, sh));

	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][m_innerSideCo [scvf.from()]] += linDefect;
		vvvLinDef[ip][_C_][m_innerSideCo [scvf.to()]  ] -= linDefect;
	}
}

///	computes the linearized defect w.r.t to the orthogonal velocity
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_ortho_velocity
(
	const LocalVector& u,
	std::vector<std::vector<number> > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ip++)
	{
	// 	Get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv(ip);

	// 	Get associated node of the element (not side!)
		const int co = m_innerSideCo [scv.node_id()];
		
	//	Get the corner values
		number orthC_f = u(_C_, co);
		number orthC_m = u(_C_, m_assCo[co]);
	
	//	Compute the derivative of the flux
		number D_flux_vel = - ((m_imOrthoVelocity[ip] >= 0)? orthC_f : orthC_m) * scv.volume ();

	//	add parts for both sides of the fracture
		vvvLinDef[ip][_C_][m_assCo[co]] += D_flux_vel;
		vvvLinDef[ip][_C_][co         ] -= D_flux_vel;
	}
}

///	computes the linearized defect w.r.t to the fracture diffusion
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_diffusion
(
	const LocalVector& u,
	std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes (true);

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < m_pFractGeo->num_scvf (); ++ip)
	{
	// get current SCVF
		const typename TFractFVGeom::SCVF& scvf = m_pFractGeo->scvf (ip);

	// 	compute gradient at ip
		MathVector<dim> grad_c;
		VecSet (grad_c, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend (grad_c, u (_C_, m_innerSideCo[sh]), scvf.global_grad (sh));

	//	compute the lin defect at this ip
		MathMatrix<dim,dim> linDefect;

	//	part coming from $-\nabla u * \vec{n}
		for(size_t k=0; k < (size_t) dim; ++k)
			for(size_t j = 0; j < (size_t)dim; ++j)
				linDefect(j,k) = (scvf.normal())[j] * grad_c[k];

	//	add contribution from convection shapes
		if(convShape.non_zero_deriv_diffusion ())
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				MatAdd (linDefect,
					convShape.D_diffusion (ip, sh), u (_C_, m_innerSideCo[sh]));

	//	add contributions
		vvvLinDef[ip][_C_][m_innerSideCo[scvf.from()]] -= linDefect;
		vvvLinDef[ip][_C_][m_innerSideCo[scvf.to()]  ] += linDefect;
	}
}

///	computes the linearized defect w.r.t to the orthogonal diffusion
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_ortho_diffusion
(
	const LocalVector& u,
	std::vector<std::vector<number> > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

	number orthC_f, orthC_m;
	
	number half_fr_width = m_imAperture[0] / 2;
	
//	loop over the corners of the inner side
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ip++)
	{
	// 	Get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv(ip);

	// 	Get associated node of the element (not side!)
		const int co = m_innerSideCo [scv.node_id()];
		
	//	Get the corner values
		orthC_f = u(_C_, co);
		orthC_m = u(_C_, m_assCo[co]);
	
		number linDefect = (orthC_m - orthC_f) * scv.volume () / half_fr_width;

	//	add contributions
		vvvLinDef[ip][_C_][m_assCo[co]] += linDefect;
		vvvLinDef[ip][_C_][co         ] -= linDefect;
    }
}

///	computes the linearized defect w.r.t to the fracture flux
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_flux
(
	const LocalVector& u,
	std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < m_pFractGeo->num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFractFVGeom::SCVF& scvf = m_pFractGeo->scvf (ip);

	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][m_innerSideCo [scvf.from()]] += scvf.normal();
		vvvLinDef[ip][_C_][m_innerSideCo [scvf.to()  ]] -= scvf.normal();
	}
}

///	computes the linearized defect w.r.t to the orthogonal flux
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_ortho_flux
(
	const LocalVector& u,
	std::vector<std::vector<number> > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//	loop over the corners of the inner side
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ip++)
	{
	// 	Get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv(ip);
	
	// 	Get associated node of the element (not side!)
		const int co = m_innerSideCo [scv.node_id()];
		
	//	add parts for both sides of the fracture
		vvvLinDef[ip][_C_][m_assCo[co]] += scv.volume ();
		vvvLinDef[ip][_C_][co         ] -= scv.volume ();
	}
}

///	computes the linearized defect w.r.t to the reaction source
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_reaction
(
	const LocalVector& u,
	std::vector<std::vector<number> > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

	number half_fr_width = m_imAperture[0] / 2;
	
//	loop over the corners of the inner side
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ip++)
	{
	// 	Get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv(ip);
	
	// 	Get associated node of the element (not side!)
		const int co = m_innerSideCo [scv.node_id()];
		
	//	add parts for both sides of the fracture
		vvvLinDef[ip][_C_][co] = scv.volume () * half_fr_width;
	}
}

///	computes the linearized defect w.r.t to the reaction rate
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_reaction_rate
(
	const LocalVector& u,
	std::vector<std::vector<number> > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

	number half_fr_width = m_imAperture[0] / 2;
	
//	loop over the corners of the inner side
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ip++)
	{
	// 	Get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv(ip);
	
	// 	Get associated node of the element (not side!)
		const int co = m_innerSideCo [scv.node_id()];
		
	//	add parts for both sides of the fracture
		vvvLinDef[ip][_C_][co] = u(_C_, co) * scv.volume () * half_fr_width;
	}
}

///	computes the linearized defect w.r.t to the source
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_source
(
	const LocalVector& u,
	std::vector<std::vector<number> > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

	number half_fr_width = m_imAperture[0] / 2;
	
//	loop over the corners of the inner side
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ip++)
	{
	// 	Get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv(ip);
	
	// 	Get associated node of the element (not side!)
		const int co = m_innerSideCo [scv.node_id()];
		
	//	add parts for both sides of the fracture
		vvvLinDef[ip][_C_][co] = scv.volume () * half_fr_width;
	}
}

///	computes the linearized defect w.r.t to the fracture flux
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_vector_source
(
	const LocalVector& u,
	std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < m_pFractGeo->num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFractFVGeom::SCVF& scvf = m_pFractGeo->scvf (ip);

	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][m_innerSideCo [scvf.from()]] -= scvf.normal();
		vvvLinDef[ip][_C_][m_innerSideCo [scvf.to()  ]] += scvf.normal();
	}
}

///	computes the linearized defect w.r.t to the orthogonal flux
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_ortho_vector_source
(
	const LocalVector& u,
	std::vector<std::vector<number> > vvvLinDef[],
	const size_t nip
)
{
//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//	loop over the corners of the inner side
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ip++)
	{
	// 	Get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv(ip);
	
	// 	Get associated node of the element (not side!)
		const int co = m_innerSideCo [scv.node_id()];
		
	//	add parts for both sides of the fracture
		vvvLinDef[ip][_C_][m_assCo[co]] -= scv.volume ();
		vvvLinDef[ip][_C_][co         ] += scv.volume ();
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_mass_scale
(
	const LocalVector& u,
	std::vector<std::vector<number> > vvvLinDef[],
	const size_t nip
)
{
	number half_fr_width = m_imAperture[0] / 2;
	
// 	loop Sub Control Volumes (SCV)
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);

	// 	get associated node
		const int co = m_innerSideCo [scv.node_id ()];

	// 	set lin defect
		vvvLinDef[co][_C_][co] = u(_C_, co) * scv.volume () * half_fr_width;
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::lin_def_mass
(
	const LocalVector& u,
	std::vector<std::vector<number> > vvvLinDef[],
	const size_t nip
)
{
	number half_fr_width = m_imAperture[0] / 2;

// 	loop Sub Control Volumes (SCV)
	for (size_t ip = 0; ip < m_pFractGeo->num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFractFVGeom::SCV& scv = m_pFractGeo->scv (ip);

	// 	get associated node
		const int co = m_innerSideCo [scv.node_id ()];

	// 	set lin defect
		vvvLinDef[co][_C_][co] = scv.volume () * half_fr_width;
	}
}

/// registers the local assembler functions for a given element
template<typename TDomain>
template<typename TElem>
void ConvectionDiffusionFractFV1<TDomain>::register_func ()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;

	this->clear_add_fct(id);
	
	this->set_prep_elem_loop_fct (id, &T::template prep_elem_loop<TElem>);
	this->set_prep_elem_fct (     id, &T::template prep_elem<TElem>);
	this->set_fsh_elem_loop_fct ( id, &T::template fsh_elem_loop<TElem>);
	this->set_add_jac_A_elem_fct (id, &T::template add_jac_A_elem<TElem>);
	this->set_add_jac_M_elem_fct (id, &T::template add_jac_M_elem<TElem>);
	this->set_add_def_A_elem_fct (id, &T::template add_def_A_elem<TElem>);
	this->set_add_def_A_expl_elem_fct(id, &T::template add_def_A_expl_elem<TElem>);
	this->set_add_def_M_elem_fct (id, &T::template add_def_M_elem<TElem>);
	this->set_add_rhs_elem_fct (  id, &T::template add_rhs_elem<TElem>);

//	set computation of linearized defect w.r.t velocity, diffusion etc.
	m_imDiffusion.set_fct (id, this, &T::template lin_def_diffusion<TElem>);
	m_imOrthoDiffusion.set_fct (id, this, &T::template lin_def_ortho_diffusion<TElem>);
	m_imVelocity.set_fct (id, this, &T::template lin_def_velocity<TElem>);
	m_imOrthoVelocity.set_fct (id, this, &T::template lin_def_ortho_velocity<TElem>);
	m_imFlux.set_fct (id, this, &T::template lin_def_flux<TElem>);
	m_imOrthoFlux.set_fct (id, this, &T::template lin_def_ortho_flux<TElem>);
	m_imReactionRate.set_fct (id, this, &T::template lin_def_reaction_rate<TElem>);
	m_imReaction.set_fct (id, this, &T::template lin_def_reaction<TElem>);
	m_imSource.set_fct (id, this, &T::template lin_def_source<TElem>);
	m_imVectorSource.set_fct (id, this, &T::template lin_def_vector_source<TElem>);
	m_imOrthoVectorSource.set_fct (id, this, &T::template lin_def_ortho_vector_source<TElem>);
	m_imMassScale.set_fct (id, this, &T::template lin_def_mass_scale<TElem>);
	m_imMass.set_fct (id, this, &T::template lin_def_mass<TElem>);
}

///	registers the interface functions for all types of the elements
template<typename TDomain>
void ConvectionDiffusionFractFV1<TDomain>::register_all_funcs ()
{
	typedef typename domain_traits<dim>::DimElemList AssembleElemList;
	
	boost::mpl::for_each<AssembleElemList> (RegisterLocalDiscr (this));
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_2
template class ConvectionDiffusionFractFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionFractFV1<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug

/* End of File */
