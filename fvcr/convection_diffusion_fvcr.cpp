/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Christian Wehner
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

#include "convection_diffusion_fvcr.h"

#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/spatial_disc/disc_util/fvcr_geom.h"
#include "lib_disc/spatial_disc/disc_util/hfvcr_geom.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusionFVCR<TDomain>::
ConvectionDiffusionFVCR(const char* functions, const char* subsets)
 : ConvectionDiffusionBase<TDomain>(functions,subsets),
   m_spConvShape(new ConvectionShapesNoUpwind<dim>),
   m_bNonRegularGrid(false)
{
	register_all_funcs(m_bNonRegularGrid);
}

template<typename TDomain>
void ConvectionDiffusionFVCR<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
//	check number
	if(vLfeID.size() != 1)
		UG_THROW("ConvectionDiffusion: Wrong number of functions given. "
				"Need exactly "<<1);

	if(vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::CROUZEIX_RAVIART)
		UG_THROW("ConvectionDiffusion FVCR Scheme only implemented for 1st order.");

//	remember
	m_bNonRegularGrid = bNonRegularGrid;

//	update assemble functions
	register_all_funcs(m_bNonRegularGrid);
}

template<typename TDomain>
bool ConvectionDiffusionFVCR<TDomain>::
use_hanging() const
{
	return true;
}


////////////////////////////////////////////////////////////////////////////////
// Assembling functions
////////////////////////////////////////////////////////////////////////////////


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	if(m_imFlux.data_given())
		UG_THROW("ConvectionDiffusionFVCR: Flux import not implemented.");

	if(	m_imSourceExpl.data_given() ||
		m_imReactionExpl.data_given() ||
		m_imReactionRateExpl.data_given())
		UG_THROW("ConvectionDiffusionFVCR: Explicit terms not implemented.");

//	set local positions
	if(!TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;
		static TFVGeom& geo = GeomProvider<TFVGeom>::get();

		geo.update_local_data();

		m_imDiffusion.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                       	                      geo.num_scvf_ips(), false);
		m_imVelocity.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                      	                      geo.num_scvf_ips(), false);
		m_imSource.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                    	                      geo.num_scv_ips(), false);
		m_imReactionRate.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                      	                      geo.num_scv_ips(), false);
		m_imReaction.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                      	                      geo.num_scv_ips(), false);
		m_imMassScale.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                       	                      geo.num_scv_ips(), false);
		m_imMass.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                       	                      geo.num_scv_ips(), false);
	}

//	check, that upwind has been set
//	if(m_spConvShape.invalid())
//		UG_THROW("ConvectionDiffusion::prep_elem_loop:"
//						" Upwind has not been set.");

//	init upwind for element type
//	if(!m_spConvShape->template set_geometry_type<TFVGeom>())
//		UG_THROW("ConvectionDiffusion::prep_elem_loop:"
//						" Cannot init upwind for element type.");
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get();

	try{
		geo.update(elem, vCornerCoords, &(this->subset_handler()));
	}
	UG_CATCH_THROW("ConvectionDiffusion::prep_elem:"
					" Cannot update Finite Volume Geometry.");

//	set local positions
	if(TFVGeom::usesHangingNodes)
	{
		static const int refDim = TElem::dim;
		m_imDiffusion.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                       	                      geo.num_scvf_ips());
		m_imVelocity.template 	set_local_ips<refDim>(geo.scvf_local_ips(),
		                      	                      geo.num_scvf_ips());
		m_imSource.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                    	                      geo.num_scv_ips());
		m_imReactionRate.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                      	                      geo.num_scv_ips());
		m_imReaction.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                      	                      geo.num_scv_ips());
		m_imMassScale.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                       	                      geo.num_scv_ips());
		m_imMass.template 	set_local_ips<refDim>(geo.scv_local_ips(),
		                       	                      geo.num_scv_ips());
/*		if(m_spConvShape.valid())
			if(!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
				UG_THROW("ConvectionDiffusion::prep_elem_loop:"
								" Cannot init upwind for element type.");*/
	}

	//	set global positions
	const MathVector<dim>* vSCVFip = geo.scvf_global_ips();
	const size_t numSCVFip = geo.num_scvf_ips();
	const MathVector<dim>* vSCVip = geo.scv_global_ips();
	const size_t numSCVip = geo.num_scv_ips();
	m_imDiffusion.			set_global_ips(vSCVFip, numSCVFip);
	m_imVelocity.			set_global_ips(vSCVFip, numSCVFip);
	// m_imFlux.				set_global_ips(vSCVFip, numSCVFip);
	m_imSource.				set_global_ips(vSCVip, numSCVip);
	// m_imVectorSource.		set_global_ips(vSCVFip, numSCVFip);
	m_imReactionRate.		set_global_ips(vSCVip, numSCVip);
	// m_imReactionRateExpl.	set_global_ips(vSCVip, numSCVip);
	// m_imReactionExpl.		set_global_ips(vSCVip, numSCVip);
	// m_imSourceExpl.			set_global_ips(vSCVip, numSCVip);
	m_imReaction.			set_global_ips(vSCVip, numSCVip);
	// m_imMassScale.			set_global_ips(vSCVip, numSCVip);
	m_imMass.				set_global_ips(vSCVip, numSCVip);
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	Diff. Tensor times Gradient
	MathVector<dim> Dgrad;

//	get conv shapes
//	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

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
			// 	loop shape functions
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				{
				// 	Compute Diffusion Tensor times Gradient
					MatVecMult(Dgrad, m_imDiffusion[ip], scvf.global_grad(sh));

				//	Compute flux at IP
					const number D_diff_flux = VecDot(Dgrad, scvf.normal());

				// 	Add flux term to local matrix
					J(_C_, scvf.from(), _C_, sh) -= D_diff_flux;
					J(_C_, scvf.to()  , _C_, sh) += D_diff_flux;
				}
			}

		////////////////////////////////////////////////////
		// Convective Term
		////////////////////////////////////////////////////
	/*		if(m_imVelocity.data_given())
			{
			//	Add Flux contribution
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
				{
					const number D_conv_flux = convShape(ip, sh);

				//	Add flux term to local matrix
					J(_C_, scvf.from(), _C_, sh) += D_conv_flux;
					J(_C_, scvf.to(),   _C_, sh) -= D_conv_flux;
				}
			}*/
		}
	}

	// handle constrained dofs
	if(TFVGeom::usesHangingNodes){
		for (size_t i=0;i<geo.num_constrained_dofs();i++){
			const typename TFVGeom::CONSTRAINED_DOF& cd = geo.constrained_dof(i);
			const size_t index = cd.index();
			J(_C_,index,_C_,index) = 1;
			for (size_t j=0;j<cd.num_constraining_dofs();j++)
				J(_C_, index, _C_, cd.constraining_dofs_index(j)) = -cd.constraining_dofs_weight(j);
			// insert interpolation equation directly for all dofs
			for (size_t j=0;j<geo.num_scv();j++){
				const size_t nodeID = geo.scv(j).node_id();
				number alpha=J(_C_,nodeID,_C_,index);
				J(_C_,nodeID,_C_,index)=0;
				for (size_t k=0;k<cd.num_constraining_dofs();k++)
					J(_C_,nodeID,_C_,cd.constraining_dofs_index(k)) += alpha*cd.constraining_dofs_weight(k);
			}
		}
	}

/*	UG_LOG("Local Matrix is: \n");
	size_t n = geo.num_scv()+geo.num_constrained_dofs();
	UG_LOG("locA=[");
	for (size_t i=0;i<n;i++){
		for (size_t j=0;j<n;j++)
			UG_LOG(J(0,i,0,j) << " ");
		UG_LOG(";\n");
	}
	UG_LOG("]\n");*/
////////////////////////////////////////////////////
// Reaction Term (using lumping)
////////////////////////////////////////////////////

//	if no data for reaction rate given, return
	if(!m_imReactionRate.data_given()) return;

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

//	reaction term does not explicitly depend on the associated unknown function
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

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
void ConvectionDiffusionFVCR<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
//	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

	if(m_imDiffusion.data_given() || m_imVelocity.data_given())
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
					VecScaleAppend(grad_c, u(_C_,sh), scvf.global_grad(sh));

			//	scale by diffusion tensor
				MatVecMult(Dgrad_c, m_imDiffusion[ip], grad_c);

			// 	Compute flux
				const number diff_flux = VecDot(Dgrad_c, scvf.normal());

			// 	Add to local defect
				d(_C_, scvf.from()) -= diff_flux;
				d(_C_, scvf.to()  ) += diff_flux;
			}

		/////////////////////////////////////////////////////
		// Convective Term
		/////////////////////////////////////////////////////
	/*		if(m_imVelocity.data_given())
			{
			//	sum up convective flux using convection shapes
				number conv_flux = 0.0;
				for(size_t sh = 0; sh < convShape.num_sh(); ++sh)
					conv_flux += u(_C_, sh) * convShape(ip, sh);

			//  add to local defect
				d(_C_, scvf.from()) += conv_flux;
				d(_C_, scvf.to()  ) -= conv_flux;
			}*/
		}
	}

//	reaction rate
	if(m_imReactionRate.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_C_, co) += u(_C_, co) * m_imReactionRate[ip] * scv.volume();
		}
	}

//	reaction rate
	if(m_imReaction.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t ip = 0; ip < geo.num_scv(); ++ip)
		{
		// 	get current SCV
			const typename TFVGeom::SCV& scv = geo.scv(ip);

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_C_, co) += m_imReaction[ip] * scv.volume();
		}
	}

	// handle constrained dofs, compute defect of interpolation equation
	if(TFVGeom::usesHangingNodes){
		for (size_t i=0;i<geo.num_constrained_dofs();i++){
			const typename TFVGeom::CONSTRAINED_DOF& cd = geo.constrained_dof(i);
			const size_t index = cd.index();
			number defect = u(_C_,index);
			for (size_t j=0;j<cd.num_constraining_dofs();j++)
				defect -= cd.constraining_dofs_weight(j) * u(_C_,cd.constraining_dofs_index(j));
			d(_C_,index) = defect;
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

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
void ConvectionDiffusionFVCR<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	if zero data given, return
	if(!m_imSource.data_given()) return;

// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local rhs
		d(_C_, co) += m_imSource[ip] * scv.volume();
	}
}


//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
lin_def_velocity(const LocalVector& u,
                     std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                     const size_t nip)
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
//	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	//	sum up contributions of convection shapes
		MathVector<dim> linDefect;
		VecSet(linDefect, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
	//		VecScaleAppend(linDefect, u(_C_,sh), convShape.D_vel(ip, sh));

	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] += linDefect;
		vvvLinDef[ip][_C_][scvf.to()] -= linDefect;
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
lin_def_diffusion(const LocalVector& u,
                      std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
                      const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
//	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
	{
	// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

	// 	compute gradient at ip
		MathVector<dim> grad_u;	VecSet(grad_u, 0.0);
		for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			VecScaleAppend(grad_u, u(_C_,sh), scvf.global_grad(sh));

	//	compute the lin defect at this ip
		MathMatrix<dim,dim> linDefect;

	//	part coming from -\nabla u * \vec{n}
		for(size_t k=0; k < (size_t)dim; ++k)
			for(size_t j = 0; j < (size_t)dim; ++j)
				linDefect(j,k) = (scvf.normal())[j] * grad_u[k];

	//	add contribution from convection shapes
	//	if(convShape.non_zero_deriv_diffusion())
		//	for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
			//	MatAdd(linDefect, convShape.D_diffusion(ip, sh), u(_C_, sh));

	//	add contributions
		vvvLinDef[ip][_C_][scvf.from()] -= linDefect;
		vvvLinDef[ip][_C_][scvf.to()  ] += linDefect;
	}
}

//	computes the linearized defect w.r.t to the reaction rate
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
lin_def_reaction_rate(const LocalVector& u,
                          std::vector<std::vector<number> > vvvLinDef[],
                          const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	set lin defect
		vvvLinDef[ip][_C_][co] = u(_C_, co) * scv.volume();
	}
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
lin_def_reaction(const LocalVector& u,
                     std::vector<std::vector<number> > vvvLinDef[],
                     const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	set lin defect
		vvvLinDef[ip][_C_][co] = scv.volume();
	}
}

//	computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
lin_def_source(const LocalVector& u,
                   std::vector<std::vector<number> > vvvLinDef[],
                   const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	// 	set lin defect
		vvvLinDef[ip][_C_][co] = scv.volume();
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
lin_def_mass_scale(const LocalVector& u,
                       std::vector<std::vector<number> > vvvLinDef[],
                       const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t co = 0; co < geo.num_scv(); ++co)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(co);

	// 	Check associated node
		UG_ASSERT(co == scv.node_id(), "Only one shape per SCV");

	// 	set lin defect
		vvvLinDef[co][_C_][co] = u(_C_, co) * scv.volume();
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
lin_def_mass(const LocalVector& u,
                       std::vector<std::vector<number> > vvvLinDef[],
                       const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t co = 0; co < geo.num_scv(); ++co)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(co);

	// 	Check associated node
		UG_ASSERT(co == scv.node_id(), "Only one shape per SCV");

	// 	set lin defect
		vvvLinDef[co][_C_][co] = scv.volume();
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
ex_value(number vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         const LocalVector& u,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<TFVGeom::dim> vLocIP[],
         const size_t nip,
         bool bDeriv,
         std::vector<std::vector<number> > vvvDeriv[])
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;

//	CRFV SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				vValue[ip] += u(_C_, sh) * scvf.shape(sh);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = scvf.shape(sh);
		}
	}
//	CRFV SCV ip
	else if(vLocIP == geo.scv_local_ips())
	{
	//	solution at ip
		for(size_t sh = 0; sh < numSH; ++sh)
			vValue[sh] = u(_C_, sh);

	//	set derivatives if needed
		if(bDeriv)
			for(size_t sh = 0; sh < numSH; ++sh)
				for(size_t sh2 = 0; sh2 < numSH; ++sh2)
					vvvDeriv[sh][_C_][sh2] = (sh==sh2) ? 1.0 : 0.0;
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type> rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		number vShape[numSH];

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.shapes(vShape, vLocIP[ip]);

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < numSH; ++sh)
				vValue[ip] += u(_C_, sh) * vShape[sh];

		//	compute derivative w.r.t. to unknowns iff needed
		//	\todo: maybe store shapes directly in vvvDeriv
			if(bDeriv)
				for(size_t sh = 0; sh < numSH; ++sh)
					vvvDeriv[ip][_C_][sh] = vShape[sh];
		}
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
ex_grad(MathVector<dim> vValue[],
        const MathVector<dim> vGlobIP[],
        number time, int si,
        const LocalVector& u,
        GridObject* elem,
        const MathVector<dim> vCornerCoords[],
        const MathVector<TFVGeom::dim> vLocIP[],
        const size_t nip,
        bool bDeriv,
        std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
// 	Get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	reference dimension
	static const int refDim = ref_elem_type::dim;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;

//	CRFV SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//	Loop Sub Control Volume Faces (SCVF)
		for(size_t ip = 0; ip < geo.num_scvf(); ++ip)
		{
		// 	Get current SCVF
			const typename TFVGeom::SCVF& scvf = geo.scvf(ip);

			VecSet(vValue[ip], 0.0);
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(vValue[ip], u(_C_, sh), scvf.global_grad(sh));

			if(bDeriv)
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = scvf.global_grad(sh);
		}
	}
// 	general case
	else
	{
	//	get trial space
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

	//	storage for shape function at ip
		MathVector<refDim> vLocGrad[numSH];
		MathVector<refDim> locGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(vCornerCoords);

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	evaluate at shapes at ip
			rTrialSpace.grads(vLocGrad, vLocIP[ip]);

		//	compute grad at ip
			VecSet(locGrad, 0.0);
			for(size_t sh = 0; sh < numSH; ++sh)
				VecScaleAppend(locGrad, u(_C_, sh), vLocGrad[sh]);

		//	compute global grad
			mapping.jacobian_transposed_inverse(JTInv, vLocIP[ip]);
			MatVecMult(vValue[ip], JTInv, locGrad);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
				for(size_t sh = 0; sh < numSH; ++sh)
					MatVecMult(vvvDeriv[ip][_C_][sh], JTInv, vLocGrad[sh]);
		}
	}
};

////////////////////////////////////////////////////////////////////////////////
//	upwind
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
void ConvectionDiffusionFVCR<TDomain>::
set_upwind(SmartPtr<IConvectionShapes<dim> > shapes) {m_spConvShape = shapes;}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
const typename ConvectionDiffusionFVCR<TDomain>::conv_shape_type&
ConvectionDiffusionFVCR<TDomain>::
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
			UG_LOG("ERROR in 'ConvectionDiffusionFV1::add_jac_A_elem': "
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
void ConvectionDiffusionFVCR<Domain1d>::
register_all_funcs(bool bHang)
{
	UG_THROW("Crouxeiz-Raviart only senseful in dimension >= 2");
}
#endif

#ifdef UG_DIM_2
template<>
void ConvectionDiffusionFVCR<Domain2d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Triangle, CRFVGeometry<Triangle, dim> >();
		register_func<Quadrilateral, CRFVGeometry<Quadrilateral, dim> >();
	}
	else
	{
		register_func<Triangle, HCRFVGeometry<Triangle, dim> >();
		register_func<Quadrilateral, HCRFVGeometry<Quadrilateral, dim> >();
	}
}
#endif

#ifdef UG_DIM_3
template<>
void ConvectionDiffusionFVCR<Domain3d>::
register_all_funcs(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Tetrahedron, CRFVGeometry<Tetrahedron, dim> >();
		register_func<Prism, CRFVGeometry<Prism, dim> >();
		register_func<Pyramid, CRFVGeometry<Pyramid, dim> >();
		register_func<Hexahedron, CRFVGeometry<Hexahedron, dim> >();
	}
	else
	{
		register_func<Tetrahedron, HCRFVGeometry<Tetrahedron, dim> >();
		register_func<Prism, HCRFVGeometry<Prism, dim> >();
		register_func<Pyramid, HCRFVGeometry<Pyramid, dim> >();
		register_func<Hexahedron, HCRFVGeometry<Hexahedron, dim> >();
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFVCR<TDomain>::
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
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFVGeom>);
	this->set_add_rhs_elem_fct(  id, &T::template add_rhs_elem<TElem, TFVGeom>);

//	set computation of linearized defect w.r.t velocity
	m_imVelocity. set_fct(id, this, &T::template lin_def_velocity<TElem, TFVGeom>);
	m_imDiffusion.set_fct(id, this, &T::template lin_def_diffusion<TElem, TFVGeom>);
	m_imReactionRate. set_fct(id, this, &T::template lin_def_reaction_rate<TElem, TFVGeom>);
	m_imReaction. set_fct(id, this, &T::template lin_def_reaction<TElem, TFVGeom>);
	m_imSource.	  set_fct(id, this, &T::template lin_def_source<TElem, TFVGeom>);
	m_imMassScale.set_fct(id, this, &T::template lin_def_mass_scale<TElem, TFVGeom>);
	m_imMass.	set_fct(id, this, &T::template lin_def_mass<TElem, TFVGeom>);

//	exports
	m_exValue->	   template set_fct<T,refDim>(id, this, &T::template ex_value<TElem, TFVGeom>);
	m_exGrad->template set_fct<T,refDim>(id, this, &T::template ex_grad<TElem, TFVGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ConvectionDiffusionFVCR<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConvectionDiffusionFVCR<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionFVCR<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug

