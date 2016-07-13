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

#include "convection_diffusion_fe.h"

#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"
#include "lib_disc/local_finite_element/lagrange/lagrangep1.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusionFE<TDomain>::
ConvectionDiffusionFE(const char* functions, const char* subsets)
 : ConvectionDiffusionBase<TDomain>(functions,subsets),
	m_bQuadOrderUserDef(false)
{
	this->clear_add_fct();
}

template<typename TDomain>
void ConvectionDiffusionFE<TDomain>::set_quad_order(size_t order)
{
	m_quadOrder = order;
	m_bQuadOrderUserDef = true;
}

template<typename TDomain>
void ConvectionDiffusionFE<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	//	check number of fcts
	if(vLfeID.size() != 1)
		UG_THROW("ConvectionDiffusion: Wrong number of functions given. "
				"Need exactly "<<1);

	//	check that not ADAPTIVE
	if(vLfeID[0].order() < 1)
		UG_THROW("ConvectionDiffusion: Adaptive order not implemented.");

	//	set order
	m_lfeID = vLfeID[0];
	if(!m_bQuadOrderUserDef) m_quadOrder = 2*m_lfeID.order()+1;

	register_all_funcs(m_lfeID, m_quadOrder);
}

template<typename TDomain>
bool ConvectionDiffusionFE<TDomain>::
use_hanging() const
{
	return false;
}

////////////////////////////////////////////////////////////////////////////////
// Assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	if(	m_imSourceExpl.data_given() ||
		m_imReactionExpl.data_given() ||
		m_imReactionRateExpl.data_given())
		UG_THROW("ConvectionDiffusionFE: Explicit terms not implemented.");

//	request geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	prepare geometry for type and order
	try{
		geo.update_local(roid, m_lfeID, m_quadOrder);
	}UG_CATCH_THROW("ConvectionDiffusion::prep_elem_loop:"
					" Cannot update Finite Element Geometry.");

//	set local positions
	static const int refDim = TElem::dim;
	m_imDiffusion.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imVelocity.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imFlux.template  		set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imSource.template    set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imVectorSource.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imReactionRate.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imReaction.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imMassScale.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imMass.template 	   set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	try{
		geo.update(elem, vCornerCoords);
	}
	UG_CATCH_THROW("ConvectionDiffusion::prep_elem:"
					" Cannot update Finite Element Geometry.");

//	set global positions for rhs
	m_imDiffusion.	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imVelocity. 	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imFlux. 		set_global_ips(geo.global_ips(), geo.num_ip());
	m_imSource.   	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imVectorSource.set_global_ips(geo.global_ips(), geo.num_ip());
	m_imReactionRate.set_global_ips(geo.global_ips(), geo.num_ip());
	m_imReaction. 	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imMassScale.	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imMass.	  	set_global_ips(geo.global_ips(), geo.num_ip());
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	MathVector<dim> v, Dgrad;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop trial space
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
		//	Diffusion
			if(m_imDiffusion.data_given())
				MatVecMult(Dgrad, m_imDiffusion[ip], geo.global_grad(ip, j));
			else
				VecSet(Dgrad, 0.0);

		//  Convection
			if(m_imVelocity.data_given())
				VecScaleAppend(Dgrad, -1*geo.shape(ip,j), m_imVelocity[ip]);

		//	no explicit dependency on m_imFlux

		//	loop test space
			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
			//	compute integrand
				number integrand = VecDot(Dgrad, geo.global_grad(ip, i));

			// 	Reaction
				if(m_imReactionRate.data_given())
					integrand += m_imReactionRate[ip] * geo.shape(ip, j) * geo.shape(ip, i);

			//	no explicit dependency on m_imReaction

			//	multiply by weight
				integrand *= geo.weight(ip);

			//	add to local matrix
				J(_C_, i, _C_, j) += integrand;
			}
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	if(!m_imMassScale.data_given()) return;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test space
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	loop trial space
			for(size_t j= 0; j < geo.num_sh(); ++j)
			{
			//	add to local matrix
				J(_C_, i, _C_, j) += geo.shape(ip, i) *geo.shape(ip, j)
									 * geo.weight(ip) * m_imMassScale[ip];

			//	no explicit dependency on m_imMass
			}
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	number integrand, shape_u;
	MathMatrix<dim,dim> D;
	MathVector<dim> v, Dgrad_u, grad_u;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	// 	get current u and grad_u
		VecSet(grad_u, 0.0);
		shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
			VecScaleAppend(grad_u, u(_C_,j), geo.global_grad(ip, j));
			shape_u += u(_C_,j) * geo.shape(ip, j);
		}

	// 	Diffusion
		if(m_imDiffusion.data_given())
			MatVecMult(Dgrad_u, m_imDiffusion[ip], grad_u);
		else
			VecSet(Dgrad_u, 0.0);

	// 	Convection
		if(m_imVelocity.data_given())
			VecScaleAppend(Dgrad_u, -1*shape_u, m_imVelocity[ip]);

	// 	Convection
		if(m_imFlux.data_given())
			VecScaleAppend(Dgrad_u, 1.0, m_imFlux[ip]);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute integrand
			integrand = VecDot(Dgrad_u, geo.global_grad(ip, i));

		// 	add Reaction Rate
			if(m_imReactionRate.data_given())
				integrand += m_imReactionRate[ip] * shape_u * geo.shape(ip, i);

		// 	add Reaction
			if(m_imReaction.data_given())
				integrand += m_imReaction[ip] * geo.shape(ip, i);

		//	multiply by integration weight
			integrand *= geo.weight(ip);

		//	add to local defect
			d(_C_, i) += integrand;
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	number shape_u = 0.0;

	if(!m_imMassScale.data_given() && !m_imMass.data_given()) return;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	compute value of current solution at ip
		if(m_imMassScale.data_given()){
			shape_u = 0.0;
			for(size_t j = 0; j < geo.num_sh(); ++j)
				shape_u += u(_C_,j) * geo.shape(ip, j);
		}

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			number val = 0.0;

		//	add MassScaling
			if(m_imMassScale.data_given())
				val += shape_u * m_imMassScale[ip];

		//	add Maxx
			if(m_imMass.data_given())
				val += m_imMass[ip];

		//	add to local defect
			d(_C_, i) +=  val * geo.shape(ip, i) * geo.weight(ip);
		}
	}
};

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	skip if no source present
	if(!m_imSource.data_given() && !m_imVectorSource.data_given()) return;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces

		// only do this if (volume) source is given
		if(m_imSource.data_given()) {
			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
			//	add contribution to local defect
				d(_C_, i) += m_imSource[ip] * geo.shape(ip, i) * geo.weight(ip);
			}
		}

		//	only do this if vector source is given
		if(m_imVectorSource.data_given()) {
			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
			//	add contribution to local defect
				d(_C_, i) += geo.weight(ip) * VecDot(m_imVectorSource[ip], geo.global_grad(ip, i));
			}
		}
	}
}

// ////////////////////////////////
//   error estimation (begin)   ///

//	prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{
	//	get the error estimator data object and check that it is of the right type
	//	we check this at this point in order to be able to dispense with this check later on
	//	(i.e. in prep_err_est_elem and compute_err_est_A_elem())
	if (this->m_spErrEstData.get() == NULL)
	{
		UG_THROW("No ErrEstData object has been given to this ElemDisc!");
	}

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (!err_est_data)
	{
		UG_THROW("Dynamic cast to SideAndElemErrEstData failed."
				<< std::endl << "Make sure you handed the correct type of ErrEstData to this discretization.");
	}

//	set local positions
	static const int refDim = TElem::dim;

	// get local IPs
	size_t numSideIPs, numElemIPs;
	const MathVector<refDim>* sideIPs;
	const MathVector<refDim>* elemIPs;
	try
	{
		numSideIPs = err_est_data->num_all_side_ips(roid);
		numElemIPs = err_est_data->num_elem_ips(roid);
		sideIPs = err_est_data->template side_local_ips<refDim>(roid);
		elemIPs = err_est_data->template elem_local_ips<refDim>(roid);

		if (!sideIPs || !elemIPs) return;	// are NULL if TElem is not of the same dim as TDomain
	}
	UG_CATCH_THROW("Integration points for error estimator cannot be set.");

	// set local IPs in imports
	m_imDiffusion.template 		set_local_ips<refDim>(sideIPs, numSideIPs, false);
	m_imVelocity.template 		set_local_ips<refDim>(sideIPs, numSideIPs, false);
	m_imFlux.template 			set_local_ips<refDim>(sideIPs, numSideIPs, false);
	m_imSource.template 		set_local_ips<refDim>(elemIPs, numElemIPs, false);
	m_imVectorSource.template 	set_local_ips<refDim>(sideIPs, numSideIPs, false);
	m_imReactionRate.template 	set_local_ips<refDim>(elemIPs, numElemIPs, false);
	m_imReaction.template 		set_local_ips<refDim>(elemIPs, numElemIPs, false);
	m_imMassScale.template 		set_local_ips<refDim>(elemIPs, numElemIPs, false);
	m_imMass.template 			set_local_ips<refDim>(elemIPs, numElemIPs, false);

	// store values of shape functions in local IPs
	LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace
				= Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

	m_shapeValues.resize(numElemIPs, numSideIPs, trialSpace.num_sh());
	for (size_t ip = 0; ip < numElemIPs; ip++)
		trialSpace.shapes(m_shapeValues.shapesAtElemIP(ip), elemIPs[ip]);
	for (size_t ip = 0; ip < numSideIPs; ip++)
		trialSpace.shapes(m_shapeValues.shapesAtSideIP(ip), sideIPs[ip]);
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

//	request geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	try{
		geo.update(elem, vCornerCoords);
	}
	UG_CATCH_THROW("ConvectionDiffusion::prep_elem:"
					" Cannot update Finite Element Geometry.");

//	roid
	ReferenceObjectID roid = elem->reference_object_id();

//	set global positions
	size_t numSideIPs, numElemIPs;
	const MathVector<dim>* sideIPs;
	const MathVector<dim>* elemIPs;

	try
	{
		numSideIPs = err_est_data->num_all_side_ips(roid);
		numElemIPs = err_est_data->num_elem_ips(roid);

		sideIPs = err_est_data->all_side_global_ips(elem, vCornerCoords);
		elemIPs = err_est_data->elem_global_ips(elem, vCornerCoords);
	}
	UG_CATCH_THROW("Global integration points for error estimator cannot be set.");

	m_imDiffusion.			set_global_ips(&sideIPs[0], numSideIPs);
	m_imVelocity.			set_global_ips(&sideIPs[0], numSideIPs);
	m_imFlux.				set_global_ips(&sideIPs[0], numSideIPs);
	m_imSource.				set_global_ips(&elemIPs[0], numElemIPs);
	m_imVectorSource.		set_global_ips(&sideIPs[0], numSideIPs);
	m_imReactionRate.		set_global_ips(&elemIPs[0], numElemIPs);
	m_imReaction.			set_global_ips(&elemIPs[0], numElemIPs);
	m_imMassScale.			set_global_ips(&elemIPs[0], numElemIPs);
	m_imMass.				set_global_ips(&elemIPs[0], numElemIPs);
}

//	computes the error estimator contribution (stiffness part) for one element
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

//	request geometry
	static const TFEGeom& geo = GeomProvider<TFEGeom>::get();


// SIDE TERMS //

//	get the sides of the element
	//	We have to cast elem to a pointer of type SideAndElemErrEstData::elem_type
	//	for the SideAndElemErrEstData::operator() to work properly.
	//	This cannot generally be achieved by casting to TElem*, since this method is also registered for
	//	lower-dimensional types TElem, and must therefore be compilable, even if it is never EVER to be executed.
	//	The way we achieve this here, is by calling associated_elements_sorted() which has an implementation for
	//	all possible types. Whatever comes out of it is of course complete nonsense if (and only if)
	//	SideAndElemErrEstData::elem_type != TElem. To be on the safe side, we throw an error if the number of
	//	entries in the list is not as it should be.

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::side_type>::secure_container side_list;
	pErrEstGrid->associated_elements_sorted(side_list, (TElem*) elem);
	if (side_list.size() != (size_t) ref_elem_type::numSides)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFE::compute_err_est_elem'");

// 	some help variables
	MathVector<dim> fluxDensity, gradC, normal;

	// FIXME: The computation of the gradient has to be reworked.
	// In the case of P1 shape functions, it is valid. For Q1 shape functions, however,
	// the gradient is not constant (but bilinear) on the element - and along the sides.
	// We cannot use the FVGeom here. Instead, we need to calculate the gradient in each IP!

	// calculate grad u as average (over scvf)
	VecSet(gradC, 0.0);
	for(size_t ii = 0; ii < geo.num_ip(); ++ii)
	{
		for (size_t j=0; j<m_shapeValues.num_sh(); j++)
				VecScaleAppend(gradC, u(_C_,j), geo.global_grad(ii, j));
	}
	VecScale(gradC, gradC, (1.0/geo.num_ip()));

// calculate flux through the sides
	size_t passedIPs = 0;
	for (size_t side=0; side < (size_t) ref_elem_type::numSides; side++)
	{
		// normal on side
		SideNormal<ref_elem_type,dim>(normal, side, vCornerCoords);
		VecNormalize(normal, normal);

		try
		{
			for (size_t sip = 0; sip < err_est_data->num_side_ips(side_list[side]); sip++)
			{
				size_t ip = passedIPs + sip;

				VecSet(fluxDensity, 0.0);

			// diffusion //
				if (m_imDiffusion.data_given())
					MatVecScaleMultAppend(fluxDensity, -1.0, m_imDiffusion[ip], gradC);

			// convection //
				if (m_imVelocity.data_given())
				{
					number val = 0.0;
					for (size_t sh = 0; sh < m_shapeValues.num_sh(); sh++)
						val += u(_C_,sh) * m_shapeValues.shapeAtSideIP(sh,sip);

					VecScaleAppend(fluxDensity, val, m_imVelocity[ip]);
				}

			// general flux //
				if (m_imFlux.data_given())
					VecAppend(fluxDensity, m_imFlux[ip]);

				(*err_est_data)(side_list[side],sip) += scale * VecDot(fluxDensity, normal);
			}

			passedIPs += err_est_data->num_side_ips(side_list[side]);
		}
		UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
				<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
	}

// VOLUME TERMS //

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFE::compute_err_est_elem'");

	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
		{
			number total = 0.0;

		// diffusion //	TODO ONLY FOR (PIECEWISE) CONSTANT DIFFUSION TENSOR SO FAR!
		// div(D*grad(c))
		// nothing to do, as c is piecewise linear and div(D*grad(c)) disappears
		// if D is diagonal and c bilinear, this should also vanish (confirm this!)

		// convection // TODO ONLY FOR (PIECEWISE) CONSTANT OR DIVERGENCE-FREE
					  //      VELOCITY FIELDS SO FAR!
		// div(v*c) = div(v)*c + v*grad(c) -- gradC has been calculated above
			if (m_imVelocity.data_given())
				total += VecDot(m_imVelocity[ip], gradC);

		// general flux // TODO ONLY FOR DIVERGENCE-FREE FLUX FIELD SO FAR!
		// nothing to do

		// reaction //
			if (m_imReactionRate.data_given())
			{
				number val = 0.0;
				for (size_t sh = 0; sh < geo.num_sh(); sh++)
					val += u(_C_,sh) * m_shapeValues.shapeAtElemIP(sh,ip);

				total += m_imReactionRate[ip] * val;
			}

			if (m_imReaction.data_given())
			{
				total += m_imReaction[ip];
			}

			(*err_est_data)(elem_list[0],ip) += scale * total;
		}
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
}

//	computes the error estimator contribution (mass part) for one element
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
// note: mass parts only enter volume term

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFE::compute_err_est_elem'");

//	request geometry
	static const TFEGeom& geo = GeomProvider<TFEGeom>::get();

// 	loop integration points
	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
		{
			number total = 0.0;

		// mass scale //
			if (m_imMassScale.data_given())
			{
				number val = 0.0;
				for (size_t sh = 0; sh < geo.num_sh(); sh++)
					val += u(_C_,sh) * m_shapeValues.shapeAtElemIP(sh,ip);

				total += m_imMassScale[ip] * val;
			}

		// mass //
			if (m_imMass.data_given())
			{
				total += m_imMass[ip];
			}

			(*err_est_data)(elem_list[0],ip) += scale * total;
		}
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
}

//	computes the error estimator contribution (rhs part) for one element
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

// SIDE TERMS //

//	get the sides of the element
	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::side_type>::secure_container side_list;
	pErrEstGrid->associated_elements_sorted(side_list, (TElem*) elem);
	if (side_list.size() != (size_t) ref_elem_type::numSides)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFE::compute_err_est_elem'");

// loop sides
	size_t passedIPs = 0;
	for (size_t side = 0; side < (size_t) ref_elem_type::numSides; side++)
	{
		// normal on side
		MathVector<dim> normal;
		SideNormal<ref_elem_type,dim>(normal, side, vCornerCoords);
		VecNormalize(normal, normal);

		try
		{
			for (size_t sip = 0; sip < err_est_data->num_side_ips(side_list[side]); sip++)
			{
				size_t ip = passedIPs + sip;

			// vector source //
				if (m_imVectorSource.data_given())
					(*err_est_data)(side_list[side],sip) += scale * VecDot(m_imVectorSource[ip], normal);
			}

			passedIPs += err_est_data->num_side_ips(side_list[side]);
		}
		UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
				<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
	}

// VOLUME TERMS //

	if (!m_imSource.data_given()) return;

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFE::compute_err_est_elem'");

// source //
	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
			(*err_est_data)(elem_list[0],ip) += scale * m_imSource[ip];
	}
	UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
			<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
}

//	postprocesses the loop over all elements of one type in the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
fsh_err_est_elem_loop()
{
//	finish the element loop in the same way as the actual discretization
	this->template fsh_elem_loop<TElem, TFEGeom> ();
};

//    error estimation (end)     ///
// /////////////////////////////////


//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
lin_def_velocity(const LocalVector& u,
                     std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                     const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	// 	get current u and grad_u
		number shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
			shape_u += u(_C_,j) * geo.shape(ip, j);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	add to local defect
			VecScale(vvvLinDef[ip][_C_][i], geo.global_grad(ip, i),
			         	 	 	 	 	 	(-1)* geo.weight(ip) * shape_u);
		}
	}
}

//	computes the linearized defect w.r.t to the flux
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
lin_def_flux(const LocalVector& u,
             std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
             const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	add to local defect
			VecScale(vvvLinDef[ip][_C_][i], geo.global_grad(ip, i),
			         	 	 	 	 	 	(-1)* geo.weight(ip));
		}
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
lin_def_diffusion(const LocalVector& u,
                      std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
                      const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	MathVector<dim> grad_u;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	// 	get current u and grad_u
		VecSet(grad_u, 0.0);
		for(size_t j = 0; j < geo.num_sh(); ++j)
			VecScaleAppend(grad_u, u(_C_,j), geo.global_grad(ip, j));

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			for(size_t k = 0; k < (size_t)dim; ++k)
				for(size_t j = 0; j < (size_t)dim; ++j)
					(vvvLinDef[ip][_C_][i])(k,j) = grad_u[j] * geo.global_grad(ip, i)[k]
												* geo.weight(ip);
		}
	}
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
lin_def_reaction(const LocalVector& u,
                     std::vector<std::vector<number> > vvvLinDef[],
                     const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			const number val = geo.shape(ip, i) * geo.weight(ip);

		//	add to local defect
			vvvLinDef[ip][_C_][i] = val;
		}
	}
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
lin_def_reaction_rate(const LocalVector& u,
                         std::vector<std::vector<number> > vvvLinDef[],
                         const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	compute value of current solution at ip
		number shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
			shape_u += u(_C_,j) * geo.shape(ip, j);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			const number val = shape_u * geo.shape(ip, i) * geo.weight(ip);

		//	add to local defect
			vvvLinDef[ip][_C_][i] = val;
		}
	}
}


//	computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
lin_def_source(const LocalVector& u,
                   std::vector<std::vector<number> > vvvLinDef[],
                   const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	add contribution to local defect
			vvvLinDef[ip][_C_][i] = geo.shape(ip, i) * geo.weight(ip);
		}
	}
}

//	computes the linearized defect w.r.t to the "vector source"
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
lin_def_vector_source(const LocalVector& u,
                           std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                           const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	add to local defect
			VecScale(vvvLinDef[ip][_C_][i], geo.global_grad(ip, i),
			         	 	 	 	 	 	 geo.weight(ip));
		}
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
lin_def_mass_scale(const LocalVector& u,
                       std::vector<std::vector<number> > vvvLinDef[],
                       const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	compute value of current solution at ip
		number shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
			shape_u += u(_C_,j) * geo.shape(ip, j);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			const number val = shape_u * geo.shape(ip, i) * geo.weight(ip);

		//	add to local defect
			vvvLinDef[ip][_C_][i] = val;
		}
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
lin_def_mass(const LocalVector& u,
                std::vector<std::vector<number> > vvvLinDef[],
                const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			const number val = geo.shape(ip, i) * geo.weight(ip);

		//	add to local defect
			vvvLinDef[ip][_C_][i] = val;
		}
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
ex_value(number vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         const LocalVector& u,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<TFEGeom::dim> vLocIP[],
         const size_t nip,
         bool bDeriv,
         std::vector<std::vector<number> > vvvDeriv[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	reference dimension
	static const int refDim = reference_element_traits<TElem>::dim;

//	reference object id
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	FE ip
	if(vLocIP == geo.local_ips())
	{
	//	Loop ips
		for(size_t ip = 0; ip < geo.num_ip(); ++ip)
		{
		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t sh = 0; sh < geo.num_sh(); ++sh)
				vValue[ip] += u(_C_, sh) * geo.shape(ip, sh);

		//	compute derivative w.r.t. to unknowns iff needed
			if(bDeriv)
				for(size_t sh = 0; sh < geo.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = geo.shape(ip, sh);
		}
	}
// 	general case
	else
	{
	//	request for trial space
		try{
		const LocalShapeFunctionSet<refDim>& rTrialSpace
			 = LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

	//	number of shape functions
		const size_t numSH = rTrialSpace.num_sh();

	//	storage for shape function at ip
		std::vector<number> vShape(numSH);

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
		UG_CATCH_THROW("ConvectionDiffusion::ex_value: trial space missing.");
	}
}

template<typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
ex_grad(MathVector<dim> vValue[],
        const MathVector<dim> vGlobIP[],
        number time, int si,
        const LocalVector& u,
        GridObject* elem,
        const MathVector<dim> vCornerCoords[],
        const MathVector<TFEGeom::dim> vLocIP[],
        const size_t nip,
        bool bDeriv,
        std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	reference dimension
	static const int refDim = reference_element_traits<TElem>::dim;

//	reference object id
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	FE
	if(vLocIP == geo.local_ips())
	{
	//	Loop ip
		for(size_t ip = 0; ip < geo.num_ip(); ++ip)
		{
			VecSet(vValue[ip], 0.0);
			for(size_t sh = 0; sh < geo.num_sh(); ++sh)
				VecScaleAppend(vValue[ip], u(_C_, sh), geo.global_grad(ip, sh));

			if(bDeriv)
				for(size_t sh = 0; sh < geo.num_sh(); ++sh)
					vvvDeriv[ip][_C_][sh] = geo.global_grad(ip, sh);
		}
	}
// 	general case
	else
	{
	//	request for trial space
		try{
		const LocalShapeFunctionSet<refDim>& rTrialSpace
			 = LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);

	//	number of shape functions
		const size_t numSH = rTrialSpace.num_sh();

	//	storage for shape function at ip
		std::vector<MathVector<refDim> > vLocGrad(numSH);
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
		UG_CATCH_THROW("ConvectionDiffusion::ex_grad: trial space missing.");
	}
};

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void ConvectionDiffusionFE<Domain1d>::
register_all_funcs(const LFEID& lfeid, const int quadOrder)
{
//	RegularEdge
	register_func<RegularEdge, DimFEGeometry<dim> >();
}
#endif

#ifdef UG_DIM_2
template<>
void ConvectionDiffusionFE<Domain2d>::
register_all_funcs(const LFEID& lfeid, const int quadOrder)
{
	register_func<RegularEdge, DimFEGeometry<dim, 1> >();

	const int order = lfeid.order();
	if(quadOrder != 2*order+1 || lfeid.type() != LFEID::LAGRANGE)
	{
		register_func<Triangle, DimFEGeometry<dim> >();
		register_func<Quadrilateral, DimFEGeometry<dim> >();
		return;
	}

//	special compiled cases

//	Triangle
	switch(order)
	{
		case 1:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 1>, GaussQuadrature<ReferenceTriangle, 3> > FEGeom;
				 register_func<Triangle, FEGeom >(); break;}
		case 2:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 2>, GaussQuadrature<ReferenceTriangle, 5> > FEGeom;
				 register_func<Triangle, FEGeom >(); break;}
		case 3:	{typedef FEGeometry<Triangle, dim, LagrangeLSFS<ReferenceTriangle, 3>, GaussQuadrature<ReferenceTriangle, 7> > FEGeom;
				 register_func<Triangle, FEGeom >(); break;}
		default: register_func<Triangle, DimFEGeometry<dim> >();  break;
	}

//	Quadrilateral
	switch(order) {
		case 1:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 1>, GaussQuadrature<ReferenceQuadrilateral, 3> > FEGeom;
				 register_func<Quadrilateral, FEGeom >(); break;}
		case 2:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 2>, GaussQuadrature<ReferenceQuadrilateral, 7> > FEGeom;
				 register_func<Quadrilateral, FEGeom >(); break;}
		case 3:	{typedef FEGeometry<Quadrilateral, dim, LagrangeLSFS<ReferenceQuadrilateral, 3>, GaussQuadrature<ReferenceQuadrilateral, 11> > FEGeom;
				 register_func<Quadrilateral, FEGeom >(); break;}
		default: register_func<Quadrilateral, DimFEGeometry<dim> >();  break;
	}
}
#endif

#ifdef UG_DIM_3
template<>
void ConvectionDiffusionFE<Domain3d>::
register_all_funcs(const LFEID& lfeid, const int quadOrder)
{
	register_func<RegularEdge, DimFEGeometry<dim, 1> >();
	register_func<Triangle, DimFEGeometry<dim, 2> >();
	register_func<Quadrilateral, DimFEGeometry<dim, 2> >();

	const int order = lfeid.order();
	if(quadOrder != 2*order+1 || lfeid.type() != LFEID::LAGRANGE)
	{
		register_func<Tetrahedron, DimFEGeometry<dim> >();
		register_func<Prism, DimFEGeometry<dim> >();
		register_func<Pyramid, DimFEGeometry<dim> >();
		register_func<Hexahedron, DimFEGeometry<dim> >();
		register_func<Octahedron, DimFEGeometry<dim> >();
		return;
	}

//	special compiled cases

//	Tetrahedron
	switch(order)
	{
		case 1:	{typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<ReferenceTetrahedron, 1>, GaussQuadrature<ReferenceTetrahedron, 3> > FEGeom;
				 register_func<Tetrahedron, FEGeom >(); break;}
		case 2:	{typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<ReferenceTetrahedron, 2>, GaussQuadrature<ReferenceTetrahedron, 5> > FEGeom;
				 register_func<Tetrahedron, FEGeom >(); break;}
		case 3:	{typedef FEGeometry<Tetrahedron, dim, LagrangeLSFS<ReferenceTetrahedron, 3>, GaussQuadrature<ReferenceTetrahedron, 7> > FEGeom;
				 register_func<Tetrahedron, FEGeom >(); break;}
		default: register_func<Tetrahedron, DimFEGeometry<dim> >();  break;
	}

//	Prism
	switch(order) {
		default: register_func<Prism, DimFEGeometry<dim> >();  break;
	}

//	Pyramid
	switch(order)
	{
		default: register_func<Pyramid, DimFEGeometry<dim> >();  break;
	}

//	Octahedron
	switch(order)
	{
		default: register_func<Octahedron, DimFEGeometry<dim> >();  break;
	}

//	Hexahedron
	switch(order)
	{
		case 1:	{typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<ReferenceHexahedron, 1>, GaussQuadrature<ReferenceHexahedron, 3> > FEGeom;
				 register_func<Hexahedron, FEGeom >(); break;}
		case 2:	{typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<ReferenceHexahedron, 2>, GaussQuadrature<ReferenceHexahedron, 7> > FEGeom;
				 register_func<Hexahedron, FEGeom >(); break;}
		case 3:	{typedef FEGeometry<Hexahedron, dim, LagrangeLSFS<ReferenceHexahedron, 3>, GaussQuadrature<ReferenceHexahedron, 11> > FEGeom;
				 register_func<Hexahedron, FEGeom >(); break;}
		default: register_func<Hexahedron, DimFEGeometry<dim> >();  break;
	}
}
#endif

template <typename TDomain>
template <typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	static const int refDim = reference_element_traits<TElem>::dim;

	this->clear_add_fct(id);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFEGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TFEGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template fsh_elem_loop<TElem, TFEGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFEGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFEGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFEGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFEGeom>);
	this->set_add_rhs_elem_fct(  id, &T::template add_rhs_elem<TElem, TFEGeom>);

// error estimator parts
	this->set_prep_err_est_elem_loop(id, &T::template prep_err_est_elem_loop<TElem, TFEGeom>);
	this->set_prep_err_est_elem(id, &T::template prep_err_est_elem<TElem, TFEGeom>);
	this->set_compute_err_est_A_elem(id, &T::template compute_err_est_A_elem<TElem, TFEGeom>);
	this->set_compute_err_est_M_elem(id, &T::template compute_err_est_M_elem<TElem, TFEGeom>);
	this->set_compute_err_est_rhs_elem(id, &T::template compute_err_est_rhs_elem<TElem, TFEGeom>);
	this->set_fsh_err_est_elem_loop(id, &T::template fsh_err_est_elem_loop<TElem, TFEGeom>);

//	set computation of linearized defect w.r.t velocity
	m_imDiffusion.		set_fct(id, this, &T::template lin_def_diffusion<TElem, TFEGeom>);
	m_imVelocity. 		set_fct(id, this, &T::template lin_def_velocity<TElem, TFEGeom>);
	m_imFlux.	 		set_fct(id, this, &T::template lin_def_flux<TElem, TFEGeom>);
	m_imReactionRate. 	set_fct(id, this, &T::template lin_def_reaction_rate<TElem, TFEGeom>);
	m_imReaction. 		set_fct(id, this, &T::template lin_def_reaction<TElem, TFEGeom>);
	m_imSource.	  		set_fct(id, this, &T::template lin_def_source<TElem, TFEGeom>);
	m_imVectorSource.	set_fct(id, this, &T::template lin_def_vector_source<TElem, TFEGeom>);
	m_imMassScale.		set_fct(id, this, &T::template lin_def_mass_scale<TElem, TFEGeom>);
	m_imMass.	  		set_fct(id, this, &T::template lin_def_mass<TElem, TFEGeom>);

//	exports
	m_exValue->	template set_fct<T,refDim>(id, this, &T::template ex_value<TElem, TFEGeom>);
	m_exGrad->	template set_fct<T,refDim>(id, this, &T::template ex_grad<TElem, TFEGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ConvectionDiffusionFE<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConvectionDiffusionFE<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionFE<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug

