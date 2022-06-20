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

#include "lib_disc/spatial_disc/disc_util/fe_geom.h"
#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"
#include "lib_disc/local_finite_element/lagrange/lagrangep1.h"
#include "lib_disc/local_finite_element/raviart-thomas/raviart_thomas.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"

#include "common/math/math_vector_matrix/math_matrix_vector_functions.h"
#include "common/math/math_vector_matrix/math_matrix_functions.h"
#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "weak_form_rt.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

static const std::string WeakFormulationID = "L2Projection:";

template<typename TDomain>
WeakFormulationFE<TDomain>::
WeakFormulationFE(const char* functions, const char* subsets)
 : ConvectionDiffusionBase<TDomain>(functions,subsets),
	m_bQuadOrderUserDef(false)
{
	this->clear_add_fct();
}

template<typename TDomain>
void WeakFormulationFE<TDomain>::set_quad_order(size_t order)
{
	m_quadOrder = order;
	m_bQuadOrderUserDef = true;
}

template<typename TDomain>
void WeakFormulationFE<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	//	check number of fcts
	UG_COND_THROW(vLfeID.size() != 1, WeakFormulationID <<
			" Wrong number of functions given. Need exactly "<<1);

	//	check that not ADAPTIVE
	UG_COND_THROW(vLfeID[0].order() < 0, WeakFormulationID <<
			" Adaptive order not implemented.");

	//	set order
	m_lfeID = vLfeID[0];
	if(!m_bQuadOrderUserDef) m_quadOrder = 2*m_lfeID.order()+1;

	register_all_funcs(m_lfeID, m_quadOrder);
}

template<typename TDomain>
bool WeakFormulationFE<TDomain>::
use_hanging() const
{
	return false;
}

////////////////////////////////////////////////////////////////////////////////
// Assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	// Check for explicit terms.
	//UG_COND_THROW(m_imSourceExpl.data_given() || m_imReactionExpl.data_given() || m_imReactionRateExpl.data_given(),
	//		"WeakFormulationFE: Explicit terms not implemented.");

	// Request geometry.
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	//	Prepare geometry for type and order.
	try{
		geo.update_local(roid, m_lfeID, m_quadOrder);
	} UG_CATCH_THROW(WeakFormulationID << ":prep_elem_loop: Cannot update Finite Element Geometry:" << m_lfeID);

	//	Set local positions.
	static const int refDim = TElem::dim;
	m_imDiffusion.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	/*m_imVelocity.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imFlux.template  		set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);*/
	m_imSource.template    set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imVectorSource.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
/*	m_imReactionRate.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imReaction.template  set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imMassScale.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);
	m_imMass.template 	   set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);*/
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
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
	// m_imVelocity. 	set_global_ips(geo.global_ips(), geo.num_ip());
	// m_imFlux. 		set_global_ips(geo.global_ips(), geo.num_ip());*/
	m_imSource.   	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imVectorSource.set_global_ips(geo.global_ips(), geo.num_ip());
	/*m_imReactionRate.set_global_ips(geo.global_ips(), geo.num_ip());
	m_imReaction. 	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imMassScale.	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imMass.	  	set_global_ips(geo.global_ips(), geo.num_ip());*/
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
	std::cerr << "NumIP=" << geo.num_ip()<< std::endl;
	MathVector<dim> DinvShapej;
	MathMatrix<dim,dim> Dinv;
//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop trial space
		for(size_t j = 0; j < geo.num_sh(); ++j)
		{
		//	Diffusion
			if(m_imDiffusion.data_given()) Inverse(Dinv, m_imDiffusion[ip]);
			else Dinv = 1.0;

			MatVecMult(DinvShapej, Dinv, geo.shape(ip,j));
		//  Convection
		//	if(m_imVelocity.data_given())
		//		VecScaleAppend(Dgrad, -1*geo.shape(ip,j), m_imVelocity[ip]);

		//	no explicit dependency on m_imFlux

		//	loop test space
			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
			//	compute integrand
				number integrand = VecDot(geo.shape(ip,i), DinvShapej);
				std::cerr << i << "," << j <<":"<< integrand << "*" << geo.weight(ip);
			// 	Reaction
				//if(m_imReactionRate.data_given())
				//	integrand += m_imReactionRate[ip] * geo.shape(ip, j) * geo.shape(ip, i);

			//	no explicit dependency on m_imReaction

			//	multiply by weight
				integrand *= geo.weight(ip);
				std::cerr << "="<< integrand << std::endl;
			//	add to local matrix
				J(_U_, i, _U_, j) += integrand*geo.signum(j)*geo.signum(i);
			}
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
}


template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	number integrand;
	MathVector<dim> ShapeU, DinvShapeU;
	MathMatrix<dim,dim> Dinv;
	MathMatrix<dim,dim> GradU;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	// 	get current u and grad_u
		VecSet(ShapeU, 0.0);
		DinvShapeU=0.0;


		for(size_t j = 0; j < geo.num_sh(); ++j)
		{ VecScaleAppend(ShapeU, u(_U_,j)*geo.signum(j), geo.shape(ip, j)); }

	// 	Diffusion
		if(m_imDiffusion.data_given()) Inverse(Dinv, m_imDiffusion[ip]);
		else Dinv = 1.0;
		MatVecMult(DinvShapeU, Dinv, ShapeU);

	// 	Convection
	//if(m_imVelocity.data_given())
	//		VecScaleAppend(Dgrad_u, -1*shape_u, m_imVelocity[ip]);

	// 	Convection
	//	if(m_imFlux.data_given())
	//		VecScaleAppend(Dgrad_u, 1.0, m_imFlux[ip]);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute integrand
			integrand = VecDot(DinvShapeU, geo.shape(ip, i))*geo.signum(i);

		// 	add Reaction Rate
		//	if(m_imReactionRate.data_given())
		//		integrand += m_imReactionRate[ip] * shape_u * geo.shape(ip, i);

		// 	add Reaction
		//	if(m_imReaction.data_given())
		//		integrand += m_imReaction[ip] * geo.shape(ip, i);


		//	add to local defect
			d(_U_, i) += integrand*geo.weight(ip);
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
};

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	skip if no source present
	if(!m_imSource.data_given() && !m_imVectorSource.data_given()) return;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	Loop test spaces

		// only do this if (volume) source is given
		if(m_imSource.data_given())
		{
			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
			//	add contribution to local defect
				const double divNi = geo.global_div(ip,i)*geo.signum(i);
				d(_U_, i) += geo.weight(ip) * m_imSource[ip] * divNi;
			}
		}

		//	only do this if vector source is given
		if(m_imVectorSource.data_given())
		{
			for(size_t i = 0; i < geo.num_sh(); ++i)
			{
			//	add contribution to local defect
				d(_U_, i) += geo.weight(ip) * geo.signum(i) * VecDot(m_imVectorSource[ip], geo.shape(ip, i));
			}
		}
	}
}

// ////////////////////////////////
//   error estimation (begin)   ///

//	prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
prep_err_est_elem_loop(const ReferenceObjectID roid, const int si)
{

}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

}

//	computes the error estimator contribution (stiffness part) for one element
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{

}

//	computes the error estimator contribution (mass part) for one element
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
}

//	computes the error estimator contribution (rhs part) for one element
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
}

//	postprocesses the loop over all elements of one type in the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
fsh_err_est_elem_loop()
{
//	finish the element loop in the same way as the actual discretization
	this->template fsh_elem_loop<TElem, TFEGeom> ();
};

//    error estimation (end)     ///
// /////////////////////////////////


//	computes the linearized defect w.r.t to the diffusion tensor
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
lin_def_diffusion(const LocalVector& u,
                      std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
                      const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	MathVector<dim> grad_u;
	/*
//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	// 	get current u and grad_u
		VecSet(grad_u, 0.0);
		for(size_t j = 0; j < geo.num_sh(); ++j)
			VecScaleAppend(grad_u, u(_U_,j), geo.global_grad(ip, j));

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			for(size_t k = 0; k < (size_t)dim; ++k)
				for(size_t j = 0; j < (size_t)dim; ++j)
					(vvvLinDef[ip][_U_][i])(k,j) = grad_u[j] * geo.global_grad(ip, i)[k]
												* geo.weight(ip);
		}
	}
	*/
}



//!	computes the linearized defect w.r.t to the source

template<typename TDomain>
template <typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
lin_def_source(const LocalVector& u,
                   std::vector<std::vector<number> > vvvLinDef[],
                   const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	add contribution to local defect
		// where d(_U_, i) += geo.weight(ip) * m_imSource[ip] * divNi;
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			//const MathMatrix<dim,dim>& GradNi=geo.global_grad(ip, i);
			//const double divNi = Trace(GradNi)*geo.signum(i);
			const double divNi = geo.global_div(ip, i)*geo.signum(i);
			vvvLinDef[ip][_U_][i] = geo.weight(ip) * divNi;
		}
	}

}


//!	computes the linearized defect w.r.t to the "vector source"
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
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
		// d(_U_, i) += geo.weight(ip) * VecDot(m_imVectorSource[ip], geo.shape(ip, i));
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
			VecScale(vvvLinDef[ip][_U_][i], geo.shape(ip, i), geo.signum(i)*geo.weight(ip));
		}
	}

}


/*
//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
lin_def_mass_scale(const LocalVector& u,
                       std::vector<std::vector<number> > vvvLinDef[],
                       const size_t nip)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
/*
//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	compute value of current solution at ip
		number shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
			shape_u += u(_U_,j) * geo.shape(ip, j);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			const number val = shape_u * geo.shape(ip, i) * geo.weight(ip);

		//	add to local defect
			vvvLinDef[ip][_U_][i] = val;
		}
	}
}
*/

/*
//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
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
			vvvLinDef[ip][_U_][i] = val;
		}
	}

}
*/
//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
ex_divergence(number  vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         const LocalVector& u,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<TFEGeom::dim> vLocIP[],
         const size_t nip,
         bool bDeriv,
         std::vector<std::vector<number > > vvvDeriv[])
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
			for(size_t j = 0; j < geo.num_sh(); ++j)
			{
				const double divNj = geo.global_div(ip, j)*geo.signum(j);
				vValue[ip] += u(_U_, j) * divNj ;

				//	compute derivative w.r.t. to unknowns iff needed
				if(bDeriv) vvvDeriv[ip][_U_][j] = divNj;
			}
		}
	}
// 	general case
	else
	{
		UG_COND_THROW( m_lfeID.order()!=1,
				"Huhh: Implementation assumes constant gradient!, but order =" << m_lfeID.order() << std::endl);

		//std::cout << "Divergence1:" << roid << ", "<< m_lfeID << ", "<< nip <<  std::endl;
		//std::cout << "Divergence1:" << vLocIP[0] << m_lfeID.order() << std::endl;
	//	request for trial space
		try{

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{

		//	compute concentration at ip
			vValue[ip] = 0.0;
			for(size_t j = 0; j < geo.num_sh(); ++j)
			{
				//const MathMatrix<dim,dim>& GradNj=geo.global_grad(ip, j);
				//const double divNj = Trace(GradNj)*geo.signum(j);
				const double divNj = geo.global_div(ip, j)*geo.signum(j);
				vValue[ip] += u(_U_, j) * divNj;

				//	compute derivative w.r.t. to unknowns iff needed
				if(bDeriv) vvvDeriv[ip][_U_][j] = divNj;
			}
		}

		}
		UG_CATCH_THROW("WeakFormulation::ex_divergence: trial space missing.");
	}

}

template<typename TDomain>
template <typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::
ex_field(MathVector<dim> vValue[],
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
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	typedef RaviartThomasLSFS<ref_elem_type> RT1;
	typedef typename RT1::shape_type shape_type;


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
			for(size_t j = 0; j < geo.num_sh(); ++j)
				VecScaleAppend(vValue[ip], u(_U_,j), geo.shape(ip, j));

			if(bDeriv)
				for(size_t j = 0; j < geo.num_sh(); ++j)
					vvvDeriv[ip][_U_][j] = geo.shape(ip, j);
		}
	}
	else
	{


	try{
		//	request for trial space
			const RT1& rTrialSpace = Provider<RT1>::get();

		//	number of shape functions
			const size_t numSH = rTrialSpace.num_sh();

		//	loop ips
			for(size_t ip = 0; ip < nip; ++ip)
			{
				VecSet(vValue[ip], 0.0);

			//	compute value at ip
				for(size_t j = 0; j < numSH; ++j)
				{
					shape_type myShape = rTrialSpace.shape(j, vLocIP[ip]);
					VecScaleAppend(vValue[ip], u(_U_, j), myShape);

					//	compute derivative w.r.t. to unknowns iff needed
					if(bDeriv) vvvDeriv[ip][_U_][j] = myShape;
				}
			}

		} UG_CATCH_THROW("WeakFormulation::ex_field: trial space missing.");
	}
};

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void WeakFormulationFE<Domain1d>::
register_all_funcs(const LFEID& lfeid, const int quadOrder)
{
//	RegularEdge
	//register_func<RegularEdge, DimFEGeometry<dim> >();
}
#endif

#ifdef UG_DIM_2
template<>
void WeakFormulationFE<Domain2d>::
register_all_funcs(const LFEID& lfeid, const int quadOrder)
{

	int order = 1;
//	Triangle
	switch(order)
	{
		case 1:
		{
			typedef FEGeometryMixed<Triangle, dim, RaviartThomasLSFS<ReferenceTriangle>, GaussQuadrature<ReferenceTriangle, 3> > FEGeomMixed;
			register_func<Triangle, FEGeomMixed>();
		}

		{
			typedef FEGeometryMixed<Quadrilateral, dim, RaviartThomasLSFS<ReferenceQuadrilateral>, GaussQuadrature<ReferenceQuadrilateral, 3> > FEGeomMixed;
			register_func<Quadrilateral, FEGeomMixed>();
		}

		break;

		/*case 2:	{typedef FEGeometry<Triangle, dim, RaviartThomasLSFS<ReferenceTriangle>, GaussQuadrature<ReferenceTriangle, 5> > FEGeomMixed;
				 register_func<Triangle, FEGeomMixed >(); break;}
		case 3:	{typedef FEGeometry<Triangle, dim, RaviartThomasLSFS<ReferenceTriangle>, GaussQuadrature<ReferenceTriangle, 7> > FEGeomMixed;
				 register_func<Triangle, FEGeomMixed>(); break;}*/
		//default: register_func<Triangle, DimFEGeometry<dim> >();  break;
	}

//	Quadrilateral
	switch(order) {

		//default: register_func<Quadrilateral, DimFEGeometry<dim> >();  break;
	}
}
#endif

#ifdef UG_DIM_3
template<>
void WeakFormulationFE<Domain3d>::
register_all_funcs(const LFEID& lfeid, const int quadOrder)
{

}
#endif

template <typename TDomain>
template <typename TElem, typename TFEGeom>
void WeakFormulationFE<TDomain>::register_func()
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
	/*m_imVelocity. 		set_fct(id, this, &T::template lin_def_velocity<TElem, TFEGeom>);
	m_imFlux.	 		set_fct(id, this, &T::template lin_def_flux<TElem, TFEGeom>);
	m_imReactionRate. 	set_fct(id, this, &T::template lin_def_reaction_rate<TElem, TFEGeom>);
	m_imReaction. 		set_fct(id, this, &T::template lin_def_reaction<TElem, TFEGeom>);*/
	m_imSource.	  		set_fct(id, this, &T::template lin_def_source<TElem, TFEGeom>);
	m_imVectorSource.	set_fct(id, this, &T::template lin_def_vector_source<TElem, TFEGeom>);
/*	m_imMassScale.		set_fct(id, this, &T::template lin_def_mass_scale<TElem, TFEGeom>);
	m_imMass.	  		set_fct(id, this, &T::template lin_def_mass<TElem, TFEGeom>);
*/
//	exports
	m_exValue->	template set_fct<T,refDim>(id, this, &T::template ex_divergence<TElem, TFEGeom>);
	m_exGrad->	template set_fct<T,refDim>(id, this, &T::template ex_field<TElem, TFEGeom>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class WeakFormulationFE<Domain1d>;
#endif
#ifdef UG_DIM_2
template class WeakFormulationFE<Domain2d>;
#endif
#ifdef UG_DIM_3
template class WeakFormulationFE<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug

