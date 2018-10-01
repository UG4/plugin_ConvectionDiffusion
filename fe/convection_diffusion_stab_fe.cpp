/*
 * Copyright (c) 2018-:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
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

#include "convection_diffusion_stab_fe.h"

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
ConvectionDiffusionStabFE<TDomain>::
ConvectionDiffusionStabFE(const char* functions, const char* subsets)
 : base_type(functions,subsets),
	m_bQuadOrderUserDef(false), m_stabParamM(0.0), m_stabParamA(0.0)
{
	this->clear_add_fct();
}


template<typename TDomain>
ConvectionDiffusionStabFE<TDomain>::
ConvectionDiffusionStabFE(const char* functions, const char* subsets, double stabM)
 : base_type(functions,subsets),
	m_bQuadOrderUserDef(false), m_stabParamM(stabM), m_stabParamA(0.0)
{
	this->clear_add_fct();
}

template<typename TDomain>
ConvectionDiffusionStabFE<TDomain>::
ConvectionDiffusionStabFE(const char* functions, const char* subsets, double stabM, double stabA)
 : base_type(functions,subsets),
	m_bQuadOrderUserDef(false), m_stabParamM(stabM), m_stabParamA(stabA)
{
	this->clear_add_fct();
}


template<typename TDomain>
void ConvectionDiffusionStabFE<TDomain>::
set_quad_order(size_t order)
{
	m_quadOrder = order;
	m_bQuadOrderUserDef = true;
}

template<typename TDomain>
void ConvectionDiffusionStabFE<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid)
{
	//	check number of fcts
	if(vLfeID.size() != 1)
		UG_THROW("ConvectionDiffusionStab: Wrong number of functions given. "
				"Need exactly "<<1);

	//	check that not ADAPTIVE
	if(vLfeID[0].order() < 1)
		UG_THROW("ConvectionDiffusionStab: Adaptive order not implemented.");

	//	set order
	m_lfeID = vLfeID[0];
	if(!m_bQuadOrderUserDef) m_quadOrder = 2*m_lfeID.order()+1;

	register_all_funcs(m_lfeID, m_quadOrder);
}

template<typename TDomain>
bool ConvectionDiffusionStabFE<TDomain>::
use_hanging() const
{
	return false;
}

////////////////////////////////////////////////////////////////////////////////
// Assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionStabFE<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{

//	request geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	prepare geometry for type and order
	try{
		geo.update_local(roid, m_lfeID, m_quadOrder);
	}UG_CATCH_THROW("ConvectionDiffusion::prep_elem_loop:"
					" Cannot update Finite Element Geometry.");

//	set local positions
//	static const int refDim = TElem::dim;
//	m_imDiffusion.template set_local_ips<refDim>(geo.local_ips(), geo.num_ip(), false);

}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionStabFE<TDomain>::
fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionStabFE<TDomain>::
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
	// m_imDiffusion.	set_global_ips(geo.global_ips(), geo.num_ip());

}


template<int dim>
void PointsBoundingBox(size_t npoints, const MathVector<dim> points[], MathVector<dim> &vMinBB, MathVector<dim> &vMaxBB)
{
	// determine bounding box
	vMinBB= points[0];
	vMaxBB = points[0];

	for(size_t ii = 1; ii < npoints; ++ii)
	{
		for(int i = 0; i < dim; ++i)
		{
			const MathVector<dim>& v = points[ii];
			if(v[i] < vMinBB[i]) vMinBB[i] = v[i];
			else if(v[i] > vMaxBB[i]) vMaxBB[i] = v[i];
		}
	}
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionStabFE<TDomain>::
add_jac_X_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], double m_stabParam)
{
	const TFEGeom& geo = //	request geometry
				GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
	/* const number scale = // Scaling factor
				m_stabParam *ElementDiameterSq<GridObject, TDomain>(*elem, *this->domain()); */


	MathVector<dim> vBoundingBox[2];
	PointsBoundingBox<dim>(NumVertices((TElem*)elem), vCornerCoords, vBoundingBox[0], vBoundingBox[1]);

	// rescale using bounding box
	for (size_t ip = 0; ip < geo.num_ip(); ++ip){
		for (size_t psh = 0; psh < geo.num_sh(); ++psh){
			for (size_t psh2 = 0; psh2 < geo.num_sh(); ++psh2){

				MathVector<dim> grad = geo.global_grad(ip, psh); // *scale
				for(size_t i = 0; i < dim; ++i)
				{
					double hi = vBoundingBox[1][i] - vBoundingBox[0][i];
					grad[i] *= (hi*hi);
				}

				J(_C_, psh, _C_, psh2) += m_stabParam* VecDot(grad, geo.global_grad(ip, psh2))
											* geo.weight(ip);
			}
		}
	}

}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionStabFE<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	if (m_stabParamM != 0.0) add_jac_X_elem<TElem,TFEGeom>(J, u, elem,vCornerCoords, m_stabParamM);
}

template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionStabFE<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	if (m_stabParamA != 0.0) add_jac_X_elem<TElem,TFEGeom>(J, u, elem,vCornerCoords, m_stabParamA);
}
/*
template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionStabFE<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	const TFEGeom& geo =   // Geometry
			GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);
	const number scale =   // Scaling factor
			m_stabParam * ElementDiameterSq<GridObject, TDomain>(*elem, *this->domain());
	MathVector<dim> gradU; // Temporary

	//	Loop integration points.
	for (size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
		//	Loop trial space
		VecSet(gradU, 0.0);
		for (size_t psh = 0; psh < geo.num_sh(); ++psh)
			VecScaleAppend(gradU, u(_C_, psh), geo.global_grad(ip, psh));

		for (size_t psh = 0; psh < geo.num_sh(); ++psh)
		{
			d(_C_, psh) += scale * VecDot(geo.global_grad(ip, psh), gradU) * geo.weight(ip);
		}
	}

}
*/




////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void ConvectionDiffusionStabFE<Domain1d>::
register_all_funcs(const LFEID& lfeid, const int quadOrder)
{
//	RegularEdge
	register_func<RegularEdge, DimFEGeometry<dim> >();
}
#endif

#ifdef UG_DIM_2
template<>
void ConvectionDiffusionStabFE<Domain2d>::
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
void ConvectionDiffusionStabFE<Domain3d>::
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
void ConvectionDiffusionStabFE<TDomain>::register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
//	static const int refDim = reference_element_traits<TElem>::dim;

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
/*	m_imDiffusion.		set_fct(id, this, &T::template lin_def_diffusion<TElem, TFEGeom>);
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
	m_exGrad->	template set_fct<T,refDim>(id, this, &T::template ex_grad<TElem, TFEGeom>);*/
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ConvectionDiffusionStabFE<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConvectionDiffusionStabFE<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionStabFE<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug

