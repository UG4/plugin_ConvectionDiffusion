/*
 * convection_diffusion_fe.cpp
 *
 *  Created on: 02.08.2010
 *      Author: andreasvogel
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
	this->enable_fast_add_elem(true);
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
prep_elem(TElem* elem, const LocalVector& u)
{
//	get corners
	m_vCornerCoords = this->template element_corners<TElem>(elem);

//	request geometry
	TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	try{
		geo.update(elem, &m_vCornerCoords[0]);
	}
	UG_CATCH_THROW("ConvectionDiffusion::prep_elem:"
					" Cannot update Finite Element Geometry.");

//	set global positions for rhs
	m_imDiffusion.	set_global_ips(geo.global_ips(), geo.num_ip());
	m_imVelocity. 	set_global_ips(geo.global_ips(), geo.num_ip());
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
add_jac_A_elem(LocalMatrix& J, const LocalVector& u)
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
add_jac_M_elem(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	loop test space
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	loop trial space
			for(size_t j= 0; j < geo.num_sh(); ++j)
			{
			//	compute integrand
				number val = geo.shape(ip, i) *geo.shape(ip, j) * geo.weight(ip);

			//	add MassScale
				if(m_imMassScale.data_given())
					val *= m_imMassScale[ip];

			//	no explicit dependency on m_imMass

			//	add to local matrix
				J(_C_, i, _C_, j) += val;
			}
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TFEGeom>
void ConvectionDiffusionFE<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u)
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
add_def_M_elem(LocalVector& d, const LocalVector& u)
{
//	request geometry
	const TFEGeom& geo = GeomProvider<TFEGeom>::get(m_lfeID, m_quadOrder);

	number shape_u;

//	loop integration points
	for(size_t ip = 0; ip < geo.num_ip(); ++ip)
	{
	//	compute value of current solution at ip
		shape_u = 0.0;
		for(size_t j = 0; j < geo.num_sh(); ++j)
			shape_u += u(_C_,j) * geo.shape(ip, j);

	//	loop test spaces
		for(size_t i = 0; i < geo.num_sh(); ++i)
		{
		//	compute contribution
			number val = shape_u;

		//	add MassScaling
			if(m_imMassScale.data_given())
				val *= m_imMassScale[ip];

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
add_rhs_elem(LocalVector& d)
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
			for(size_t k=0; k < (size_t)dim; ++k)
				for(size_t j = 0; j < (size_t)dim; ++j)
					(vvvLinDef[ip][_C_][i])(k,j) += grad_u[j] * geo.global_grad(ip, i)[k]
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
         GeometricObject* elem,
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

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;

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
		const LocalShapeFunctionSet<dim>& rTrialSpace
			 = LocalShapeFunctionSetProvider::get<dim>(roid, m_lfeID);

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
        GeometricObject* elem,
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

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;

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
		const LocalShapeFunctionSet<dim>& rTrialSpace
			 = LocalShapeFunctionSetProvider::get<dim>(roid, m_lfeID);

	//	storage for shape function at ip
		MathVector<refDim> vLocGrad[numSH];
		MathVector<refDim> locGrad;

	//	Reference Mapping
		MathMatrix<dim, refDim> JTInv;
		ReferenceMapping<ref_elem_type, dim> mapping(&m_vCornerCoords[0]);

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
//	Edge
	register_func<Edge, DimFEGeometry<dim> >();
}
#endif

#ifdef UG_DIM_2
template<>
void ConvectionDiffusionFE<Domain2d>::
register_all_funcs(const LFEID& lfeid, const int quadOrder)
{
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
	const int order = lfeid.order();
	if(quadOrder != 2*order+1 || lfeid.type() != LFEID::LAGRANGE)
	{
		register_func<Tetrahedron, DimFEGeometry<dim> >();
		register_func<Prism, DimFEGeometry<dim> >();
		register_func<Pyramid, DimFEGeometry<dim> >();
		register_func<Hexahedron, DimFEGeometry<dim> >();
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
		case 1:	{typedef FEGeometry<Prism, dim, LagrangeLSFS<ReferencePrism, 1>, GaussQuadrature<ReferencePrism, 2> > FEGeom;
				 register_func<Prism, FEGeom >(); break;}
		default: register_func<Prism, DimFEGeometry<dim> >();  break;
	}

//	Pyramid
	switch(order)
	{
		default: register_func<Pyramid, DimFEGeometry<dim> >();  break;
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

	this->enable_fast_add_elem(true);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TFEGeom>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TFEGeom>);
	this->set_fsh_elem_loop_fct( id, &T::template fsh_elem_loop<TElem, TFEGeom>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TFEGeom>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TFEGeom>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TFEGeom>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TFEGeom>);
	this->set_add_rhs_elem_fct(  id, &T::template add_rhs_elem<TElem, TFEGeom>);

//	set computation of linearized defect w.r.t velocity
	m_imVelocity. 		set_fct(id, this, &T::template lin_def_velocity<TElem, TFEGeom>);
	m_imDiffusion.		set_fct(id, this, &T::template lin_def_diffusion<TElem, TFEGeom>);
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

