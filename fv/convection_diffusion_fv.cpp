/*
 * convection_diffusion_fv.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#include "convection_diffusion_fv.h"

#include "common/util/provider.h"
#include "lib_disc/spatial_disc/disc_util/fvho_geom.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusionFV<TDomain>::
ConvectionDiffusionFV(const char* functions, const char* subsets)
 : ConvectionDiffusionBase<TDomain>(functions,subsets),
   m_order(1), m_lfeID(LFEID::LAGRANGE, m_order),
   m_bQuadOrderUserDef(false), m_quadOrderSCV(m_order+1), m_quadOrderSCVF(m_order+1)
{
	register_all_funcs(m_order, m_quadOrderSCV, m_quadOrderSCVF);
}

template<typename TDomain>
bool ConvectionDiffusionFV<TDomain>::
request_finite_element_id(const std::vector<LFEID>& vLfeID)
{
//	check number
	if(vLfeID.size() != 1)
	{
		UG_LOG("ERROR in 'ConvectionDiffusion::request_finite_element_id':"
				" Wrong number of functions given. Need exactly "<<1<<"\n");
		return false;
	}

	if(vLfeID[0].type() != LFEID::LAGRANGE)
	{
		UG_LOG("ERROR in 'ConvectionDiffusion::request_finite_element_id':"
				" FV Scheme only implemented for 1st order.\n");
		return false;
	}

//	check that not ADAPTIVE
	if(vLfeID[0].order() < 1)
	{
		UG_LOG("ERROR in 'ConvectionDiffusion::request_finite_element_id':"
				" Adaptive or invalid order not implemented.\n");
		return false;
	}

//	set order
	m_lfeID = vLfeID[0];
	m_order = vLfeID[0].order();
	if(!m_bQuadOrderUserDef) {
		m_quadOrderSCV = m_order+1;
		m_quadOrderSCVF = m_order+1;
	}

	register_all_funcs(m_order, m_quadOrderSCV, m_quadOrderSCVF);

//	is supported
	return true;
}

template<typename TDomain>
bool ConvectionDiffusionFV<TDomain>::
request_non_regular_grid(bool bNonRegular)
{
//	non-regular not supported
	if(bNonRegular) return false;
	else return true;
}

template<typename TDomain>
bool ConvectionDiffusionFV<TDomain>::
use_hanging() const
{
	return true;
}

template<typename TDomain>
void ConvectionDiffusionFV<TDomain>::
set_quad_order(size_t order)
{
	set_quad_order_scv(order);
	set_quad_order_scvf(order);
}

template<typename TDomain>
void ConvectionDiffusionFV<TDomain>::
set_quad_order_scv(size_t order)
{
	m_quadOrderSCV = order; m_bQuadOrderUserDef = true;
	register_all_funcs(m_order, m_quadOrderSCV, m_quadOrderSCVF);
}

template<typename TDomain>
void ConvectionDiffusionFV<TDomain>::
set_quad_order_scvf(size_t order)
{
	m_quadOrderSCVF = order; m_bQuadOrderUserDef = true;
	register_all_funcs(m_order, m_quadOrderSCV, m_quadOrderSCVF);
}

////////////////////////////////////////////////////////////////////////////////
// Assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
prep_elem_loop()
{
//	reference dimension
	static const int refDim = reference_element_traits<TElem>::dim;
	typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
	static const ReferenceObjectID roid = reference_element_type::REFERENCE_OBJECT_ID;

//	request geometry
	static typename TGeomProvider::Type& geo = TGeomProvider::get();

	try{
		geo.update_local(roid, m_order, m_quadOrderSCVF, m_quadOrderSCV);
	}
	UG_CATCH_THROW("ConvectionDiffusion::prep_elem_loop:"
						" Cannot update Finite Volume Geometry.");

//	set local positions
	if(!TGeomProvider::Type::usesHangingNodes)
	{
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
		m_imMass.template	 	set_local_ips<refDim>(geo.scv_local_ips(),
		                       	                      geo.num_scv_ips(), false);
	}
}

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
prep_elem(TElem* elem, const LocalVector& u)
{
//	get reference elements
	static const int refDim = reference_element_traits<TElem>::dim;

//	get corners
	m_vCornerCoords = this->template element_corners<TElem>(elem);

//	request geometry
	static typename TGeomProvider::Type& geo = TGeomProvider::get();

	try{
		geo.update(elem, &m_vCornerCoords[0], &(this->subset_handler()));
	}
	UG_CATCH_THROW("ConvectionDiffusion::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	set local positions
	if(TGeomProvider::Type::usesHangingNodes)
	{
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
		m_imMass.template 		set_local_ips<refDim>(geo.scv_local_ips(),
		                       	                      geo.num_scv_ips());
	}

//	set global positions
	m_imDiffusion.	set_global_ips(geo.scvf_global_ips(), geo.num_scvf_ips());
	m_imVelocity.	set_global_ips(geo.scvf_global_ips(), geo.num_scvf_ips());
	m_imSource.		set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
	m_imReactionRate.	set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
	m_imReaction.	set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
	m_imMassScale.	set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
	m_imMass.		set_global_ips(geo.scv_global_ips(), geo.num_scv_ips());
}

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

//	Diff. Tensor times Gradient
	MathVector<dim> Dgrad;

//	Diffusion and Velocity Term
	if(m_imDiffusion.data_given() || m_imVelocity.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t i = 0, ipCnt = 0; i < geo.num_scvf(); ++i)
		{
		// 	get current SCVF
			const typename TGeomProvider::Type::SCVF& scvf = geo.scvf(i);

		//	loop integration points
			for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
			{
			////////////////////////////////////////////////////
			// Diffusive Term
			////////////////////////////////////////////////////
				if(m_imDiffusion.data_given())
				{
				// 	loop shape functions
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					{
					// 	Compute Diffusion Tensor times Gradient
						MatVecMult(Dgrad, m_imDiffusion[ipCnt], scvf.global_grad(ip, sh));

					//	Compute flux at IP
						const number D_diff_flux = VecDot(Dgrad, scvf.normal());

					// 	Add flux term to local matrix
						J(_C_, scvf.from(), _C_, sh) -= D_diff_flux * scvf.weight(ip);
						J(_C_, scvf.to()  , _C_, sh) += D_diff_flux * scvf.weight(ip);
					}
				}

			////////////////////////////////////////////////////
			// Convective Term
			////////////////////////////////////////////////////
				if(m_imVelocity.data_given())
				{
				//	Add Flux contribution
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					{
						const number D_conv_flux = VecDot(m_imVelocity[ipCnt], scvf.normal())
													* scvf.shape(ip, sh);

					//	Add fkux term to local matrix
						J(_C_, scvf.from(), _C_, sh) += D_conv_flux * scvf.weight(ip);
						J(_C_, scvf.to(),   _C_, sh) -= D_conv_flux * scvf.weight(ip);
					}
				}

				ipCnt++;
			} // end loop ip
		} // end loop scvf
	} // end data given

////////////////////////////////////////////////////
// Reaction Term
////////////////////////////////////////////////////

//	if no data for reaction, return
	if(!m_imReactionRate.data_given()) return;

// 	loop Sub Control Volume (SCV)
	for(size_t i = 0, ipOffset = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

	// 	get associated node
		const int co = scv.node_id();

	//	loop shapes
		for(size_t sh = 0; sh < scv.num_sh(); ++sh)
		{
		//	reset integral
			number integral = 0;

		//	loop integration points
			for(size_t ip = 0; ip < scv.num_ip(); ++ip)
			{
				integral += m_imReactionRate[ipOffset+ip] * scv.shape(ip, sh) * scv.weight(ip);
			}

		// 	Add to local matrix
			J(_C_, co, _C_, sh) += integral;
		}

	//	increase ip offset
		ipOffset += scv.num_ip();
	}

//	no explicit dependency in m_imReaction
}


template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0, ipOffset = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TGeomProvider::Type::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int co = scv.node_id();

	//	loop shapes
		for(size_t sh = 0; sh < scv.num_sh(); ++sh)
		{
		//	reset integral
			number integral = 0;

		//	loop integration points
			for(size_t ip = 0; ip < scv.num_ip(); ++ip)
			{
				if(m_imMassScale.data_given())
					integral += scv.shape(ip, sh) * scv.weight(ip) * m_imMassScale[ipOffset+ip];
				else
					integral += scv.shape(ip, sh) * scv.weight(ip);

				//	no explicit dependency in m_imMass
			}

		// 	Add to local matrix
			J(_C_, co, _C_, sh) += integral;
		}

	//	increase ip offset
		ipOffset += scv.num_ip();
	}
}


template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

	if(m_imDiffusion.data_given() || m_imVelocity.data_given())
	{
	// 	loop Sub Control Volume Faces (SCVF)
		for(size_t i = 0, ipCnt = 0; i < geo.num_scvf(); ++i)
		{
		// 	get current SCVF
			const typename TGeomProvider::Type::SCVF& scvf = geo.scvf(i);

		//	the flux of the scvf
			number flux = 0;

		//	loop integration points
			for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
			{
				number fluxIP = 0;

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
						VecScaleAppend(grad_c, u(_C_,sh), scvf.global_grad(ip, sh));

				//	scale by diffusion tensor
					MatVecMult(Dgrad_c, m_imDiffusion[ipCnt], grad_c);

				// 	Compute flux
					fluxIP = -VecDot(Dgrad_c, scvf.normal());
				}

			/////////////////////////////////////////////////////
			// Convective Term
			/////////////////////////////////////////////////////
				if(m_imVelocity.data_given())
				{
				//	sum up solution
					number solIP = 0;
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						solIP += u(_C_, sh) * scvf.shape(ip, sh);

				//	add convective flux
					fluxIP += solIP * VecDot(m_imVelocity[ipCnt++], scvf.normal());
				}

			//	sum flux
				flux += fluxIP * scvf.weight(ip);
			} // end loop ip

		//	no multiplication with volume is needed, since already contained
		//	in the normal

		//  add to local defect
			d(_C_, scvf.from()) += flux;
			d(_C_, scvf.to()  ) -= flux;

		} // end loop scvf
	} // end data switch

//	reaction rate
	if(m_imReactionRate.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t i = 0, ipCnt = 0; i < geo.num_scv(); ++i)
		{
		// 	get current SCV
			const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

		//	reset integral
			number integral = 0;

		//	loop integration points
			for(size_t ip = 0; ip < scv.num_ip(); ++ip)
			{
			//	compute solution at ip
				number solIP = 0;
				for(size_t sh = 0; sh < scv.num_sh(); ++sh)
					solIP += u(_C_, sh) * scv.shape(ip, sh);

			//	add to integral-sum
				integral += m_imReactionRate[ipCnt++] * solIP * scv.weight(ip);
			}

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_C_, co) += integral;
		}
	}

//	reaction rate
	if(m_imReaction.data_given())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t i = 0, ipCnt = 0; i < geo.num_scv(); ++i)
		{
		// 	get current SCV
			const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

		//	reset integral
			number integral = 0;

		//	loop integration points
			for(size_t ip = 0; ip < scv.num_ip(); ++ip)
			{
			//	add to integral-sum
				integral += m_imReaction[ipCnt++] * scv.weight(ip);
			}

		// 	get associated node
			const int co = scv.node_id();

		// 	Add to local defect
			d(_C_, co) += integral;
		}
	}
}


template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0, ipCnt = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

	//	reset integral
		number integral = 0;

	//	loop integration points
		for(size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
		//	compute solution at ip
			number solIP = 0;
			for(size_t sh = 0; sh < scv.num_sh(); ++sh)
				solIP += u(_C_, sh) * scv.shape(ip, sh);

		//	compute value
			number val = solIP;

		//	multiply by scaling
			if(m_imMassScale.data_given())
				val *= m_imMassScale[ipCnt];

		//	add mass
			if(m_imMass.data_given())
				val += m_imMass[ipCnt];

		//	next ip
			ipCnt++;

		//	add to integral-sum
			integral += val * scv.weight(ip);
		}

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		d(_C_, co) +=  integral;
	}
}


template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
add_rhs_elem(LocalVector& d)
{
//	if zero data given, return
	if(!m_imSource.data_given()) return;

//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0, ipCnt = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

	//	reset integral
		number integral = 0;

	//	loop integration points
		for(size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
		//	add to integral-sum
			integral += m_imSource[ipCnt++] * scv.weight(ip);
		}

	// 	get associated node
		const int co = scv.node_id();

	// 	Add to local defect
		d(_C_, co) +=  integral;
	}
}


//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
lin_def_velocity(const LocalVector& u,
                     std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                     const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0, ipCnt = 0; i < geo.num_scvf(); ++i)
	{
	// 	get current SCVF
		const typename TGeomProvider::Type::SCVF& scvf = geo.scvf(i);

		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
		// 	compute shape at ip
			number shape_u = 0.0;
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				shape_u += u(_C_,sh) * scvf.shape(ip, sh);

		//	add parts for both sides of scvf
			VecScale(vvvLinDef[ip][_C_][scvf.from()], scvf.normal(), shape_u);
			VecScale(vvvLinDef[ip][_C_][scvf.to()  ], scvf.normal(), -shape_u);

			++ipCnt;
		}
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
lin_def_diffusion(const LocalVector& u,
                      std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
                      const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

//  loop Sub Control Volume Faces (SCVF)
	for(size_t i = 0, ipCnt = 0; i < geo.num_scvf(); ++i)
	{
	// 	get current SCVF
		const typename TGeomProvider::Type::SCVF& scvf = geo.scvf(i);

		for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
		{
		// 	compute gradient at ip
			MathVector<dim> grad_u;	VecSet(grad_u, 0.0);
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				VecScaleAppend(grad_u, u(_C_,sh), scvf.global_grad(ip, sh));

		//	part coming from -\nabla u * \vec{n}
			for(size_t k=0; k < (size_t)dim; ++k)
				for(size_t j = 0; j < (size_t)dim; ++j)
				{
					const number val = (scvf.normal())[j] * grad_u[k];

					vvvLinDef[ipCnt][_C_][scvf.from()](k,j) = val;
					vvvLinDef[ipCnt][_C_][scvf.to()  ](k,j) -= val;
				}

			++ipCnt;
		}
	}
}

//	computes the linearized defect w.r.t to the reaction rate
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
lin_def_reaction_rate(const LocalVector& u,
                           std::vector<std::vector<number> > vvvLinDef[],
                           const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0, ipCnt = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

	// 	get associated node
		const int co = scv.node_id();

	//	loop integration points
		for(size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
		//	compute solution at ip
			number solIP = 0;
			for(size_t sh = 0; sh < scv.num_sh(); ++sh)
				solIP += u(_C_, sh) * scv.shape(ip, sh);

		// 	set lin defect
			vvvLinDef[ipCnt++][_C_][co] = solIP * scv.weight(ip);
		}
	}
}

//	computes the linearized defect w.r.t to the reaction
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
lin_def_reaction(const LocalVector& u,
                     std::vector<std::vector<number> > vvvLinDef[],
                     const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0, ipCnt = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

	// 	get associated node
		const int co = scv.node_id();

	//	loop integration points
		for(size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
		// 	set lin defect
			vvvLinDef[ipCnt++][_C_][co] = scv.weight(ip);
		}
	}
}

//	computes the linearized defect w.r.t to the source
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
lin_def_source(const LocalVector& u,
                   std::vector<std::vector<number> > vvvLinDef[],
                   const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0, ipCnt = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

	// 	get associated node
		const int co = scv.node_id();

	//	loop integration points
		for(size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
		// 	set lin defect
			vvvLinDef[ipCnt++][_C_][co] = scv.weight(ip);
		}
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
lin_def_mass_scale(const LocalVector& u,
                       std::vector<std::vector<number> > vvvLinDef[],
                       const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0, ipCnt = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

	// 	get associated node
		const int co = scv.node_id();

	//	loop integration points
		for(size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
		//	compute solution at ip
			number solIP = 0;
			for(size_t sh = 0; sh < scv.num_sh(); ++sh)
				solIP += u(_C_, sh) * scv.shape(ip, sh);

		// 	set lin defect
			vvvLinDef[ipCnt++][_C_][co] = solIP * scv.weight(ip);
		}
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
lin_def_mass(const LocalVector& u,
                  std::vector<std::vector<number> > vvvLinDef[],
                  const size_t nip)
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

// 	loop Sub Control Volumes (SCV)
	for(size_t i = 0, ipCnt = 0; i < geo.num_scv(); ++i)
	{
	// 	get current SCV
		const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

	// 	get associated node
		const int co = scv.node_id();

	//	loop integration points
		for(size_t ip = 0; ip < scv.num_ip(); ++ip)
		{
		// 	set lin defect
			vvvLinDef[ipCnt++][_C_][co] = scv.weight(ip);
		}
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
ex_value(const LocalVector& u,
         const MathVector<dim> vGlobIP[],
         const MathVector<TGeomProvider::Type::dim> vLocIP[],
         const size_t nip,
         number vValue[],
         bool bDeriv,
         std::vector<std::vector<number> > vvvDeriv[])
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;

//	reference object id
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	FV SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//  loop Sub Control Volume Faces (SCVF)
		for(size_t i = 0, ipCnt = 0; i < geo.num_scvf(); ++i)
		{
		// 	get current SCVF
			const typename TGeomProvider::Type::SCVF& scvf = geo.scvf(i);

			for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
			{
			//	compute concentration at ip
				vValue[ipCnt] = 0.0;
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					vValue[ipCnt] += u(_C_, sh) * scvf.shape(ip, sh);

			//	compute derivative w.r.t. to unknowns iff needed
				if(bDeriv)
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						vvvDeriv[ipCnt][_C_][sh] = scvf.shape(ip, sh);

				++ipCnt;
			}
		}
	}
//	FV SCV ip
	else if(vLocIP == geo.scv_local_ips())
	{
	// 	loop Sub Control Volumes (SCV)
		for(size_t i = 0, ipCnt = 0; i < geo.num_scv(); ++i)
		{
		// 	get current SCV
			const typename TGeomProvider::Type::SCV& scv = geo.scv(i);

		//	loop integration points
			for(size_t ip = 0; ip < scv.num_ip(); ++ip)
			{
			//	compute solution at ip
				vValue[ipCnt] = 0.0;
				for(size_t sh = 0; sh < scv.num_sh(); ++sh)
					vValue[ipCnt] += u(_C_, sh) * scv.shape(ip, sh);

			//	compute derivative w.r.t. to unknowns iff needed
				if(bDeriv)
					for(size_t sh = 0; sh < scv.num_sh(); ++sh)
						vvvDeriv[ipCnt][_C_][sh] = scv.shape(ip, sh);

				++ipCnt;
			}
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

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
ex_grad(const LocalVector& u,
        const MathVector<dim> vGlobIP[],
        const MathVector<TGeomProvider::Type::dim> vLocIP[],
        const size_t nip,
        MathVector<dim> vValue[],
        bool bDeriv,
        std::vector<std::vector<MathVector<dim> > > vvvDeriv[])
{
//	request geometry
	static const typename TGeomProvider::Type& geo = TGeomProvider::get();

//	reference element
	typedef typename reference_element_traits<TElem>::reference_element_type
			ref_elem_type;

//	reference dimension
	static const int refDim = ref_elem_type::dim;

//	number of shape functions
	static const size_t numSH =	ref_elem_type::numCorners;

//	reference object id
	static const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	FV SCVF ip
	if(vLocIP == geo.scvf_local_ips())
	{
	//  loop Sub Control Volume Faces (SCVF)
		for(size_t i = 0, ipCnt = 0; i < geo.num_scvf(); ++i)
		{
		// 	get current SCVF
			const typename TGeomProvider::Type::SCVF& scvf = geo.scvf(i);

			for(size_t ip = 0; ip < scvf.num_ip(); ++ip)
			{
				VecSet(vValue[ipCnt], 0.0);
				for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
					VecScaleAppend(vValue[ipCnt], u(_C_, sh), scvf.global_grad(ip, sh));

				if(bDeriv)
					for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
						vvvDeriv[ipCnt][_C_][sh] = scvf.global_grad(ip, sh);

				++ipCnt;
			}
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
		UG_CATCH_THROW("ConvectionDiffusion::ex_grad: trial space missing");
	}
};

////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

// register for all dim
template<>
void ConvectionDiffusionFV<Domain1d>::
register_all_funcs(int order, int quadOrderSCV, int quadOrderSCVF)
{
//	Edge
	switch(order)
	{
/*		case 1:	{typedef FVGeometry<1, Edge, dim> FVGeom;
				 register_func<Edge, Provider<FVGeom> >(); break;}
		case 2:	{typedef FVGeometry<2, Edge, dim> FVGeom;
				 register_func<Edge, Provider<FVGeom> >(); break;}
		case 3:	{typedef FVGeometry<3, Edge, dim> FVGeom;
				 register_func<Edge, Provider<FVGeom> >(); break;}
		default: {typedef DimFVGeometry<1, dim> FVGeom;
		 	 	 register_func<Edge, Provider<FVGeom> >(); break;}
*/	}
}

// register for all dim
template<>
void ConvectionDiffusionFV<Domain2d>::
register_all_funcs(int order, int quadOrderSCV, int quadOrderSCVF)
{
	if(quadOrderSCV == order+1 && quadOrderSCVF ==  order+1)
	{
	//	Triangle
		switch(order)
		{
			case 1:	{typedef FVGeometry<1, Triangle, dim> FVGeom;
					 register_func<Triangle, Provider<FVGeom> >(); break;}
			case 2:	{typedef FVGeometry<2, Triangle, dim> FVGeom;
					 register_func<Triangle, Provider<FVGeom> >(); break;}
			case 3:	{typedef FVGeometry<3, Triangle, dim> FVGeom;
					 register_func<Triangle, Provider<FVGeom> >(); break;}
			default: {typedef DimFVGeometry<2, dim> FVGeom;
					 register_func<Triangle, FlexGeomProvider<FVGeom> >(); break;}
		}

	//	Quadrilateral
		switch(order) {
			case 1:	{typedef FVGeometry<1, Quadrilateral, dim> FVGeom;
					 register_func<Quadrilateral, Provider<FVGeom> >(); break;}
			case 2:	{typedef FVGeometry<2, Quadrilateral, dim> FVGeom;
					 register_func<Quadrilateral, Provider<FVGeom> >(); break;}
			case 3:	{typedef FVGeometry<3, Quadrilateral, dim> FVGeom;
					 register_func<Quadrilateral, Provider<FVGeom> >(); break;}
			default: {typedef DimFVGeometry<2, dim> FVGeom;
					  register_func<Quadrilateral, FlexGeomProvider<FVGeom> >(); break;}
		}
	}
	else
	{
		typedef DimFVGeometry<2, dim> FVGeom;
		register_func<Triangle, FlexGeomProvider<FVGeom> >();
		register_func<Quadrilateral, FlexGeomProvider<FVGeom> >();
	}
}

// register for all dim
template<>
void ConvectionDiffusionFV<Domain3d>::
register_all_funcs(int order, int quadOrderSCV, int quadOrderSCVF)
{
	if(quadOrderSCV == order+1 && quadOrderSCVF ==  order+1)
	{
	//	Tetrahedron
		switch(order)
		{
			case 1:	{typedef FVGeometry<1, Tetrahedron, dim> FVGeom;
					 register_func<Tetrahedron, Provider<FVGeom> >(); break;}
			case 2:	{typedef FVGeometry<2, Tetrahedron, dim> FVGeom;
					 register_func<Tetrahedron, Provider<FVGeom> >(); break;}
			case 3:	{typedef FVGeometry<3, Tetrahedron, dim> FVGeom;
					 register_func<Tetrahedron, Provider<FVGeom> >(); break;}
			default: {typedef DimFVGeometry<3, dim> FVGeom;
					  register_func<Tetrahedron, FlexGeomProvider<FVGeom> >(); break;}
		}

	//	Prism
		switch(order) {
			case 1:	{typedef FVGeometry<1, Prism, dim> FVGeom;
					 register_func<Prism, Provider<FVGeom> >(); break;}
			default: {typedef DimFVGeometry<3, dim> FVGeom;
					  register_func<Prism, FlexGeomProvider<FVGeom> >(); break;}
		}

	//	Hexahedron
		switch(order)
		{
			case 1:	{typedef FVGeometry<1, Hexahedron, dim> FVGeom;
					 register_func<Hexahedron, Provider<FVGeom> >(); break;}
			case 2:	{typedef FVGeometry<2, Hexahedron, dim> FVGeom;
					 register_func<Hexahedron, Provider<FVGeom> >(); break;}
			case 3:	{typedef FVGeometry<3, Hexahedron, dim> FVGeom;
					 register_func<Hexahedron, Provider<FVGeom> >(); break;}
			default: {typedef DimFVGeometry<3, dim> FVGeom;
					  register_func<Hexahedron, FlexGeomProvider<FVGeom> >(); break;}
		}
	}
	else
	{
		typedef DimFVGeometry<3, dim> FVGeom;
		register_func<Tetrahedron, FlexGeomProvider<FVGeom> >();
		register_func<Prism, FlexGeomProvider<FVGeom> >();
		register_func<Hexahedron, FlexGeomProvider<FVGeom> >();
	}

}

template<typename TDomain>
template<typename TElem, typename TGeomProvider>
void ConvectionDiffusionFV<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef this_type T;
	static const int refDim = reference_element_traits<TElem>::dim;

	this->enable_fast_add_elem(true);
	this->set_prep_elem_loop_fct(id, &T::template prep_elem_loop<TElem, TGeomProvider>);
	this->set_prep_elem_fct(	 id, &T::template prep_elem<TElem, TGeomProvider>);
	this->set_fsh_elem_loop_fct( id, &T::template fsh_elem_loop<TElem, TGeomProvider>);
	this->set_add_jac_A_elem_fct(id, &T::template add_jac_A_elem<TElem, TGeomProvider>);
	this->set_add_jac_M_elem_fct(id, &T::template add_jac_M_elem<TElem, TGeomProvider>);
	this->set_add_def_A_elem_fct(id, &T::template add_def_A_elem<TElem, TGeomProvider>);
	this->set_add_def_M_elem_fct(id, &T::template add_def_M_elem<TElem, TGeomProvider>);
	this->set_add_rhs_elem_fct(  id, &T::template add_rhs_elem<TElem, TGeomProvider>);

//	set computation of linearized defect w.r.t velocity
	m_imVelocity. set_fct(id, this, &T::template lin_def_velocity<TElem, TGeomProvider>);
	m_imDiffusion.set_fct(id, this, &T::template lin_def_diffusion<TElem, TGeomProvider>);
	m_imReactionRate. set_fct(id, this, &T::template lin_def_reaction_rate<TElem, TGeomProvider>);
	m_imReaction. set_fct(id, this, &T::template lin_def_reaction<TElem, TGeomProvider>);
	m_imSource.	  set_fct(id, this, &T::template lin_def_source<TElem, TGeomProvider>);
	m_imMassScale.set_fct(id, this, &T::template lin_def_mass_scale<TElem, TGeomProvider>);
	m_imMass.	  set_fct(id, this, &T::template lin_def_mass<TElem, TGeomProvider>);

//	exports
	m_exValue->	   template set_fct<T,refDim>(id, this, &T::template ex_value<TElem, TGeomProvider>);
	m_exGrad->template set_fct<T,refDim>(id, this, &T::template ex_grad<TElem, TGeomProvider>);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ConvectionDiffusionFV<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConvectionDiffusionFV<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionFV<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug

