/*
 * convection_diffusion_fv1.cpp
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#include "convection_diffusion_fv1.h"

#include "lib_disc/spatial_disc/disc_util/fv1Cut_geom.h"

#include "lib_disc/spatial_disc/disc_util/geom_provider.h"
//#include "lib_disc/spatial_disc/disc_util/fv1_geom.h"
//#include "lib_disc/spatial_disc/disc_util/hfv1_geom.h"
#include "lib_disc/spatial_disc/disc_util/conv_shape.h"

namespace ug{
namespace ConvectionDiffusionPlugin{

////////////////////////////////////////////////////////////////////////////////
//	general
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusionFV1<TDomain>::
ConvectionDiffusionFV1(const char* functions, const char* subsets)
 : ConvectionDiffusionBase<TDomain>(functions,subsets),
   m_spConvShape(new ConvectionShapesNoUpwind<dim>),
   m_bNonRegularGrid(false)
{
	register_all_funcs(m_bNonRegularGrid);
}


template<typename TDomain>
void ConvectionDiffusionFV1<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid, bool bStaticRoid)
{
//	check number
	if(vLfeID.size() != 1)
		UG_THROW("ConvectionDiffusion: Wrong number of functions given. "
				"Need exactly "<<1);

	if(vLfeID[0].order() != 1 || vLfeID[0].type() != LFEID::LAGRANGE)
		UG_THROW("ConvectionDiffusion FV Scheme only implemented for 1st order.");

//	remember
	m_bNonRegularGrid = bNonRegularGrid;

	m_LFEID = vLfeID[0];

//	update assemble functions
	register_all_funcs(m_bNonRegularGrid);
}

template<typename TDomain>
bool ConvectionDiffusionFV1<TDomain>::
use_hanging() const
{
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// Assembling functions
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
prep_elem_loop(const ReferenceObjectID roid, const int si)
{
	// 	Only first order implementation
	if(!(TFVGeom::order == 1))
		UG_THROW("Only first order implementation, but other Finite Volume"
				" Geometry set.");

//	check, that upwind has been set
	if(m_spConvShape.invalid())
		UG_THROW("ConvectionDiffusionFV1::prep_elem_loop:"
						" Upwind has not been set.");

//	set local positions
	if(!TFVGeom::usesHangingNodes && TFVGeom::staticLocalData)
	{
		static const int refDim = TElem::dim;
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
//		TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);
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
			UG_THROW("ConvectionDiffusionFV1::prep_elem_loop:"
						" Cannot init upwind for element type.");
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
fsh_elem_loop()
{}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{
// 	Update Geometry for this element
	//static TFVGeom& geo = GeomProvider<TFVGeom>::get();
//	TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID,1);
    TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1),1);

//	TFVGeom& geo = GeomProvider<TFVGeom>::get();
//	TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);

    // fix: set orientation initially globally!
	geo.set_orientation(1);

	try{
		geo.update(elem, roid, vCornerCoords, &(this->subset_handler()));
	}UG_CATCH_THROW("ConvectionDiffusionFV1::prep_elem:"
						" Cannot update Finite Volume Geometry.");

//	set local positions
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
				UG_THROW("ConvectionDiffusionFV1::prep_elem_loop:"
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

void LU(MathMatrix<3, 3>& R, MathMatrix<3, 3>& L, MathMatrix<3, 3> A)
{
	// n-1 Iterationsschritte
	for ( size_t i = 1; i < 2; ++i)
	{
		for ( size_t k = i+1; i < 3; ++i)
		{
			L[k][i] = R[k][i] / R[i][i];

			for ( size_t j = i; i < 3; ++i)
				R[k][j] = R[k][j] - L[k][i] * R[i][j];
		}
	}

}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{

	bool debug = false;
	bool boundary = false;

	// get finite volume geometry
	//	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID,1);
	const bool bElementIsCut = geo.get_element_modus();
	const bool bElementIsOutside = geo.get_boolian_for_diffusion();

	// First call with orientation = 1:
	int orientation = 1;
	geo.set_orientation(orientation);

	geo.init_integral();

	UG_LOG("------------------> jac = " << geo.get_integral() << "\n");

// normal assembling if not cut by interface:
	if ( !bElementIsCut )
	{
		LocalVector dummyU;
		LocalIndices ind = u.get_indices();
		dummyU.resize(ind);
		dummyU = 0;

		this->template add_jac_A_elem_local<TElem,TFVGeom> (geo, J, u, dummyU, elem, vCornerCoords, bElementIsOutside);
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

		this->template add_jac_A_elem_local<TElem,TFVGeom> (geo, locJ_tri, locU_tri, jump_tri, elem, vCornerCoords, false);

		if ( boundary )
		{
			geo.reset_jacobian_on_interface(locJ_tri, 3);
			this->template add_jac_A_elem_boundary<TElem,TFVGeom> (geo, locJ_tri, locU_tri, elem, vCornerCoords, false);
		}

		geo.set_jacobian_tri(locJ_tri);
		geo.set_DoF_tag_tri(false);
 	}
	if ( roidCheck == ROID_QUADRILATERAL )
	{
		geo.set_local_sol(locU_quad, 4, u, orientation);
		LocalVector jump_quad = geo.set_jump_values(ind, 4);

		this->template add_jac_A_elem_local<TElem,TFVGeom> (geo, locJ_quad, locU_quad, jump_quad, elem, vCornerCoords, false);

		if ( boundary )
		{
			geo.reset_jacobian_on_interface(locJ_quad, 4);
			this->template add_jac_A_elem_boundary<TElem,TFVGeom> (geo, locJ_quad, locU_quad, elem, vCornerCoords, false);
		}

		geo.set_jacobian_quad(locJ_quad);
		geo.set_DoF_tag_quad(false);
	}


// Second call with orientation = 1:
	bool shiftTag = geo.get_bScaleDoFs(); // shiftTag = true in case of double DoFs on interface!

	orientation *= -1;
	geo.set_orientation(orientation);
	try{
		geo.update(elem, ROID_UNKNOWN, vCornerCoords, &(this->subset_handler()));
	}UG_CATCH_THROW("ConvectionDiffusionFV1::update:"
						" Cannot update Finite Volume Geometry.");

	if ( debug ) geo.print_InterfaceIDdata();

	roidCheck = geo.get_roid();
	if ( roidCheck == ROID_TRIANGLE )
	{
		geo.set_local_sol(locU_tri, 3, u, orientation);
		LocalVector jump_tri = geo.set_jump_values(ind, 3);

		this->template add_jac_A_elem_local<TElem,TFVGeom> (geo, locJ_tri, locU_tri, jump_tri, elem, vCornerCoords, true);

		if ( boundary )
		{
			geo.reset_jacobian_on_interface(locJ_tri, 3);
			this->template add_jac_A_elem_boundary<TElem,TFVGeom> (geo, locJ_tri, locU_tri, elem, vCornerCoords, true);
		}

		geo.set_jacobian_tri(locJ_tri);
		geo.set_DoF_tag_tri(shiftTag);

	}
	if ( roidCheck == ROID_QUADRILATERAL )
	{
		geo.set_local_sol(locU_quad, 4, u, orientation);
		LocalVector jump_quad = geo.set_jump_values(ind, 4);

		this->template add_jac_A_elem_local<TElem,TFVGeom> (geo, locJ_quad, locU_quad, jump_quad, elem, vCornerCoords, true);

		if ( boundary )
		{
			geo.reset_jacobian_on_interface(locJ_quad, 4);
			this->template add_jac_A_elem_boundary<TElem,TFVGeom> (geo, locJ_quad, locU_quad, elem, vCornerCoords, true);
		}

        geo.set_jacobian_quad(locJ_quad);
		geo.set_DoF_tag_quad(shiftTag);
 	}


}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
add_jac_A_elem_local(TFVGeom& geo, LocalMatrix& J, const LocalVector& u, const LocalVector& jump, GridObject* elem, const MathVector<dim> vCornerCoords[], const bool bElementIsOutside)
{
	ug::MathMatrix<dim, dim, double> diffusion = m_imDiffusion[0];
	diffusion *= 1.0;

	if ( bElementIsOutside ) // = inside circle line!!
		diffusion *= 1.0;

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
		// 	get current SCV
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
					MatVecMult(Dgrad, diffusion, scvf.global_grad(sh));

				//	Compute flux at IP
					const number D_diff_flux = VecDot(Dgrad, scvf.normal());

					J(_C_, scvf.from(), _C_, sh) -= D_diff_flux;
					J(_C_, scvf.to()  , _C_, sh) += D_diff_flux;
				}
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
void ConvectionDiffusionFV1<TDomain>::
add_jac_A_elem_boundary(TFVGeom& geo, LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const bool bElementIsOutside)
{
 	double diffusion = 10.0;

	if ( bElementIsOutside )
		diffusion = 1.0;

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
void ConvectionDiffusionFV1<TDomain>::
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
void ConvectionDiffusionFV1<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
    bool output = false;
    bool output_integral = true;

	bool debug = false;
	bool boundary = false;
	bool add = true;

	// get finite volume geometry
	//	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();
	TFVGeom& geo = GeomProvider<TFVGeom>::get(m_LFEID,1);
	const bool bElementIsCut = geo.get_element_modus();
	const bool bElementIsOutside = geo.get_boolian_for_diffusion();

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
		LocalVector source = geo.set_source(imSource, ind, 3, true);

		if ( output )
		{
        for ( size_t i = 0; i < 3; ++i)
            UG_LOG("*** corner" << vCornerCoords[i][0] << " and " << vCornerCoords[i][1] << "\n" );
        UG_LOG("\n" );
		}

		this->template add_def_A_elem_local<TElem,TFVGeom> (geo, d, u, dummyU, dummyU, source, elem, vCornerCoords, bElementIsOutside);

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
		LocalVector source_tri = geo.set_source(imSource, ind, 3, false);

		if ( output ) UG_LOG(" tri 1: orientaten: " << orientation << "\n");

		this->template add_def_A_elem_local<TElem,TFVGeom> (geo, locD_tri, locU_tri, jump_tri, jump_tri_grad, source_tri, elem, vCornerCoords, false);

		if ( boundary )
		{
			geo.reset_defect_on_interface(locD_tri, 3);
			this->template add_def_A_elem_boundary<TElem,TFVGeom> (geo, locD_tri, locU_tri, elem, vCornerCoords, false);
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
		LocalVector source_quad = geo.set_source(imSource, ind, 4, false);

		if ( output ) UG_LOG(" quad 1: orientaten: " << orientation << "\n");
 		this->template add_def_A_elem_local<TElem,TFVGeom> (geo, locD_quad, locU_quad, jump_quad, jump_quad_grad, source_quad, elem, vCornerCoords, false);

		if ( boundary )
		{
			geo.reset_defect_on_interface(locD_quad, 4);
	 		this->template add_def_A_elem_boundary<TElem,TFVGeom> (geo, locD_quad, locU_quad, elem, vCornerCoords, false);
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
		geo.update(elem, ROID_UNKNOWN, vCornerCoords, &(this->subset_handler()));
	}UG_CATCH_THROW("ConvectionDiffusionFV1::update:"
						" Cannot update Finite Volume Geometry.");

	if ( debug ) geo.print_InterfaceIDdata();

	roidCheck = geo.get_roid();
	if ( roidCheck == ROID_TRIANGLE )
	{
		geo.set_local_sol(locU_tri, 3, u, orientation);
		LocalVector jump_tri = geo.set_jump_values(ind, 3);
		LocalVector jump_tri_grad = geo.set_jump_grad_values(ind, 3);
		LocalVector source_tri = geo.set_source(imSource, ind, 3, false);

		if ( output ) UG_LOG(" tri 2: orientaten: " << orientation << "\n");

		this->template add_def_A_elem_local<TElem,TFVGeom> (geo, locD_tri, locU_tri, jump_tri, jump_tri_grad, source_tri, elem, vCornerCoords, true);

		if ( boundary )
		{
			geo.reset_defect_on_interface(locD_tri, 3);
			this->template add_def_A_elem_boundary<TElem,TFVGeom> (geo, locD_tri, locU_tri, elem, vCornerCoords, true);
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
		LocalVector source_quad = geo.set_source(imSource, ind, 4, false);

		if ( output ) UG_LOG(" quad 2: orientaten: " << orientation << "\n");

 		this->template add_def_A_elem_local<TElem,TFVGeom> (geo, locD_quad, locU_quad, jump_quad, jump_quad_grad, source_quad, elem, vCornerCoords, true);

 		if ( boundary )
		{
			geo.reset_defect_on_interface(locD_quad, 4);
			this->template add_def_A_elem_boundary<TElem,TFVGeom> (geo, locD_quad, locU_quad, elem, vCornerCoords, true);
		}

 		number intValElem = this->template add_l2error_A_elem<TElem,TFVGeom> (geo, ROID_QUADRILATERAL, locD_quad, locU_quad, elem);
 		if ( add ) geo.add_to_integral(intValElem);

 		if ( output_integral ) UG_LOG("------------------> quad2: integral = " << sqrt(geo.get_integral()) << "\n");

 		geo.set_defect_quad(locD_quad);
		geo.set_DoF_tag_quad(shiftTag);

 	}
}

template <int dim>
number get_exact_sol_test(MathVector<dim> position)
{
	return sin(2*M_PI*position[0]) + sin(2*M_PI*position[1]);
}

template <int dim>
number get_exact_sol_Gangl(MathVector<dim> position)
{
	double kappa_2 = 10.0;
	double dist_x = position[0] - 0.1;
	double dist_y = position[1] - 0.2;
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
	double dist_x = position[0] - 0.1;
	double dist_y = position[1] - 0.2;
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
number ConvectionDiffusionFV1<TDomain>::
add_l2error_A_elem(TFVGeom& geo, ReferenceObjectID roid, LocalVector& d, const LocalVector& u, GridObject* elem)
{
	bool output = false;

	number integral = 0;

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
			number exactSolIP = get_exact_sol_FedkiwEx5<dim>(vGlobIP[ip]);

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
void ConvectionDiffusionFV1<TDomain>::
add_def_A_elem_boundary(TFVGeom& geo, LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const bool bElementIsOutside)
{
	double diffusion = 10.0;

	if ( bElementIsOutside )
		diffusion = 1.0;

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


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
add_def_A_elem_local(TFVGeom& geo, LocalVector& d, const LocalVector& u, const LocalVector& jump, const LocalVector& jump_grad, const LocalVector& source, GridObject* elem, const MathVector<dim> vCornerCoords[], const bool bElementIsOutside)
{
	ug::MathMatrix<dim, dim, double> diffusion = m_imDiffusion[0];
	diffusion *= 1.0;
	double diffCoeff = 1.0;

	if ( bElementIsOutside ) // = inside circle line!!
	{	diffusion *= 1.0; diffCoeff = 1.0;}

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

}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
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
void ConvectionDiffusionFV1<TDomain>::
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
void ConvectionDiffusionFV1<TDomain>::
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


////////////////////////////////////
///   error estimation (begin)   ///

//	prepares the loop over all elements of one type for the computation of the error estimator
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
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


//	check that upwind has been set
	if (m_spConvShape.invalid())
		UG_THROW("ConvectionDiffusionFV1::prep_err_est_elem_loop: "
				 "Upwind has not been set.");

//	set local positions
	if (!TFVGeom::usesHangingNodes)
	{
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

		//	init upwind for element type
		TFVGeom& geo = GeomProvider<TFVGeom>::get();
		if (!m_spConvShape->template set_geometry_type<TFVGeom>(geo))
			UG_THROW("ConvectionDiffusionFV1::prep_err_est_elem_loop: "
					 "Cannot init upwind for element type.");

		// store values of shape functions in local IPs
		LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> trialSpace
					= Provider<LagrangeP1<typename reference_element_traits<TElem>::reference_element_type> >::get();

		m_shapeValues.resize(numElemIPs, numSideIPs, trialSpace.num_sh());
		for (size_t ip = 0; ip < numElemIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtElemIP(ip), elemIPs[ip]);
		for (size_t ip = 0; ip < numSideIPs; ip++)
			trialSpace.shapes(m_shapeValues.shapesAtSideIP(ip), sideIPs[ip]);
	}
}

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{}

//	computes the error estimator contribution (stiffness part) for one element
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

//	request geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

////////////////
// SIDE TERMS //
////////////////

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
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

// 	some help variables
	MathVector<dim> fluxDensity, gradC, normal;

// calculate grad u (take grad from first scvf ip (grad u is constant on the entire element))
	if (geo.num_scvf() < 1) {UG_THROW("Element has no SCVFs!");}
	const typename TFVGeom::SCVF& scvf = geo.scvf(0);

	VecSet(gradC, 0.0);
	for (size_t j=0; j<m_shapeValues.num_sh(); j++)
		VecScaleAppend(gradC, u(_C_,j), scvf.global_grad(j));

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

			////// diffusion //////
				if (m_imDiffusion.data_given())
					MatVecScaleMultAppend(fluxDensity, -1.0, m_imDiffusion[0], gradC);

			////// convection //////
				if (m_imVelocity.data_given())
				{
					number val = 0.0;
					for (size_t sh = 0; sh < m_shapeValues.num_sh(); sh++)
						val += u(_C_,sh) * m_shapeValues.shapeAtSideIP(sh,sip);

					VecScaleAppend(fluxDensity, val, m_imVelocity[ip]);
				}

			////// general flux //////
				if (m_imFlux.data_given())
					VecAppend(fluxDensity, m_imFlux[ip]);

				(*err_est_data)(side_list[side],sip) += scale * VecDot(fluxDensity, normal);
			}

			passedIPs += err_est_data->num_side_ips(side_list[side]);
		}
		UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
				<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
	}

//////////////////
// VOLUME TERMS //
//////////////////

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
		{
			number total = 0.0;

		////// diffusion //////	TODO ONLY FOR (PIECEWISE) CONSTANT DIFFUSION TENSOR SO FAR!
		// div(D*grad(c)) = div(v)*u + v*grad(c)
		// nothing to do, as u is piecewise linear and div(D*grad(c)) disappears

		////// convection ////// TODO ONLY FOR CONSTANT VELOCITY FIELDS SO FAR!
		// div(v*c) = div(v)*u + v*grad(c) -- gradC has been calculated above
			if (m_imVelocity.data_given())
				total += VecDot(m_imVelocity[ip], gradC);

		////// general flux ////// TODO ONLY FOR DIVERGENCE-FREE FLUX FIELD SO FAR!
		// nothing to do

		////// reaction //////
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
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
// note: mass parts only enter volume term

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

//	request geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

// 	loop integration points
	try
	{
		for (size_t ip = 0; ip < err_est_data->num_elem_ips(elem->reference_object_id()); ip++)
		{
			number total = 0.0;

		////// mass scale //////
			if (m_imMassScale.data_given())
			{
				number val = 0.0;
				for (size_t sh = 0; sh < geo.num_sh(); sh++)
					val += u(_C_,sh) * m_shapeValues.shapeAtElemIP(sh,ip);

				total += m_imMassScale[ip] * val;
			}

		////// mass //////
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
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	err_est_type* err_est_data = dynamic_cast<err_est_type*>(this->m_spErrEstData.get());

	if (err_est_data->surface_view().get() == NULL) {UG_THROW("Error estimator has NULL surface view.");}
	MultiGrid* pErrEstGrid = (MultiGrid*) (err_est_data->surface_view()->subset_handler()->multi_grid());

////////////////
// SIDE TERMS //
////////////////
//	get the sides of the element
	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::side_type>::secure_container side_list;
	pErrEstGrid->associated_elements_sorted(side_list, (TElem*) elem);
	if (side_list.size() != (size_t) ref_elem_type::numSides)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

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

			////// vector source //////
				if (m_imVectorSource.data_given())
					(*err_est_data)(side_list[side],sip) += scale * VecDot(m_imVectorSource[ip], normal);
			}

			passedIPs += err_est_data->num_side_ips(side_list[side]);
		}
		UG_CATCH_THROW("Values for the error estimator could not be assembled at every IP." << std::endl
				<< "Maybe wrong type of ErrEstData object? This implementation needs: SideAndElemErrEstData.");
	}

//////////////////
// VOLUME TERMS //
//////////////////
	if (!m_imSource.data_given()) return;

	typename MultiGrid::traits<typename SideAndElemErrEstData<TDomain>::elem_type>::secure_container elem_list;
	pErrEstGrid->associated_elements_sorted(elem_list, (TElem*) elem);
	if (elem_list.size() != 1)
		UG_THROW ("Mismatch of numbers of sides in 'ConvectionDiffusionFV1::compute_err_est_elem'");

////// source //////
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
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
fsh_err_est_elem_loop()
{
//	finish the element loop in the same way as the actual discretization
	this->template fsh_elem_loop<TElem, TFVGeom> ();
};

///   error estimation (end)     ///
////////////////////////////////////

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
lin_def_velocity(const LocalVector& u,
                 std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                 const size_t nip)
{
// 	get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

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
			VecScaleAppend(linDefect, u(_C_,sh), convShape.D_vel(ip, sh));

	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] += linDefect;
		vvvLinDef[ip][_C_][scvf.to()] -= linDefect;
	}
}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
lin_def_diffusion(const LocalVector& u,
                  std::vector<std::vector<MathMatrix<dim,dim> > > vvvLinDef[],
                  const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	get conv shapes
	const IConvectionShapes<dim>& convShape = get_updated_conv_shapes(geo);

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
		if(convShape.non_zero_deriv_diffusion())
			for(size_t sh = 0; sh < scvf.num_sh(); ++sh)
				MatAdd(linDefect, convShape.D_diffusion(ip, sh), u(_C_, sh));

	//	add contributions
		vvvLinDef[ip][_C_][scvf.from()] -= linDefect;
		vvvLinDef[ip][_C_][scvf.to()  ] += linDefect;
	}
}

template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
lin_def_flux(const LocalVector& u,
             std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
             const size_t nip)
{
//  get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

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

	//	add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] += scvf.normal();
		vvvLinDef[ip][_C_][scvf.to()] -= scvf.normal();
	}
}
//	computes the linearized defect w.r.t to the reaction rate
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
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
void ConvectionDiffusionFV1<TDomain>::
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
void ConvectionDiffusionFV1<TDomain>::
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

//	computes the linearized defect w.r.t to the vector source
//	(in analogy to velocity)
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
lin_def_vector_source(const LocalVector& u,
                      std::vector<std::vector<MathVector<dim> > > vvvLinDef[],
                      const size_t nip)
{
	// get finite volume geometry
	static const TFVGeom& geo = GeomProvider<TFVGeom>::get();

//	reset the values for the linearized defect
	for(size_t ip = 0; ip < nip; ++ip)
		for(size_t c = 0; c < vvvLinDef[ip].size(); ++c)
			for(size_t sh = 0; sh < vvvLinDef[ip][c].size(); ++sh)
				vvvLinDef[ip][c][sh] = 0.0;

	// loop Sub Control Volumes Faces (SCVF)
	for ( size_t ip = 0; ip < geo.num_scvf(); ++ip ) {
		// get current SCVF
		const typename TFVGeom::SCVF& scvf = geo.scvf( ip );

		// add parts for both sides of scvf
		vvvLinDef[ip][_C_][scvf.from()] -= scvf.normal();
		vvvLinDef[ip][_C_][scvf.to()] += scvf.normal();
	}
}

//	computes the linearized defect w.r.t to the mass scale
template<typename TDomain>
template <typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
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
void ConvectionDiffusionFV1<TDomain>::
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
void ConvectionDiffusionFV1<TDomain>::
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

//	FV1 SCVF ip
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
//	FV1 SCV ip
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
		LagrangeP1<ref_elem_type>& rTrialSpace = Provider<LagrangeP1<ref_elem_type> >::get();

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
void ConvectionDiffusionFV1<TDomain>::
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

//	FV1 SCVF ip
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
void ConvectionDiffusionFV1<TDomain>::
set_upwind(SmartPtr<IConvectionShapes<dim> > shapes) {m_spConvShape = shapes;}

//	computes the linearized defect w.r.t to the velocity
template<typename TDomain>
const typename ConvectionDiffusionFV1<TDomain>::conv_shape_type&
ConvectionDiffusionFV1<TDomain>::
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
void ConvectionDiffusionFV1<Domain1d>::
register_all_funcs(bool bHang)
{
    register_func<RegularEdge, DimFV1CutGeometry<dim, dim, MovingInterfaceBase::InterfaceHandlerLocalDiffusion<dim> > >();


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
void ConvectionDiffusionFV1<Domain2d>::
register_all_funcs(bool bHang)
{
     register_func<Triangle, DimFV1CutGeometry<dim, dim, MovingInterfaceBase::InterfaceHandlerLocalDiffusion<dim> > >();

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
void ConvectionDiffusionFV1<Domain3d>::
register_all_funcs(bool bHang)
{
     register_func<Tetrahedron, DimFV1CutGeometry<dim, dim, MovingInterfaceBase::InterfaceHandlerLocalDiffusion<dim> > >();

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
	}
	*/
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ConvectionDiffusionFV1<TDomain>::
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

// error estimator parts
/*	this->set_prep_err_est_elem_loop(id, &T::template prep_err_est_elem_loop<TElem, TFVGeom>);
	this->set_prep_err_est_elem(id, &T::template prep_err_est_elem<TElem, TFVGeom>);
	this->set_compute_err_est_A_elem(id, &T::template compute_err_est_A_elem<TElem, TFVGeom>);
	this->set_compute_err_est_M_elem(id, &T::template compute_err_est_M_elem<TElem, TFVGeom>);
	this->set_compute_err_est_rhs_elem(id, &T::template compute_err_est_rhs_elem<TElem, TFVGeom>);
	this->set_fsh_err_est_elem_loop(id, &T::template fsh_err_est_elem_loop<TElem, TFVGeom>);
*/
//	set computation of linearized defect w.r.t velocity
	m_imDiffusion.set_fct(id, this, &T::template lin_def_diffusion<TElem, TFVGeom>);
	m_imVelocity. set_fct(id, this, &T::template lin_def_velocity<TElem, TFVGeom>);
	m_imFlux.set_fct(id, this, &T::template lin_def_flux<TElem, TFVGeom>);
	m_imReactionRate. set_fct(id, this, &T::template lin_def_reaction_rate<TElem, TFVGeom>);
	m_imReaction. set_fct(id, this, &T::template lin_def_reaction<TElem, TFVGeom>);
	m_imSource.	  set_fct(id, this, &T::template lin_def_source<TElem, TFVGeom>);
	m_imVectorSource.set_fct(id, this, &T::template lin_def_vector_source<TElem, TFVGeom>);
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
template class ConvectionDiffusionFV1<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConvectionDiffusionFV1<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionFV1<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug

