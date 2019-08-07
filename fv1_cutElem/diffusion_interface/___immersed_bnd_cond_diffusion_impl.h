/*
 * particle_bnd_cond_impl.h
 *
 *  Created on: 30.01.2015
 *      Author: suze
 */

#ifndef DIFFUSION_BND_COND_IMPL_H_
#define DIFFUSION_BND_COND_IMPL_H_



namespace ug{
namespace MovingInterfaceDiffusion{

////////////////////////////////////////////////////////////////////////////////
//	Constructor - set default values
////////////////////////////////////////////////////////////////////////////////
// see 'no_normal_stress_outflow.cpp':
template<typename TDomain>
ParticleBndCondDiffusion<TDomain>::
ParticleBndCondDiffusion(SmartPtr<NavierStokes::NavierStokesFV1<TDomain> > spMaster,
				SmartPtr<ParticleFlatTopHandler<dim> > localHandler)
				: MovingInterface::IInterfaceBndCond<TDomain>(spMaster->symb_fcts(), spMaster->symb_subsets(), localHandler),
				  m_spMaster(spMaster),
				  m_spInterfaceHandlerLocal(localHandler)
{
//	update assemble functions
	register_all(false);
}

template<typename TDomain>
void ParticleBndCondDiffusion<TDomain>::
prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid, bool bStaticRoid)
{
// write delta t to data
	number dt = 0.13;
	if(this->is_time_dependent())
	{
		UG_LOG("ParticleBndCondDiffusion::prepare_setting(): vorher: m_dt = " << get_time_step() << "\n");

		set_time_dependent(true);

	//	get and check current and old solution
		const LocalVectorTimeSeries* vLocSol = this->local_time_solutions();
		dt = vLocSol->time(0) - vLocSol->time(1);

		set_time_step(dt);
 	}

}

// see 'NavierStokesNoNormalStressOutflowFV1::prep_elem()'
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ParticleBndCondDiffusion<TDomain>::
prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[])
{

// 	Update Geometry for this element
	static TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);

	geo.update(elem, roid, vCornerCoords, &(this->subset_handler()));

}


////////////////////////////////////////////////////////////////////////////////
//	jacobian methods
////////////////////////////////////////////////////////////////////////////////
// removes the equations with IDs in 'vFctID' locally assembled in the corners of 'vDofID'
template<typename TDomain>
void ParticleBndCondDiffusion<TDomain>::
remove_equations(LocalMatrix& J, std::vector<size_t> vFctID, std::vector<size_t> vDofID)
{

    UG_THROW("in ParticleBndCondDiffusion::remove_equations():  vFctID.size() = " << vFctID.size() << "\n");
    
	for(size_t i = 0; i < vDofID.size(); ++i)
		for(size_t j = 0; j < vFctID.size(); ++j)
			for(size_t fct = 0; fct < J.num_all_row_fct(); ++fct)
				for(size_t dof = 0; dof < J.num_all_row_dof(fct); ++dof)
					J(vFctID[j], vDofID[i], fct, dof) = 0.0;

}




// see 'NavierStokesNoNormalStressOutflowFV1::diffusive_flux_Jac()'
template<typename TDomain>
void ParticleBndCondDiffusion<TDomain>::
diffusive_flux_Jac(const size_t ip, const interfaceBF& bf, LocalMatrix& J, const LocalVector& u)
{

	MathMatrix<dim,dim> diffFlux, tang_diffFlux;
	MathVector<dim> normalStress;

	for(size_t sh = 0; sh < bf.num_sh(); ++sh) // loop shape functions
	{
	//	1. Compute the total flux
	//	- add \nabla u
		MatSet (diffFlux, 0);
		MatDiagSet (diffFlux, VecDot (bf.global_grad(sh), bf.normal()));


	//	- add (\nabla u)^T
		if( !m_spMaster->laplace())
			for (size_t d1 = 0; d1 < (size_t)dim; ++d1)
				for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
					diffFlux(d1,d2) += bf.global_grad(sh)[d1] * bf.normal()[d2];

  	//	2. Scale by viscosity
		diffFlux *= get_kinVisc_fluid();

	//	3. Change sign, since force acts onto particle in inverse direction:
		//diffFlux *= -1.0;

	//	4. Add flux to local Jacobian
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			for(size_t d2 = 0; d2 < (size_t)dim; ++d2){
				J(d1, bf.node_id(), d2, sh) += diffFlux (d1, d2);
			// 	if ( sh < 2 ) UG_LOG("bf.node_id(): " << bf.node_id() << "sh: " << sh << ": diff_flux added: " << diffFlux (d1, d2) << "\n");
			}

	//	5. Add pressure term to local Jacobian
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			J(d1, bf.node_id(), _P_, sh) -= bf.shape(sh) * bf.normal()[d1];

	}

}



void copy_and_reset(LocalMatrix& cpJ, LocalMatrix& J)
{

    UG_THROW("oha!!...in ParticleBndCondDiffusion: copy_and_reset()\n");
    
	for(size_t fct1 = 0; fct1 < J.num_all_row_fct(); ++fct1)
		for(size_t dof1 = 0; dof1 < J.num_all_row_dof(fct1); ++dof1)
			for(size_t fct2 = 0; fct2 < J.num_all_row_fct(); ++fct2)
				for(size_t dof2 = 0; dof2 < J.num_all_row_dof(fct2); ++dof2)
				{
					cpJ(fct1, dof1, fct2, dof2) = J(fct1, dof1, fct2, dof2);
					J(fct1, dof1, fct2, dof2) = 0.0;
				}
}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ParticleBndCondDiffusion<TDomain>::
add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	ElementModus elemModus = m_spInterfaceHandlerLocal->get_element_modus(elem); // computed via 'compute_element_modus()' during 'update_marker()'

    const LocalIndices& rowInd = J.get_row_indices();

	////////////////////////////////////////////////////////////////////////////////
	//  Remove local impuls equations for interface-corners
	if ( ! m_spInterfaceHandlerLocal->StdFV_assembling() ) // only remove, if NOT StdFVAssembling
	{

		if ( elemModus != INSIDE_DOM )
		{
		// all impulse equations will be removed
			std::vector<size_t> vFctID(dim);
			for(size_t i = 0; i < dim; ++i) vFctID[i] = i;

		// equations for interface nodes will be removed
			std::vector<size_t> vDoFID; vDoFID.clear();
			for(size_t i = 0; i < m_spInterfaceHandlerLocal->m_vInterfaceID.size(); ++i)
			{
				size_t interfaceID = m_spInterfaceHandlerLocal->m_vInterfaceID[i];
				size_t indexToRemove = interfaceID; //m_spInterfaceHandlerLocal->m_vOriginalCornerID[interfaceID];
				vDoFID.push_back(indexToRemove);
			}

			remove_equations(J, vFctID, vDoFID);
		}

	}
	

////////////////////////////////////////////////////////////////////////////////
// see 'ParticleFlatTop::add_jac_A_elem_interface()'
////////////////////////////////////////////////////////////////////////////////

// 	Loop the boundary faces for new impuls equations
// --> IFF INSIDE_DOM: vBF.size() = 0 ;-)
	SmartPtr<CplUserData<number, dim> > fluidDensity = m_spMaster->density();
	std::vector<interfaceBF>& vBF = m_spInterfaceHandlerLocal->get_boundary_faces();

// initialize data (written during call of 'diffusive_flux_Jac()':
	LocalIndices ind = u.get_indices();
	rotJ_ind.resize(ind); rotJ_ind = 0.0;
	rotJ_rot.resize(ind); rotJ_rot = 0.0;

 	for(size_t ip = 0; ip < vBF.size(); ++ip)
	{

		interfaceBF bf = vBF[ip];

		if ( fluidDensity->value(0, ip) != fluidDensity->value(1, ip) )
			UG_THROW("ParticleBndCondDiffusion::add_jac_M_elem(): density different for series 0 and 1: "
					<< fluidDensity->value(0, ip) << " != " << fluidDensity->value(1, ip) << "\n");

		number importDensity = fluidDensity->value(0, ip);

	// The momentum equation:
		if ( ! m_spInterfaceHandlerLocal->StdFV_assembling() ) // only remove, if NOT StdFVAssembling
			diffusive_flux_Jac(ip, bf, J, u);
 
	// scale with deltaT ( = 1.0 for non-time-dependent)
		// buffJ *= deltaT;
		// ToDo --> NOT necessary, since add_def_A_elem() will be multiplied by 'dt'
		// during elem_dis_assemble_util.h ???

	//	The continuity equation
		if ( 1 ) //! m_spInterfaceHandlerLocal->StdFV_assembling() ) // only remove, if NOT StdFVAssembling
		{
			for(size_t sh = 0; sh < bf.num_sh(); ++sh) // loop shape functions
				for (size_t d2 = 0; d2 < (size_t)dim; ++d2)
					J(_P_, bf.node_id(), d2, sh) += bf.shape(sh) * bf.normal()[d2]
					              * importDensity;
		}

	} // end vBF

  // copy rotJ_ind/rotJ_rot to m_spInterfaceHandlerLocal data for access from mapper
 	if ( vBF.size() > 0 )
 		copy_local_couplings_jac();


}


////////////////////////////////////////////////////////////////////////////////
//	defect methods
////////////////////////////////////////////////////////////////////////////////

// removes the defects with IDs in 'vFctID' locally assembled in the corners of 'vDofID'
template<typename TDomain>
void ParticleBndCondDiffusion<TDomain>::
remove_equations(LocalVector& d, std::vector<size_t> vFctID, std::vector<size_t> vDofID)
{
 	for(size_t i = 0; i < vDofID.size(); ++i)
		for(size_t j = 0; j < vFctID.size(); ++j)
  			d(vFctID[j], vDofID[i]) = 0.0;

}



// see 'NavierStokesNoNormalStressOutflowFV1::diffusive_flux_defect()'
template<typename TDomain>
void ParticleBndCondDiffusion<TDomain>::
diffusive_flux_defect(const size_t ip, const interfaceBF& bf, LocalVector& d, const LocalVector& u)
{
  	MathMatrix<dim, dim> gradVel;
	MathVector<dim> diffFlux;

// 	1. Get the gradient of the velocity at ip
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		for(size_t d2 = 0; d2 < (size_t)dim; ++d2)
		{
		//	sum up contributions of each shape
			gradVel(d1, d2) = 0.0;
			for(size_t sh = 0; sh < bf.num_sh(); ++sh)
				gradVel(d1, d2) += bf.global_grad(sh)[d2] * u(d1, sh);
		}

//	2. Compute the total flux

//	- add (\nabla u) \cdot \vec{n}
	MatVecMult(diffFlux, gradVel, bf.normal());

//	- add (\nabla u)^T \cdot \vec{n}
  	if( !m_spMaster->laplace())
		TransposedMatVecMultAdd(diffFlux, gradVel, bf.normal());

//	3. Scale by viscosity
	diffFlux *= get_kinVisc_fluid();

//	4. Change sign, since force acts onto particle in inverse direction:
	//diffFlux *= -1.0;

//	5. Add flux to local defect
	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
 		d(d1, bf.node_id()) += diffFlux[d1];


//	6. Add pressure term to local defect
	number pressure = 0.0;
	for(size_t sh = 0; sh < bf.num_sh(); ++sh)
 			pressure += bf.shape(sh) * u(_P_, sh);

	for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
		d(d1, bf.node_id()) -= pressure * bf.normal()[d1];

}

template<typename TDomain>
void ParticleBndCondDiffusion<TDomain>::
add_def_A_elem_Quadri_for2(LocalVector& locD, const LocalVector locU)
{
//  Loop the boundary faces to assemble impulse equations
	SmartPtr<CplUserData<number, dim> > fluidDensity = m_spMaster->density();
	std::vector<interfaceBF>& vBF = m_spInterfaceHandlerLocal->get_boundary_faces();

	if ( vBF.size() != 4 )
		UG_THROW("in 'ParticleBndCondDiffusion::add_def_A_elem_Quadri_for2(): vBF.size() should be 4 but is " << vBF.size() << "\n");

	UG_LOG("vBF.size(): " << vBF.size() << "\n");

// get density
	number importDensity = fluidDensity->value(0, 0);
	if ( fluidDensity->value(0, 0) != fluidDensity->value(1, 0) )
		UG_THROW("ParticleBndCondDiffusion::add_def_A_elem_Quadri_for2(): density different for series 0 and 1: "
				<< fluidDensity->value(0, 0) << " != " << fluidDensity->value(1, 0) << "\n");
	if ( importDensity != get_density_fluid() )
		UG_THROW("ParticleBndCondDiffusion::add_def_A_elem_Quadri_for2(): importDensity = " << importDensity << " = "
				"fluidDensity = " << get_density_fluid() << "\n");

// loop all boundary faces
	for(size_t ip = 0; ip < vBF.size(); ++ip)
	{
		UG_LOG("---- ip = " << ip << "\n");

		interfaceBF bf = vBF[ip];

	// Compute Velocity at ip
		MathVector<dim> stdVel(0.0);
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			for(size_t sh = 0; sh < bf.num_sh(); ++sh)
				stdVel[d1] += locU(d1, sh) * bf.shape(sh);

	// Momentum equation:
		if ( ! m_spInterfaceHandlerLocal->StdFV_assembling() ) // only remove, if NOT StdFVAssembling
			diffusive_flux_defect(ip, bf, locD, locU);

	// Continuity equation:
		if ( ! m_spInterfaceHandlerLocal->StdFV_assembling() ) // only remove, if NOT StdFVAssembling
			locD(_P_, bf.node_id()) += VecDot(stdVel, bf.normal()) * importDensity;

	} // end vBF-loop

}

template<typename TDomain>
void ParticleBndCondDiffusion<TDomain>::
write_QuadriSol(const LocalVector origU)
{
	UG_LOG("start write_QuadriSol\n");

// initialize data
	LocalIndices ind = origU.get_indices();
 	LocalVector quadriU;

	UG_LOG("1 start write_QuadriSol\n");

	if ( m_spInterfaceHandlerLocal->interface_id_all().size() == 2 )
	{
	  	for(size_t fct = 0; fct < ind.num_fct(); ++fct)
			ind.resize_dof(fct, 5);

	  	if ( ind.num_dof(0) != 5 )
	  		UG_THROW("hmm: ind.num_dof(0) = " << ind.num_dof(0) << "\n");
	}
 	quadriU.resize(ind);

	UG_LOG("2 start write_QuadriSol\n");

 // A. remap solution (velocity AND pressure!) of inside node
 	if ( m_spInterfaceHandlerLocal->interface_id_all().size() == 2 )
 	{
		size_t copyID = m_spInterfaceHandlerLocal->m_vNOInterfaceID[0];
 		for(size_t fct = 0; fct < dim+1; ++fct)
 			{quadriU(fct,4) = origU(fct,copyID);
 			UG_LOG("origU(fct,copyID) = " << origU(fct,copyID) << "\n");}
 	}

	UG_LOG("3 start write_QuadriSol: quadriU \n" << quadriU << "\n");

 	if ( m_spInterfaceHandlerLocal->interface_id_all().size() == 2 )
 	{
	for ( size_t i = 0; i < m_spInterfaceHandlerLocal->m_vNOInterfaceID.size(); ++i )
		UG_THROW("m_vNOInterfaceID[" << i << "]: " << m_spInterfaceHandlerLocal->m_vNOInterfaceID[i] << "\n");

 	}
 	for ( size_t i = 0; i < m_spInterfaceHandlerLocal->m_vInterfaceID.size(); ++i )
		UG_LOG("m_vInterfaceID[" << i << "]: " << m_spInterfaceHandlerLocal->m_vInterfaceID[i] << "\n");
	for ( size_t i = 0; i < m_spInterfaceHandlerLocal->m_vQuadriOrigID.size(); ++i )
		UG_LOG("m_vQuadriOrigID[" << i << "]: " << m_spInterfaceHandlerLocal->m_vQuadriOrigID[i] << "\n");

	std::vector<number> testV;
	for(size_t dof = 0; dof < 4; ++dof)
		testV.push_back(1.0+0.1*dof);

 // B. remap pressure solution in all interface nodes
	for(size_t fct = 0; fct < dim+1; ++fct)
		for(size_t dof = 0; dof < 4; ++dof)
		{
			quadriU(fct,dof) = testV[m_spInterfaceHandlerLocal->m_vQuadriOrigID[dof]]; //origU(fct,m_spInterfaceHandlerLocal->m_vQuadriOrigID[dof]);
			UG_LOG("quadriU(fct,dof): " << quadriU(fct,dof) << "\n");
			UG_LOG("testV[" << m_spInterfaceHandlerLocal->m_vQuadriOrigID[dof] << "]: " << testV[m_spInterfaceHandlerLocal->m_vQuadriOrigID[dof]] << "\n");
		}

 	UG_LOG("start write_QuadriSol: quadriU = \n" << quadriU << "\n");
	UG_THROW("4 start write_QuadriSol\n");

// C. write velocities of the 2 particles:
 	MathVector<dim> transSol1 = m_spInterfaceHandlerLocal->get_transSol(0, 0);
	MathVector<dim> rotSol1 = m_spInterfaceHandlerLocal->get_rotSol(0, 0);
 	MathVector<dim> transSol2 = m_spInterfaceHandlerLocal->get_transSol(1, 0);
	MathVector<dim> rotSol2 = m_spInterfaceHandlerLocal->get_rotSol(1, 0);

	UG_LOG("transSol1 = " << transSol1 << "\n");
	UG_LOG("transSol2 = " << transSol2 << "\n");
	UG_LOG("rotSol1 = " << rotSol1 << "\n");
	UG_LOG("rotSol2 = " << rotSol2 << "\n");

// resize local data as done during 'modify_LocalSol'
//	 --> but there with num_co = 3 for the Triangle with inside node!
	for(size_t fct = 0; fct < ind.num_fct(); ++fct)
		for(size_t dof = 0; dof < 4; ++dof)
		{
			MathMatrix<dim,dim> rotationMatCo = m_spInterfaceHandlerLocal->get_rotationMat(m_spInterfaceHandlerLocal->radial_at_co(dof));
			UG_LOG("write_QuadriSol: radial(" << dof << ": " << m_spInterfaceHandlerLocal->radial_at_co(dof) << "\n");

		// write solution of particle with prtIndex = 0:
			if ( dof < 2 )
			{
				quadriU(fct,dof) = transSol1[fct];
				for ( int d = 0; d < dim; ++d )
					quadriU(fct,dof) += rotationMatCo[fct][d]*rotSol1[d];
			}
		// write solution of particle with prtIndex = 1:
			else
			{
				quadriU(fct,dof) = transSol2[fct];
				for ( int d = 0; d < dim; ++d )
					quadriU(fct,dof) += rotationMatCo[fct][d]*rotSol2[d];
			}
		}


	UG_LOG("after write_QuadriSol: quadriU = \n" << quadriU << "\n");

}

template<typename TDomain>
void ParticleBndCondDiffusion<TDomain>::
add_quadri_to_defect(LocalVector& d, const LocalVector& quadriD)
{
	const LocalIndices& ind = quadriD.get_indices();
	UG_LOG("START add_quadri_to_defect \n");

	for ( size_t i = 0; i < m_spInterfaceHandlerLocal->m_vQuadriOrigID.size(); ++i )
		UG_LOG("m_vOrig[" << i << "] = " << m_spInterfaceHandlerLocal->m_vQuadriOrigID[i] << "\n");

	UG_LOG("d: \n" << d << "\n");
	UG_LOG("quadriD: \n" << quadriD << "\n");
	for(size_t fct=0; fct < quadriD.num_all_fct(); ++fct)
	{
		UG_LOG("fct = " << fct << "quadriD.num_all_dof(fct) = " << quadriD.num_all_dof(fct) << "\n");

		for(size_t dof=0; dof < quadriD.num_all_dof(fct); ++dof)
		{
			UG_LOG("0 -> dof = " << dof << "\n");

			if ( quadriD.value(fct,dof) != quadriD.value(fct,dof))
				UG_THROW("NAN in 'add_quadri_to_defect()'!...\n");

			bool isPrtNode = m_spInterfaceHandlerLocal->lies_onInterface(dof);
		// usual assembling for fluid-dof and pressure-fct
			if ( isPrtNode )
			{
				size_t _dof = m_spInterfaceHandlerLocal->m_vQuadriOrigID[dof];
				UG_LOG("dof = " << dof << ", _dof = " << _dof << "\n");

				d.value(fct,_dof) += quadriD.value(fct,dof);
				UG_LOG("d = \n" << d << "\n");

			}

		}
	}

}



template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ParticleBndCondDiffusion<TDomain>::
add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
////////////////////////////////////////////////////////////////////////////////
//  Remove local impulse equations for interface-corners

 	ElementModus elemModus = m_spInterfaceHandlerLocal->get_element_modus(elem); // computed via 'compute_element_modus()' during 'update_marker()'
	if ( ! m_spInterfaceHandlerLocal->StdFV_assembling() ) // only remove, if NOT StdFVAssembling
	{
		if ( elemModus != INSIDE_DOM )
		{
		// all impulse equations will be removed
			std::vector<size_t> vFctID(dim);
			for(size_t i = 0; i < dim; ++i) vFctID[i] = i;

		// equations for interface nodes will be removed
			std::vector<size_t> vDoFID; vDoFID.clear();
			for(size_t i = 0; i < m_spInterfaceHandlerLocal->m_vInterfaceID.size(); ++i)
			{
				size_t interfaceID = m_spInterfaceHandlerLocal->m_vInterfaceID[i];
				size_t indexToRemove = interfaceID; //m_spInterfaceHandlerLocal->m_vOriginalCornerID[interfaceID];
				vDoFID.push_back(indexToRemove);
			}

			remove_equations(d, vFctID, vDoFID);

		}
	}

 
	if ( elemModus == CUT_BY_2_INTERFACE && !m_spInterfaceHandlerLocal->StdFV_assembling())
	{
	 	write_QuadriSol(u);

		if ( m_spInterfaceHandlerLocal->interface_id_all().size() == 2 )

	 // case == 4: all entries (also for pressure eq) will be written newly, since only boundary faces are relevant
	 //		=> during remove_equations() only velocity equations were removed!
	 	if ( m_spInterfaceHandlerLocal->interface_id_all().size() == 4 )
	 		d = 0.0;

 	// locU was written during 'ParticleMapper::modify_LocalData()':
		if ( m_spInterfaceHandlerLocal->interface_id_all().size() == 4 || m_spInterfaceHandlerLocal->interface_id_all().size() == 2 )
   			add_def_A_elem_Quadri_for2(d, u);
		else
			UG_THROW("in ParticleBndCondDiffusion::add_def_A_elem(): wrong amount of interface.size() for CUT_BY_2_INTERFACE: " << m_spInterfaceHandlerLocal->interface_id_all().size() << " ( should be 2!)\n");

   		return;
	}
 
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// initialize data (written during call of 'diffusive_flux_Jac()':
	LocalIndices ind = u.get_indices();
	rotD.resize(ind); rotD = 0.0;


	// ToDo: if ( CUT_BY_2 && m_spInterfaceHandlerLocal->m_vInterfaceID.size() == 2 )
		//	=> TRIANGLE/5-Eck! => write u( , ) neu!!

//  Loop the boundary faces to assemble impulse equations
	SmartPtr<CplUserData<number, dim> > fluidDensity = m_spMaster->density();
	std::vector<interfaceBF>& vBF = m_spInterfaceHandlerLocal->get_boundary_faces();

 
	if ( dim == 2 && vBF.size() > 2 )
		UG_THROW("add_def_A_elem(): vBF.size() is greater than 2: " << vBF.size() << "\n");

	for(size_t ip = 0; ip < vBF.size(); ++ip)
 	{
 		interfaceBF bf = vBF[ip];

	// Compute Velocity at ip
		MathVector<dim> stdVel(0.0);
		for(size_t d1 = 0; d1 < (size_t)dim; ++d1)
			for(size_t sh = 0; sh < bf.num_sh(); ++sh)
 				stdVel[d1] += u(d1, sh) * bf.shape(sh);

	// Momentum equation:
		if ( ! m_spInterfaceHandlerLocal->StdFV_assembling() ) // only remove, if NOT StdFVAssembling
			diffusive_flux_defect(ip, bf, d, u);

	// scale with deltaT ( = 1.0 for non-time-dependent)
	//	d *= deltaT;
		// ToDo --> NOT necessary, since add_def_A_elem() will be multiplied by 'dt'
		// during elem_dis_assemble_util.h ???

		if ( fluidDensity->value(0, ip) != fluidDensity->value(1, ip) )
			UG_THROW("ParticleBndCondDiffusion::add_jac_M_elem(): density different for series 0 and 1: "
					<< fluidDensity->value(0, ip) << " != " << fluidDensity->value(1, ip) << "\n");
		number importDensity = fluidDensity->value(0, ip);
		if ( importDensity != get_density_fluid() )
			UG_THROW("ParticleBndCondDiffusion::add_def_A_elem(): importDensity = " << importDensity << " != "
					"fluidDensity = " << get_density_fluid() << "\n");

 
	// Continuity equation:
		if ( 1 ) //! m_spInterfaceHandlerLocal->StdFV_assembling() ) // only remove, if NOT StdFVAssembling
			d(_P_, bf.node_id()) += VecDot(stdVel, bf.normal()) * importDensity;
		// ToDo m_massDefect += VecDot(stdVel, bf.normal());
	}

/*	if ( elemModus == CUT_BY_INTERFACE )
		if ( m_spInterfaceHandlerLocal->m_vInterfaceID.size() == 2 )
			UG_THROW("vBF.size() = " << vBF.size() << "\n");
*/
 // copy rotJ_ind/rotJ_rot to m_spInterfaceHandlerLocal data for access from mapper
 
 	copy_local_couplings_def();


////////////////////////////////////////////////////////////////////////////////
// instead of 'set_gravity()' during 'add_local_vec_to_global()'
////////////////////////////////////////////////////////////////////////////////

 
	//ToDo...


}


// instead of calling 'set_mass_and_inertia()' during 'add_local_mat_to_global_interface()'
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ParticleBndCondDiffusion<TDomain>::
add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
 	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

// A. do nothing for inside elements
	if ( modus == INSIDE_DOM )
		return;

// B. for outside elements: geo.num_scv() = 0
//		=> nothing added during NavierStokesFV1::add_jac_M_elem
//		 	=> nothing will be subtracted here :)


///////////////////////////////////////////////////////////////////////////////
// substract part added by NavierStokesFV1::add_jac_M_elem():
// 		(instead of setting scv.volume() to zero)
////////////////////////////////////////////////////////////////////////////////

// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
 	static TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);

	SmartPtr<CplUserData<number, dim> > fluidDensity = m_spMaster->density();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

		if ( !m_spInterfaceHandlerLocal->lies_onInterface(sh) )
			continue;

		if ( fluidDensity->value(0, ip) != fluidDensity->value(1, ip) )
			UG_THROW("ParticleBndCondDiffusion::add_jac_M_elem(): density different for series 0 and 1: "
					<< fluidDensity->value(0, ip) << " != " << fluidDensity->value(1, ip) << "\n");

		number importDensity = fluidDensity->value(0, ip);

	// 	loop velocity components
		for(int d1 = 0; d1 < dim; ++d1)
		{
		// 	Add to local matrix
			J(d1, sh, d1, sh) -= scv.volume() * importDensity;
		}
	}

 	return;


 	const int prtIndex = get_prtIndex();
	const number prtDensity = get_density(prtIndex);

// if INSIDE_DOM: do nothing!
	if ( modus != INSIDE_DOM )
	{
	// get data
		const ReferenceObjectID roid = elem->reference_object_id();
		const DimReferenceElement<dim>& rRefElem
			= ReferenceElementProvider::get<dim>(roid);

	// OUTSIDE_DOM
		if ( modus == OUTSIDE_DOM )
		{
		// loop corners of reference element
			for(size_t sh = 0; sh < rRefElem.num(0); ++sh)
			{
			// 	loop velocity components
				for(int d1 = 0; d1 < dim; ++d1)
				{
				// 	Add to local matrix
					J(d1, sh, d1, sh) += geo.volume_fem_elem() * prtDensity;
				// rescale volume fraction
					J(d1, sh, d1, sh) *= 1.0/rRefElem.num(0);
				}
			}
		}
	// CUT_BY_INTERFACE
		else if ( modus == CUT_BY_INTERFACE )
		{
		// loop corners of reference element
			for(size_t sh = 0; sh < rRefElem.num(0); ++sh)
			{
			// 	loop velocity components
				for(int d1 = 0; d1 < dim; ++d1)
				{
				// 	Add to local matrix
					J(d1, sh, d1, sh) += geo.volume_fem_elem() * prtDensity;
					UG_THROW("geo.volume_fem_elem() = " << geo.volume_fem_elem() << "\n");
				// rescale volume fraction
					J(d1, sh, d1, sh) *= 1.0/rRefElem.num(0);
				}
			}
		}
		else
			UG_THROW("ParticleBndCondDiffusion::add_jac_M_elem()...");

	} // end 'if ( !is_inside_elem() )'


}


template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ParticleBndCondDiffusion<TDomain>::
add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
 	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

// A. do nothing for inside elements
	if ( modus == INSIDE_DOM )
		return;

// B. for outside elements: geo.num_scv() = 0
//		=> nothing added during NavierStokesFV1::add_jac_M_elem
//		 	=> nothing will be subtracted here :)


///////////////////////////////////////////////////////////////////////////////
// substract part added by NavierStokesFV1::add_jac_M_elem():
// 		(instead of setting scv.volume() to zero)
////////////////////////////////////////////////////////////////////////////////

// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
 	static TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);

	SmartPtr<CplUserData<number, dim> > fluidDensity = m_spMaster->density();

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

		if ( !m_spInterfaceHandlerLocal->lies_onInterface(sh) )
			continue;

		if ( fluidDensity->value(0, ip) != fluidDensity->value(1, ip) )
			UG_THROW("ParticleBndCondDiffusion::add_jac_M_elem(): density different for series 0 and 1: "
					<< fluidDensity->value(0, ip) << " != " << fluidDensity->value(1, ip) << "\n");

		number importDensity = fluidDensity->value(0, ip);

	// 	loop velocity components
		for(int d1 = 0; d1 < dim; ++d1)
		{
		// 	Add to local matrix
			d(d1, sh) -= u(d1, sh) * scv.volume() * importDensity;
		}
	}


}


// here gravity force -> independent of velocity solution!
template<typename TDomain>
template<typename TElem, typename TFVGeom>
void ParticleBndCondDiffusion<TDomain>::
add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[])
{
 	return;


 	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

// A. do nothing for inside elements
	if ( modus == INSIDE_DOM )
		return;

// B. for outside elements: geo.num_scv() = 0
//		=> nothing added during NavierStokesFV1::add_jac_M_elem
//		 	=> nothing will be subtracted here :)


///////////////////////////////////////////////////////////////////////////////
// substract part added by NavierStokesFV1::add_rgh_elem():
// 		(instead of setting scv.volume() to zero)
////////////////////////////////////////////////////////////////////////////////

// 	Only first order implementation
	UG_ASSERT((TFVGeom::order == 1), "Only first order implemented.");

// 	get finite volume geometry
 	static TFVGeom& geo = GeomProvider<TFVGeom>::get(LFEID(LFEID::LAGRANGE, dim, 1), 1);

// 	loop Sub Control Volumes (SCV)
	for(size_t ip = 0; ip < geo.num_scv(); ++ip)
	{
	// 	get current SCV
		const typename TFVGeom::SCV& scv = geo.scv(ip);

	// 	get associated node
		const int sh = scv.node_id();

		if ( !m_spInterfaceHandlerLocal->lies_onInterface(sh) )
			continue;

	// 	Add to local matrix
		d(0, sh) -= scv.volume() * (-9.81);

	}


}



////////////////////////////////////////////////////////////////////////////////
//	register assemble functions
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template<>
void ParticleBndCondDiffusion<Domain1d>::
register_all(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<RegularEdge, MovingParticle::DimFV1FTGeometry<dim, MovingInterface::FlatTopHandler<dim> > >();

//		register_func<RegularEdge, FV1Geometry<RegularEdge, dim> >();
	}
	else
	{
		UG_THROW("ParticleBndCondDiffusion: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_2
template<>
void ParticleBndCondDiffusion<Domain2d>::
register_all(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Triangle,  MovingParticle::DimFV1FTGeometry<dim, MovingInterface::FlatTopHandler<dim> > >();
/*
		register_func<Triangle, FV1Geometry<Triangle, dim> >();
		register_func<Quadrilateral, FV1Geometry<Quadrilateral, dim> >();
		*/
	}
	else
	{
		UG_THROW("ParticleBndCondDiffusion: Hanging Nodes not implemented.")
	}
}
#endif

#ifdef UG_DIM_3
template<>
void ParticleBndCondDiffusion<Domain3d>::
register_all(bool bHang)
{
//	switch assemble functions
	if(!bHang)
	{
		register_func<Tetrahedron,  MovingParticle::DimFV1FTGeometry<dim, MovingInterface::FlatTopHandler<dim> > >();
/*
		register_func<Tetrahedron, FV1Geometry<Tetrahedron, dim> >();
		register_func<Prism, FV1Geometry<Prism, dim> >();
		register_func<Pyramid, FV1Geometry<Pyramid, dim> >();
		register_func<Hexahedron, FV1Geometry<Hexahedron, dim> >();
		*/
	}
	else
	{
		UG_THROW("ParticleBndCondDiffusion: Hanging Nodes not implemented.")
	}
}
#endif

template<typename TDomain>
template<typename TElem, typename TFVGeom>
void
ParticleBndCondDiffusion<TDomain>::
register_func()
{
	ReferenceObjectID id = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
	typedef ParticleBndCondDiffusion<TDomain> T;

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
template class ParticleBndCondDiffusion<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ParticleBndCondDiffusion<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ParticleBndCondDiffusion<Domain3d>;
#endif


} // end namespace MovingInterfaceDiffusion
} // end namespace ug


#endif /* DIFFUSION_BND_COND_IMPL_H_ */
