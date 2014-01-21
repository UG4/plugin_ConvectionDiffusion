/*
 * convection_diffusion_base.cpp
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#include "convection_diffusion_base.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/user_data/data_export.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{
namespace ConvectionDiffusionPlugin{

////////////////////////////////////////////////////////////////////////////////
//	user data
////////////////////////////////////////////////////////////////////////////////

//////// Diffusion

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> > user)
{
	m_imDiffusion.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_diffusion(number val)
{
	if(val == 0.0) set_diffusion(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> >());
	else set_diffusion(make_sp(new ConstUserMatrix<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_diffusion(const char* fctName)
{
	set_diffusion(LuaUserDataFactory<MathMatrix<dim,dim>, dim>::create(fctName));
}
#endif

//////// Velocity

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imVelocity.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_velocity(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_velocity(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_velocity(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_velocity(const char* fctName)
{
	set_velocity(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
#endif

//////// Flux

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imFlux.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_flux(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_flux(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_flux(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_flux(const char* fctName)
{
	set_flux(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
#endif

//////// Reaction Rate

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReactionRate.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate(number val)
{
	if(val == 0.0) set_reaction_rate(SmartPtr<CplUserData<number, dim> >());
	else set_reaction_rate(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate(const char* fctName)
{
	set_reaction_rate(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Reaction Rate Explicit

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate_explicit(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReactionRateExpl.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate_explicit(number val)
{
	if(val == 0.0) set_reaction_rate_explicit(SmartPtr<CplUserData<number, dim> >());
	else set_reaction_rate_explicit(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_rate_explicit(const char* fctName)
{
	set_reaction_rate_explicit(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Reaction

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReaction.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction(number val)
{
	if(val == 0.0) set_reaction(SmartPtr<CplUserData<number, dim> >());
	else set_reaction(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction(const char* fctName)
{
	set_reaction(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Reaction Explicit

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_explicit(SmartPtr<CplUserData<number, dim> > user)
{
	m_imReactionExpl.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_explicit(number val)
{
	if(val == 0.0) set_reaction_explicit(SmartPtr<CplUserData<number, dim> >());
	else set_reaction_explicit(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_reaction_explicit(const char* fctName)
{
	set_reaction_explicit(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif


//////// Source

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source(SmartPtr<CplUserData<number, dim> > user)
{
	m_imSource.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source(number val)
{
	if(val == 0.0) set_source(SmartPtr<CplUserData<number, dim> >());
	else set_source(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source(const char* fctName)
{
	set_source(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Source explicit

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source_explicit(SmartPtr<CplUserData<number, dim> > user)
{
	m_imSourceExpl.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source_explicit(number val)
{
	if(val == 0.0) set_source_explicit(SmartPtr<CplUserData<number, dim> >());
	else set_source_explicit(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_source_explicit(const char* fctName)
{
	set_source_explicit(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif


//////// Vector Source

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_vector_source(SmartPtr<CplUserData<MathVector<dim>, dim> > user)
{
	m_imVectorSource.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_vector_source(const std::vector<number>& vVel)
{
	bool bZero = true;
	for(size_t i = 0; i < vVel.size(); ++i){
		if(vVel[i] != 0.0) bZero = false;
	}

	if(bZero) set_vector_source(SmartPtr<CplUserData<MathVector<dim>, dim> >());
	else set_vector_source(SmartPtr<ConstUserVector<dim> >(new ConstUserVector<dim>(vVel)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::set_vector_source(const char* fctName)
{
	set_vector_source(LuaUserDataFactory<MathVector<dim>,dim>::create(fctName));
}
#endif

//////// Mass Scale

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass_scale(SmartPtr<CplUserData<number, dim> > user)
{
	m_imMassScale.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass_scale(number val)
{
	if(val == 0.0) set_mass_scale(SmartPtr<CplUserData<number, dim> >());
	else set_mass_scale(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass_scale(const char* fctName)
{
	set_mass_scale(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

//////// Mass

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass(SmartPtr<CplUserData<number, dim> > user)
{
	m_imMass.set_data(user);
}

template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass(number val)
{
	if(val == 0.0) set_mass(SmartPtr<CplUserData<number, dim> >());
	else set_mass(make_sp(new ConstUserNumber<dim>(val)));
}

#ifdef UG_FOR_LUA
template<typename TDomain>
void ConvectionDiffusionBase<TDomain>::
set_mass(const char* fctName)
{
	set_mass(LuaUserDataFactory<number,dim>::create(fctName));
}
#endif

////////////////////////////////////////////////////////////////////////////////
//	Exports
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain>
typename ConvectionDiffusionBase<TDomain>::NumberExport
ConvectionDiffusionBase<TDomain>::
value() {return m_exValue;}


template <typename TDomain>
typename ConvectionDiffusionBase<TDomain>::GradExport
ConvectionDiffusionBase<TDomain>::
gradient() {return m_exGrad;}

////////////////////////////////////////////////////////////////////////////////
//	Constructor
////////////////////////////////////////////////////////////////////////////////

template<typename TDomain>
ConvectionDiffusionBase<TDomain>::
ConvectionDiffusionBase(const char* functions, const char* subsets)
 : IElemDisc<TDomain>(functions,subsets),
   m_exValue(new DataExport<number, dim>(functions)),
   m_exGrad(new DataExport<MathVector<dim>, dim>(functions))
{
//	check number of functions
	if(this->num_fct() != 1)
		UG_THROW("Wrong number of functions: The ElemDisc 'ConvectionDiffusion'"
					   " needs exactly "<<1<<" symbolic function.");

//	register imports
	this->register_import(m_imDiffusion);
	this->register_import(m_imVelocity);
	this->register_import(m_imFlux);
	this->register_import(m_imReactionRate);
	this->register_import(m_imReaction);
	this->register_import(m_imReactionRateExpl);
	this->register_import(m_imReactionExpl);
	this->register_import(m_imSourceExpl);
	this->register_import(m_imSource);
	this->register_import(m_imVectorSource);
	this->register_import(m_imMassScale);
	this->register_import(m_imMass);

	m_imMassScale.set_mass_part();
	m_imMass.set_mass_part();
	m_imSource.set_rhs_part();
	m_imVectorSource.set_rhs_part();
	m_imSourceExpl.set_expl_part();
	m_imReactionExpl.set_expl_part();
	m_imReactionRateExpl.set_expl_part();

//	default value for mass scale
	set_mass_scale(1.0);
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template class ConvectionDiffusionBase<Domain1d>;
#endif
#ifdef UG_DIM_2
template class ConvectionDiffusionBase<Domain2d>;
#endif
#ifdef UG_DIM_3
template class ConvectionDiffusionBase<Domain3d>;
#endif

} // end namespace ConvectionDiffusionPlugin
} // namespace ug
