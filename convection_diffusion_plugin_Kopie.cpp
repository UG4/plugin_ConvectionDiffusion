/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

/**
 * File for registration of ConvectionDiffusion routines.
 */

#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "convection_diffusion_base.h"
#include "fv1/convection_diffusion_fv1.h"
#include "fe/convection_diffusion_fe.h"
#include "fvcr/convection_diffusion_fvcr.h"
#include "fv/convection_diffusion_fv.h"

#include "lib_disc/spatial_disc/elem_disc/sss.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace ConvectionDiffusionPlugin{

/** 
 *  \defgroup convection_diffusion Convection Diffusion
 *  \ingroup plugins
 *  This plugin provides the discretization of convection and diffusion problems.
 *  \{
 */

/**
 * Class exporting the functionality of the plugin. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain dependent parts
 * of the plugin. All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	static const int dim = TDomain::dim;
	string suffix = GetDomainSuffix<TDomain>();
	string tag = GetDomainTag<TDomain>();

//	Convection Diffusion Base
	{
		typedef ConvectionDiffusionBase<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("ConvectionDiffusionBase").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_method("set_diffusion", static_cast<void (T::*)(SmartPtr<CplUserData<MathMatrix<dim, dim>, dim> >)>(&T::set_diffusion), "", "Diffusion")
			.add_method("set_diffusion", static_cast<void (T::*)(number)>(&T::set_diffusion), "", "Diagonal Diffusion")
#ifdef UG_FOR_LUA
			.add_method("set_diffusion", static_cast<void (T::*)(const char*)>(&T::set_diffusion), "", "Diffusion")
			.add_method("set_diffusion", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_diffusion), "", "Diffusion")
#endif

			.add_method("set_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(const std::vector<number>&)>(&T::set_velocity), "", "Velocity Field")
#ifdef UG_FOR_LUA
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_velocity), "", "Velocity Field")
#endif

			.add_method("set_flux", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_flux), "", "Flux")
			.add_method("set_flux", static_cast<void (T::*)(const std::vector<number>&)>(&T::set_flux), "", "Flux")
#ifdef UG_FOR_LUA
			.add_method("set_flux", static_cast<void (T::*)(const char*)>(&T::set_flux), "", "Flux")
			.add_method("set_flux", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_flux), "", "Flux")
#endif

			.add_method("set_reaction_rate", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_reaction_rate), "", "Reaction Rate")
			.add_method("set_reaction_rate", static_cast<void (T::*)(number)>(&T::set_reaction_rate), "", "Reaction Rate")
#ifdef UG_FOR_LUA
			.add_method("set_reaction_rate", static_cast<void (T::*)(const char*)>(&T::set_reaction_rate), "", "Reaction Rate")
			.add_method("set_reaction_rate", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_reaction_rate), "", "Reaction Rate")
#endif

			.add_method("set_reaction", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_reaction), "", "Reaction")
			.add_method("set_reaction", static_cast<void (T::*)(number)>(&T::set_reaction), "", "Reaction")
#ifdef UG_FOR_LUA
			.add_method("set_reaction", static_cast<void (T::*)(const char*)>(&T::set_reaction), "", "Reaction")
			.add_method("set_reaction", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_reaction), "", "Reaction")
#endif

			.add_method("set_reaction_rate_explicit", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_reaction_rate_explicit), "", "Reaction Rate Explicit")
			.add_method("set_reaction_rate_explicit", static_cast<void (T::*)(number)>(&T::set_reaction_rate_explicit), "", "Reaction Rate Explicit")
#ifdef UG_FOR_LUA
			.add_method("set_reaction_rate_explicit", static_cast<void (T::*)(const char*)>(&T::set_reaction_rate_explicit), "", "Reaction Rate Explicit")
#endif

			.add_method("set_reaction_explicit", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_reaction_explicit), "", "Reaction Explicit")
			.add_method("set_reaction_explicit", static_cast<void (T::*)(number)>(&T::set_reaction_explicit), "", "Reaction Explicit")
#ifdef UG_FOR_LUA
			.add_method("set_reaction_explicit", static_cast<void (T::*)(const char*)>(&T::set_reaction_explicit), "", "Reaction Explicit")
#endif

			.add_method("set_source_explicit", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_source_explicit), "", "Source Explicit")
			.add_method("set_source_explicit", static_cast<void (T::*)(number)>(&T::set_source_explicit), "", "Source Explicit")
			#ifdef UG_FOR_LUA
			.add_method("set_source_explicit", static_cast<void (T::*)(const char*)>(&T::set_source_explicit), "", "Source Explicit")
			#endif

			.add_method("set_source", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_source), "", "Source")
#endif

			.add_method("set_vector_source", static_cast<void (T::*)(SmartPtr<CplUserData<MathVector<dim>, dim> >)>(&T::set_vector_source), "", "Vector Source")
#ifdef UG_FOR_LUA
			.add_method("set_vector_source", static_cast<void (T::*)(const char*)>(&T::set_vector_source), "", "Vector Source")
			.add_method("set_vector_source", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_vector_source), "", "Vector Source")
#endif

			.add_method("set_mass_scale", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_mass_scale), "", "Mass Scale")
			.add_method("set_mass_scale", static_cast<void (T::*)(number)>(&T::set_mass_scale), "", "Mass Scale")
#ifdef UG_FOR_LUA
			.add_method("set_mass_scale", static_cast<void (T::*)(const char*)>(&T::set_mass_scale), "", "Mass Scale")
			.add_method("set_mass_scale", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_mass_scale), "", "Mass Scale")
#endif

			.add_method("set_mass", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_mass), "", "Mass")
			.add_method("set_mass", static_cast<void (T::*)(number)>(&T::set_mass), "", "Mass")
#ifdef UG_FOR_LUA
			.add_method("set_mass", static_cast<void (T::*)(const char*)>(&T::set_mass), "", "Mass")
			.add_method("set_mass", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_mass), "", "Mass")
#endif

			.add_method("value", &T::value)
			.add_method("gradient", &T::gradient);
		reg.add_class_to_group(name, "ConvectionDiffusionBase", tag);
	}

//	Convection Diffusion FV1
	{
		typedef ConvectionDiffusionFV1<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_upwind", &T::set_upwind)
			.add_method("set_singular_sources_and_sinks", &T::set_sss, "", "Singular Sources and Sinks")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFV1", tag);
	}

//	Convection Diffusion FE
	{
		typedef ConvectionDiffusionFE<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFE").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_quad_order", &T::set_quad_order)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFE", tag);
	}

//	Convection Diffusion FVCR
	{
		typedef ConvectionDiffusionFVCR<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFVCR").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_upwind", &T::set_upwind)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFVCR", tag);
	}

//	Convection Diffusion FV
	{
		typedef ConvectionDiffusionFV<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFV").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_quad_order", &T::set_quad_order)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFV", tag);
	}
}

template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	//	singular sources and sinks 
	{
		typedef SingularSourcesAndSinks<dim, 1> T;
		string name = string("CdSingularSourcesAndSinks").append(dimSuffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("addps", static_cast<void (T::*)(const std::vector<number>&, const std::vector<number>&)>(&T::addps))
			.add_method("addls", static_cast<void (T::*)(const std::vector<number>&, const std::vector<number>&, const std::vector<number>&)>(&T::addls))
#ifdef UG_FOR_LUA
			.add_method("addps", static_cast<void (T::*)(const std::vector<number>&, LuaFunctionHandle)>(&T::addps))
			.add_method("addls", static_cast<void (T::*)(const std::vector<number>&, const std::vector<number>&, LuaFunctionHandle)>(&T::addls))
#endif
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CdSingularSourcesAndSinks", dimTag);
	}
}

}; // end Functionality

// end group convection_diffusion
/// \}

} // end namespace ConvectionDiffusionPlugin


/**
 * This function is called when the plugin is loaded.
 */
extern "C" void
InitUGPlugin_ConvectionDiffusion(Registry* reg, string grp)
{
	grp.append("/SpatialDisc/ElemDisc");
	typedef ConvectionDiffusionPlugin::Functionality Functionality;

	try{
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
