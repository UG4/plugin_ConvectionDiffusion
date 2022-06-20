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
#include "bridge/util_domain_algebra_dependent.h"
#include "convection_diffusion_base.h"
#include "convection_diffusion_sss.h"
#include "fv1/convection_diffusion_fv1.h"
#include "fv1_cutElem/convection_diffusion_fv1_cutElem.h"
#include "fe/convection_diffusion_fe.h"
#include "fe/convection_diffusion_stab_fe.h"
#include "fvcr/convection_diffusion_fvcr.h"
#include "fv/convection_diffusion_fv.h"
#include "fv1_cutElem/diffusion_interface/diffusion_interface.h"
#include "lib_disc/spatial_disc/elem_disc/sss.h"
#include "fractfv1/convection_diffusion_fractfv1.h"

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
		  /*
			.add_method("set_partial_velocity", &T::set_partial_velocity)
			.add_method("set_partial_flux", &T::set_partial_flux)
			.add_method("set_partial_mass", &T::set_partial_mass);
		  */
		reg.add_class_to_group(name, "ConvectionDiffusionBase", tag);
	}

//	Convection Diffusion FV1
	{
		typedef ConvectionDiffusionFV1<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_condensed_FV", &T::set_condensed_FV, "", "[De-]Activates the condensed FV scvf ip's")
			.add_method("set_upwind", &T::set_upwind, "", "Sets the upwind type for the convective terms")
			.add_method("set_singular_sources_and_sinks", &T::set_sss_manager, "", "Sets the singular sources and sinks manager")
			.add_method("singular_sources_and_sinks", &T::sss_manager, "", "Returns the singular sources and sinks manager")
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFV1", tag);
	}
    
//	Convection Diffusion FV1 immersed Boundary
    {
        typedef ConvectionDiffusionFV1_cutElem<TDomain> T;
        typedef ConvectionDiffusionBase<TDomain> TBase;
        string name = string("ConvectionDiffusionFV1_cutElem").append(suffix);
        reg.add_class_<T, TBase >(name, grp)
        .template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
        .add_method("set_upwind", &T::set_upwind)
        .add_method("set_testCase", &T::set_testCase)
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "ConvectionDiffusionFV1_cutElem", tag);
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

//	Convection Diffusion (FE) stabilization
	{
		typedef ConvectionDiffusionStabFE<TDomain> T;
		typedef IElemDisc<TDomain> TBase;
		string name = string("ConvectionDiffusionStabFE").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.template add_constructor<void (*)(const char*,const char*,number)>("Function(s)#Subset(s)#stabilization")
			.template add_constructor<void (*)(const char*,const char*,number,number)>("Function(s)#Subset(s)#stabilization")
			.add_method("set_quad_order", &T::set_quad_order)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionStabFE", tag);
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

    /**
     * Function called for the registration of Domain and Algebra dependent parts
     * of the plugin. All Functions and Classes depending on both Domain and Algebra
     * are to be placed here when registering. The method is called for all
     * available Domain and Algebra types, based on the current build options.
     *
     * @param reg				registry
     * @param parentGroup		group for sorting of functionality
     */
    template <typename TDomain, typename TAlgebra>
    static void DomainAlgebra(Registry& reg, string grp)
    {
        string suffix = GetDomainAlgebraSuffix<TDomain,TAlgebra>();
        string tag = GetDomainAlgebraTag<TDomain,TAlgebra>();
        
        typedef ApproximationSpace<TDomain> approximation_space_type;
        typedef GridFunction<TDomain, TAlgebra> function_type;
        
        
        //	ImmersedInterfaceDiffusion
        {
            
            typedef ImmersedInterfaceDiffusion<TDomain, TAlgebra> T;
            typedef IImmersedInterface<TDomain, TAlgebra> TBase;
            string name = string("ImmersedInterfaceDiffusion").append(suffix);
            reg.add_class_<T, TBase>(name, grp)
            .template add_constructor<void (*)(SmartPtr<IAssemble<TAlgebra> > ass,
                                               SmartPtr<ConvectionDiffusionPlugin::ConvectionDiffusionFV1_cutElem<TDomain> > spMaster,
                                               SmartPtr<DiffusionInterfaceProvider<TDomain::dim> > interfaceProvider,
                                               SmartPtr<CutElementHandler_TwoSided<TDomain::dim> > cutElementHandler)>("domain disc, global handler")
            .add_method("init", &T::init)
            .add_method("set_source_data_lua", &T::set_source_data_lua)
            .add_method("set_jump_data_lua", &T::set_jump_data_lua)
            .add_method("set_diffusion_data_lua", &T::set_diffusion_data_lua)
            .add_method("set_jump_grad_data_lua", &T::set_jump_grad_data_lua)
            .add_method("get_L2Error", &T::get_L2Error)
            .add_method("get_numDoFs", &T::get_numDoFs)
            .add_method("set_Nitsche", &T::set_Nitsche)
            .add_method("set_print_cutElemData", &T::set_print_cutElemData)
            .add_method("get_numCutElements", &T::get_numCutElements)
            .add_method("adjust_for_error", &T::adjust_for_error)
            .add_method("initialize_threshold", &T::initialize_threshold)
            .add_method("set_threshold", &T::set_threshold)
            .add_method("set_analytic_solution", &T::set_analytic_solution)
            .set_construct_as_smart_pointer(true);
            reg.add_class_to_group(name, "ImmersedInterfaceDiffusion", tag);
        }
    }
    
        
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();
    
	//	singular sources and sinks
	{
		typedef CDSingularSourcesAndSinks<dim> T;
		typedef typename T::point_sss_type TPointSSS;
		typedef typename T::line_sss_type TLineSSS;

		string point_name = string("CDPointSourcesSink").append(dimSuffix);
		reg.add_class_<TPointSSS>(point_name, grp)
			.template add_constructor<void (*) (const std::vector<number>&)> ()
			.add_method ("set", static_cast<void (TPointSSS::*) (number)> (&TPointSSS::set))
			.add_method ("set", static_cast<void (TPointSSS::*) (typename TPointSSS::user_data_type)> (&TPointSSS::set))
#ifdef UG_FOR_LUA
			.add_method ("set", static_cast<void (TPointSSS::*) (LuaFunctionHandle)> (&TPointSSS::set))
#endif
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(point_name, "CDPointSourcesSink", dimTag);

		string line_name = string("CDLineSourcesSink").append(dimSuffix);
		reg.add_class_<TLineSSS>(line_name, grp)
			.template add_constructor<void (*) (const std::vector<number>&, const std::vector<number>&)> ()
			.add_method ("set", static_cast<void (TLineSSS::*) (number)> (&TLineSSS::set))
			.add_method ("set", static_cast<void (TLineSSS::*) (LuaFunctionHandle)> (&TLineSSS::set))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(line_name, "CDLineSourcesSink", dimTag);

		string name = string("CDSingularSourcesAndSinks").append(dimSuffix);
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method ("add_point", static_cast<void (T::*) (SmartPtr<TPointSSS>)> (&T::add_point))
			.add_method ("add_line", static_cast<void (T::*) (SmartPtr<TLineSSS>)> (&T::add_line))
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "CDSingularSourcesAndSinks", dimTag);
	}
}

}; // end Functionality

/**
 * Class exporting the functionality of the plugin restricted to 2 and 3 spatial
 * dimensions. All functionality that is to be used in scripts or visualization
 * only in 2d and 3d must be registered here.
 */
struct Functionality2d3d
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

//	Convection Diffusion FV1 for the low-dimensional fractures
	{
		typedef ConvectionDiffusionFractFV1<TDomain> T;
		typedef ConvectionDiffusionBase<TDomain> TBase;
		string name = string("ConvectionDiffusionFractFV1").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_fract_manager", static_cast<void (T::*)(SmartPtr<DegeneratedLayerManager<dim> >)>(&T::set_fract_manager), "Sets the fracture manager", "Deg. fracture manager")
			.add_method("set_upwind", &T::set_upwind)
			.add_method("set_aperture", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_aperture), "", "Fract. aperture")
			.add_method("set_aperture", static_cast<void (T::*)(number)>(&T::set_aperture), "", "Fract. aperture")
#ifdef UG_FOR_LUA
			.add_method("set_aperture", static_cast<void (T::*)(const char*)>(&T::set_aperture), "", "Fract. aperture")
			.add_method("set_aperture", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_aperture), "", "Fract. aperture")
#endif
			.add_method("set_ortho_velocity", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_ortho_velocity), "", "Orthogonal Velocity Field")
			.add_method("set_ortho_velocity", static_cast<void (T::*)(number)>(&T::set_ortho_velocity), "", "Orthogonal Velocity Field")
#ifdef UG_FOR_LUA
			.add_method("set_ortho_velocity", static_cast<void (T::*)(const char*)>(&T::set_ortho_velocity), "", "Orthogonal Velocity Field")
			.add_method("set_ortho_velocity", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_ortho_velocity), "", "Orthogonal Velocity Field")
#endif
			.add_method("set_ortho_diffusion", static_cast<void (T::*)(SmartPtr<CplUserData<number, dim> >)>(&T::set_ortho_diffusion), "", "Orthogonal Diffusion")
			.add_method("set_ortho_diffusion", static_cast<void (T::*)(number)>(&T::set_ortho_diffusion), "", "Orthogonal Diffusion")
#ifdef UG_FOR_LUA
			.add_method("set_ortho_diffusion", static_cast<void (T::*)(const char*)>(&T::set_ortho_diffusion), "", "Orthogonal Diffusion")
			.add_method("set_ortho_diffusion", static_cast<void (T::*)(LuaFunctionHandle)>(&T::set_ortho_diffusion), "", "Orthogonal Diffusion")
#endif
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusionFractFV1", tag);
	}
    
    // CutElementHandler_TwoSided
    {
        typedef CutElementHandler_TwoSided<dim> T;
        string name = string("CutElementHandler_TwoSided").append(suffix);
        reg.add_class_<T>(name, grp)
        .template add_constructor<void (*)(SmartPtr<MultiGrid> mg, const char*, SmartPtr<DiffusionInterfaceProvider<dim> >)>("multigrid, fct names")
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "CutElementHandler_TwoSided", tag);
    }
    
    // DiffusionInterfaceProvider
    {
        typedef DiffusionInterfaceProvider<dim> T;
        string name = string("DiffusionInterfaceProvider").append(suffix);
        reg.add_class_<T>(name, grp)
        .template add_constructor<void (*)( )>("")
        .add_method("print", &T::print)
        .add_method("add", &T::add)
        .set_construct_as_smart_pointer(true);
        reg.add_class_to_group(name, "DiffusionInterfaceProvider", tag);
    }
    
    
}

}; // end Functionality2d3d

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
	typedef ConvectionDiffusionPlugin::Functionality2d3d Functionality2d3d;

	try{
		RegisterDimensionDependent<Functionality>(*reg,grp);
		RegisterDomainDependent<Functionality>(*reg,grp);
        RegisterDomainAlgebraDependent<Functionality>(*reg,grp);
		RegisterDomain2d3dDependent<Functionality2d3d>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
