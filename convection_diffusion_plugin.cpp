/**
 * author: andreasvogel
 *
 * File for registration of ConvectionDiffusion routines.
 *
 */

#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"
#include "convection_diffusion.h"

using namespace std;
using namespace ug::bridge;

namespace ug{
namespace ConvectionDiffusionPlugin{

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

//	Convection Diffusion
	{
		typedef ConvectionDiffusion<TDomain> T;
		typedef IDomainElemDisc<TDomain> TBase;
		string name = string("ConvectionDiffusion").append(suffix);
		reg.add_class_<T, TBase >(name, grp)
			.template add_constructor<void (*)(const char*,const char*)>("Function(s)#Subset(s)")
			.add_method("set_disc_scheme", &T::set_disc_scheme, "", "Disc Scheme|selection|value=[\"fe\",\"fv\",\"fv1\",\"fvcr\"]")
			.add_method("set_quad_order", &T::set_quad_order)
			.add_method("set_quad_order_scvf", &T::set_quad_order_scvf)
			.add_method("set_quad_order_scv", &T::set_quad_order_scv)

			.add_method("set_diffusion", static_cast<void (T::*)(SmartPtr<UserData<MathMatrix<dim, dim>, dim> >)>(&T::set_diffusion), "", "Diffusion")
			.add_method("set_diffusion", static_cast<void (T::*)(number)>(&T::set_diffusion), "", "Diagonal Diffusion")
#ifdef UG_FOR_LUA
			.add_method("set_diffusion", static_cast<void (T::*)(const char*)>(&T::set_diffusion), "", "Diffusion")
#endif

			.add_method("set_velocity", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(const std::vector<number>&)>(&T::set_velocity), "", "Velocity Field")
			.add_method("set_velocity", static_cast<void (T::*)(number)>(&T::set_velocity), "", "Vel_x")
			.add_method("set_velocity", static_cast<void (T::*)(number,number)>(&T::set_velocity), "", "Vel_x, Vel_y")
			.add_method("set_velocity", static_cast<void (T::*)(number,number,number)>(&T::set_velocity), "", "Vel_x, Vel_y, Vel_z")
#ifdef UG_FOR_LUA
			.add_method("set_velocity", static_cast<void (T::*)(const char*)>(&T::set_velocity), "", "Velocity Field")
#endif

			.add_method("set_reaction_rate", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_reaction_rate), "", "Reaction Rate")
			.add_method("set_reaction_rate", static_cast<void (T::*)(number)>(&T::set_reaction_rate), "", "Reaction Rate")
#ifdef UG_FOR_LUA
			.add_method("set_reaction_rate", static_cast<void (T::*)(const char*)>(&T::set_reaction_rate), "", "Reaction Rate")
#endif

			.add_method("set_reaction", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_reaction), "", "Reaction")
			.add_method("set_reaction", static_cast<void (T::*)(number)>(&T::set_reaction), "", "Reaction")
#ifdef UG_FOR_LUA
			.add_method("set_reaction", static_cast<void (T::*)(const char*)>(&T::set_reaction), "", "Reaction")
#endif

			.add_method("set_source", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_source), "", "Source")
			.add_method("set_source", static_cast<void (T::*)(number)>(&T::set_source), "", "Source")
#ifdef UG_FOR_LUA
			.add_method("set_source", static_cast<void (T::*)(const char*)>(&T::set_source), "", "Source")
#endif

			.add_method("set_vector_source", static_cast<void (T::*)(SmartPtr<UserData<MathVector<dim>, dim> >)>(&T::set_vector_source), "", "Vector Source")
			.add_method("set_vector_source", static_cast<void (T::*)(number)>(&T::set_vector_source), "", "vectorSource_x")
			.add_method("set_vector_source", static_cast<void (T::*)(number,number)>(&T::set_vector_source), "", "vectorSource_x, vectorSource_y")
			.add_method("set_vector_source", static_cast<void (T::*)(number,number,number)>(&T::set_vector_source), "", "vectorSource_x, vectorSource_y, vectorSource_z")
#ifdef UG_FOR_LUA
			.add_method("set_vector_source", static_cast<void (T::*)(const char*)>(&T::set_vector_source), "", "Vector Source")
#endif

			.add_method("set_mass_scale", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_mass_scale), "", "Mass Scale")
			.add_method("set_mass_scale", static_cast<void (T::*)(number)>(&T::set_mass_scale), "", "Mass Scale")
#ifdef UG_FOR_LUA
			.add_method("set_mass_scale", static_cast<void (T::*)(const char*)>(&T::set_mass_scale), "", "Mass Scale")
#endif

			.add_method("set_mass", static_cast<void (T::*)(SmartPtr<UserData<number, dim> >)>(&T::set_mass), "", "Mass")
			.add_method("set_mass", static_cast<void (T::*)(number)>(&T::set_mass), "", "Mass")
#ifdef UG_FOR_LUA
			.add_method("set_mass", static_cast<void (T::*)(const char*)>(&T::set_mass), "", "Mass")
#endif

			.add_method("set_upwind", &T::set_upwind)
			.add_method("value", &T::value)
			.add_method("gradient", &T::gradient)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConvectionDiffusion", tag);
	}

}

}; // end Functionality
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
		RegisterDomainDependent<Functionality>(*reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// namespace ug
