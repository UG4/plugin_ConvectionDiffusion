#include "convection_diffusion_plugin.h"

#ifdef UG_USE_PYBIND11
PYBIND11_MODULE(pyconvectiondiffusion, m)
{
	m.doc() = "Convection diffusion module";
	m.attr("__name__") = "ug4py.convection_diffusion";

	ug::pybind::Registry registry(m);
	std::string name("ConvDiff");

	ug::ConvectionDiffusionPlugin::InitUGPlugin_ConvectionDiffusion(&registry, name);
}
#endif
