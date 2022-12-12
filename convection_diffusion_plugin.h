#pragma once

#include <string>
#include "bridge/util.h"

extern "C" void InitUGPlugin_ConvectionDiffusion(ug::bridge::Registry* reg, std::string grp);

#ifdef UG_USE_PYBIND11
namespace ug {
namespace ConvectionDiffusionPlugin{

	void InitUGPlugin_ConvectionDiffusion(ug::pybind::Registry* reg, std::string grp);
}
}
#endif
