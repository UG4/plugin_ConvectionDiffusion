# Copyright (c) 2012-2013:  G-CSC, Goethe University Frankfurt
# Author: Andreas Vogel
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

################################################################################
# ConvectionDiffusion
#
# The ConvectionDiffusion plugin
################################################################################

# Name of your plugin and sources.
set(pluginName	ConvectionDiffusion)
set(SOURCES		convection_diffusion_base.cpp
				fv1/convection_diffusion_fv1.cpp
				fe/convection_diffusion_fe.cpp
				fe/convection_diffusion_stab_fe.cpp
				fvcr/convection_diffusion_fvcr.cpp
				fv/convection_diffusion_fv.cpp
				fractfv1/convection_diffusion_fractfv1.cpp
				convection_diffusion_plugin.cpp)


################################################################################
# The code below doesn't have to be changed (usually)
################################################################################
cmake_minimum_required(VERSION 2.8.12...3.27.1)

project(UG_PLUGIN_${pluginName})

################################################################################
# Classic binding.
################################################################################
if(NOT USE_PYBIND11)
if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else(buildEmbeddedPlugins)

	# Create a SHARED library from the sources and link it against ug4.
	add_library(${pluginName} SHARED ${SOURCES})
	target_link_libraries (${pluginName} ug4)
	
	# CPack
	set_target_properties(${pluginName} PROPERTIES INSTALL_RPATH "$ORIGIN/../../lib/")
	install(TARGETS ${pluginName} LIBRARY DESTINATION bin/plugins COMPONENT plugins)
	
endif(buildEmbeddedPlugins)
endif(NOT USE_PYBIND11)


################################################################################
# Python binding (plugin is static, python interface dynamic).
################################################################################
if(USE_PYBIND11)

	# Create a STATIC library from the sources.
	add_library(${pluginName} STATIC ${SOURCES})

	SET(myPluginSources convection_diffusion_pybind.cpp)
	SET(myLibraries ${pluginName} ${targetLibraryName}) # targetLibraryName is either ug4 or ug4_s
	
	# First argument must match module name in PYBIND11_MODULE call
	ug4pybind_add_module(pyconvectiondiffusion ${myPluginSources} ${myLibraries})
	target_link_libraries (pyconvectiondiffusion PRIVATE ${pluginName}) # TODO: should be part of ug4pybind_add_module
	
	# CPack
	set_target_properties(pyconvectiondiffusion PROPERTIES INSTALL_RPATH "$ORIGIN/..:$ORIGIN/../../../lib")
	install(TARGETS pyconvectiondiffusion LIBRARY DESTINATION python COMPONENT pymodules)
endif(USE_PYBIND11)