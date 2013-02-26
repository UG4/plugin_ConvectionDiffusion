################################################################################
# ConvectionDiffusion
#
# The ConvectionDiffusion plugin
################################################################################

cmake_minimum_required(VERSION 2.6)

project(ConvectionDiffusion)

# include the definitions and dependencies for ug-plugins.
include("../../../ug_plugin_cmake_includes.txt")


# sources
set(SOURCES	convection_diffusion_base.cpp
            fv1/convection_diffusion_fv1.cpp
            fe/convection_diffusion_fe.cpp
            fvcr/convection_diffusion_fvcr.cpp
            fv/convection_diffusion_fv.cpp
            convection_diffusion_plugin.cpp)



if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/core/convection_diffusion ${SOURCES})
else(buildEmbeddedPlugins)
	add_library(ConvectionDiffusion SHARED ${SOURCES})
	target_link_libraries (ConvectionDiffusion ug4)
endif(buildEmbeddedPlugins)
