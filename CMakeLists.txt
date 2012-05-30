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
set(SOURCES	convection_diffusion_plugin.cpp
            convection_diffusion.cpp)



if(STATIC)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/core/convection_diffusion ${SOURCES})
else(STATIC)
	add_library(ConvectionDiffusion SHARED ${SOURCES})
	target_link_libraries (ConvectionDiffusion ug4)
endif(STATIC)
