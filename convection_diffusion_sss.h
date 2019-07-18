/*
 * Copyright (c) 2019:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko / Michael Lampe
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

/*
 * Singular (point and line) sources and sinks in ConvectionDiffusion.
 */
#ifndef __H__UG__PLUGINS__CD__SINGULAR_SOURCES_AND_SINKS__
#define __H__UG__PLUGINS__CD__SINGULAR_SOURCES_AND_SINKS__

#include <vector>

// ug4 headers
#include "lib_disc/spatial_disc/disc_util/fv1_sss.h"
#include "lib_disc/spatial_disc/user_data/user_data.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug {
namespace ConvectionDiffusionPlugin {

/// class for data for all the CD plugin sources and sinks
template <int dim>
class cd_sss_data
{
	/** the data for the source/sink:
	 * [0]: total contaminant flux through the point
	 */
	MathVector<1> m_values;
	
	SmartPtr<UserData<MathVector<1>, dim> > m_spData; ///< an alternative method to specify the data
	
public:
	
///	class construction (there must exist a 'dummy' constructor!)
	cd_sss_data () {m_values [0] = 0;}
	
///	returns the flux
	number flux () {return m_values [0];}
	
///	computes the data from the user data object
	void compute
	(
		const MathVector<dim>& x, ///< point where to evaluate
		number time, ///< time argument for the evaluation
		int si ///< subset where to evaluate
	)
	{
		if (m_spData.valid ())
			(* m_spData) (m_values, x, time, si);
	}
	
///	sets the data
	void set (number flux)
	{
		m_values[0] = flux;
		m_spData = SPNULL;
	}
	
///	sets the data by an object
	void set (SmartPtr<UserData<MathVector<1>, dim> > spData)
	{
		m_spData = spData;
	}
	
///	set as a LUA function
	void set (LuaFunctionHandle func)
	{
		m_spData = make_sp (new LuaUserData<MathVector<1>, dim> (func));
	}
};

/** Class for markers of the point sources and sinks
 *
 * Note that there the point sinks are only used for full-dimensional subdomains.
 */
class point_sss_marker
{
	GridObject * m_elem; ///< grid element for the source/sink (not to take it into account twice)
	size_t m_co; ///< corner of the element (not to take it into account twice inside of the element)
	
public:

///	class constructor
	point_sss_marker () : m_elem (NULL), m_co (0) {};
	
///	resets the mark
	void init () {m_elem = NULL; m_co = 0;}
	
///	check and set the element mark
	bool marked_for (GridObject * elem, size_t co)
	{
		if (m_elem == NULL)
		{
			m_elem = elem; m_co = co;
			return true;
		}
		return m_elem == elem && m_co == co;
	}
};

/** Class for markers of the line sources and sinks
 *
 * Note that for fractures, line sources/sinks are point sources/sinks.
 * For full-dimensional subdomains, there is up to now no special markers.
 */
class line_sss_marker
{
// All the members are used in the fractures only!

///	a special structure to identify the element and its corner in a fracture
	struct t_fract_elem
	{
		IVertexGroup * fract_face;
		size_t fract_co;
		
		t_fract_elem (IVertexGroup * face, size_t co) : fract_face (face), fract_co (co) {};
	};
	
///	array keeping the elements from different(!) fractures
	std::vector<t_fract_elem> m_intersections;
	
public:

///	class constructor
	line_sss_marker () {};
	
///	reset the mark
	void init () {m_intersections.clear ();}
	
///	check and set the element mark (use it only for fractures!)
	bool marked_for (IVertexGroup * elem, size_t co)
	{
	//	is the fracture already processed?
		for (size_t i = 0; i < m_intersections.size (); i++)
		{
		//	check if we have already registered this mark
			t_fract_elem& intersection = m_intersections[i];
			if (intersection.fract_face == elem)
				return intersection.fract_co == co;
			
		//	check if a different corner is meant (i.e. this is the same fracture)
			Vertex * vrt = elem->vertex (co);
			for (size_t j = 0; j < intersection.fract_face->num_vertices (); j++)
				if (vrt == intersection.fract_face->vertex (j))
					return false; // in this fracture, we use a different corner
		}
	// no, register this fracture, too
		m_intersections.push_back (t_fract_elem (elem, co));
		return true;
	}
};

template <int dim> class cd_point_sss_data : public cd_sss_data<dim>, public point_sss_marker {};
template <int dim> class cd_line_sss_data : public cd_sss_data<dim>, public line_sss_marker {};
template <int dim>
class CDSingularSourcesAndSinks
	: public FVSingularSourcesAndSinks<dim, cd_point_sss_data<dim>, cd_line_sss_data<dim> >
{};

} // namespace ConvectionDiffusionPlugin
} // end namespace ug

#endif // __H__UG__PLUGINS__CD__SINGULAR_SOURCES_AND_SINKS__

/* End of File */
