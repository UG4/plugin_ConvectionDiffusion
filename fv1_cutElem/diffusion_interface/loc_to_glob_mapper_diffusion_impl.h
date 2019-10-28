/*
 * diffusion_interface_handler_local_impl.h
 *
 *  Created on: 15.01.2015
 *      Author: suze
 */

#ifndef DIFFUSION_INTERFACE_MAPPER_IMPL_H_
#define DIFFUSION_INTERFACE_MAPPER_IMPL_H_


namespace ug{
namespace ConvectionDiffusionPlugin{


template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
set_identity_mat_constraint(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{

	const LocalIndices& rowInd = lmat.get_row_indices();

	for (size_t fct1 = 0; fct1 < lmat.num_all_row_fct(); ++fct1)
		for (size_t dof1 = 0; dof1 < lmat.num_all_row_dof(fct1); ++dof1)
		{
 			const size_t rowIndex = rowInd.index(fct1, dof1);
			const size_t rowComp = rowInd.comp(fct1, dof1);

			if ( fabs(lmat.value(fct1, dof1, fct1, dof1)) < 0.000001
              && BlockRef(mat(rowIndex, rowIndex), rowComp, rowComp) == 0.0 )
				BlockRef(mat(rowIndex, rowIndex), rowComp, rowComp) = 1.0;
		}

}


template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
adjust_mat_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
	DoFIndex index;

	for (size_t i = 0; i < m_numGridNodes; ++i)
	{
      // set all entries, which were not assembled by the ElemDisc to the identity:
		if ( BlockRef(mat(i, i), 0, 0) == 0.0 )
			BlockRef(mat(i, i), 0, 0) = 1.0;
	}
}


template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
modify_GlobalSol(vector_type& vecMod, const vector_type& vec, ConstSmartPtr<DoFDistribution> dd)
{
	size_t numDoFs = vecMod.size();
	const size_t numNewDoFs = numDoFs - m_numGridNodes;

	UG_LOG("---------------modify_GlobalSol--------------\n");
	UG_LOG(" vecMod.size(): " << numDoFs << "\n");
 	UG_LOG(" m_numGridNodes: " << m_numGridNodes << "\n");
	UG_LOG(" computed numNewDoFs: " << numNewDoFs << "\n");
	UG_LOG(" m_numNewDoFs: " << m_numNewDoFs << "\n");

	DoFIndex index;
    std::vector<double > verticesValue;
    verticesValue.clear();

 // loop all interface nodes, i.e new DoFs, and store the values in
 // member 'm_verticesValue' of class 'InterfaceHandlerLocalDiffusion':
	for (size_t i = 0; i < numNewDoFs; ++i)
	{
    // get GLOBAL index of the interface nodes
        index = DoFIndex(m_numGridNodes + i,0);
    // write value to local data
		verticesValue.push_back(DoFRef(vec, index));
	}

// write solution values to member 'm_verticesValue' of class 'InterfaceHandlerLocalDiffusion':
 	set_interface_values(verticesValue);

}

template<typename TDomain, typename TAlgebra>
const size_t DiffusionInterfaceMapper<TDomain, TAlgebra>::
get_resize_measure(const size_t numDoFs)
{
    size_t numAllDoFs = m_numGridNodes + m_numNewDoFs;
    
  // if m_scaleDoFs = true: 2 new DoFs will be placed for each interface node
    if ( m_scaleDoFs )
        numAllDoFs = m_numGridNodes + 2*m_numNewDoFs;
    
    return numAllDoFs - numDoFs;
    
}
    
template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
// reset print modus for cut element data to 'false'
//  => only printed during assemling of defect, i.e. ONE loop over all element
    m_spInterfaceHandlerLocal->set_print_cutElemData(false);
    
    size_t numAllDoFs = m_numGridNodes + m_numNewDoFs;

// resize global matrix dynamically:
    const size_t diffDoFs = get_resize_measure(mat.num_rows());
    
// resize global defect only ONCE:
	if ( diffDoFs > 0 )
	{
		mat.resize_and_keep_values(numAllDoFs, numAllDoFs);
	}
	else if (  diffDoFs == 0  )
	{  // no resizing
    }
	else if ( diffDoFs < 0 )
	{	UG_THROW("error in add_local_mat_to_global: diffDofs < 0\n"); }


	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

	switch(modus)
	{
		case OUTSIDE_DOM: // call usual local-to-global-mapping
			AddLocalMatrixToGlobal(mat, lmat);
			break;
		case INSIDE_DOM: // call usual local-to-global-mapping
			AddLocalMatrixToGlobal(mat, lmat);
			break;
		case CUT_BY_INTERFACE: // call adapted local-to-global-mapping
			if ( m_spInterfaceHandlerLocal->get_Nitsche() )
				add_local_mat_to_global_Nitsche(mat, lmat, dd);
			else
				add_local_mat_to_global_interface(mat, lmat, dd);
			break;
 		default:
			throw(UGError("Error in IInterfaceMapper::add_local_mat_to_global()!"));

	}

}


template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
add_local_vec_to_global(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
{
    size_t numAllDoFs = m_numGridNodes + m_numNewDoFs;

// resize global vector dynamically:
    const size_t diffDoFs = get_resize_measure(vec.size());
   
// resize global defect ONCE:
	if ( diffDoFs > 0 )
	{
		vec.resize(numAllDoFs, true);
        vec.set(0.0);
	}
	else if (  diffDoFs == 0  )
	{  // no resizing
    }
	else if ( diffDoFs < 0 )
	{  UG_THROW("error in add_local_vec_to_global: diffDofs < 0\n"); }


	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

	switch(modus)
	{
		case OUTSIDE_DOM: // call usual local-to-global-mapping
			AddLocalVector(vec, lvec);
			break;
		case INSIDE_DOM: // call usual local-to-global-mapping
			AddLocalVector(vec, lvec);
			break;
		case CUT_BY_INTERFACE: // call adapted local-to-global-mapping
			if ( m_spInterfaceHandlerLocal->get_Nitsche() )
				add_local_vec_to_global_Nitsche(vec, lvec, dd);
			else
				add_local_vec_to_global_interface(vec, lvec, dd);
			break;
 		default:
			throw(UGError("Error in IInterfaceMapper::add_local_vec_to_global()!"));

	}
}

template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
add_local_mat_to_global_Nitsche(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	DoFIndex indexRow, indexCol;

///////////////////////////////////////////////////////////////
/// FIRST: add loc to glob for locJ_tri:
///////////////////////////////////////////////////////////////

    LocalMatrix locJ_tri = get_local_jacobian(ROID_TRIANGLE);
 	const bool shift_global_index_tri = m_spInterfaceHandlerLocal->get_index_shift_tri();

	if ( lmat.num_all_row_dof(0) != locJ_tri.num_all_row_dof(0) )
		UG_THROW("in 'add_local_mat_to_global_Nitsche': non-consistent sizees!\n");

 	size_t numAllDoFs = m_numGridNodes;

 	if ( shift_global_index_tri )
 		numAllDoFs = m_numGridNodes + m_numNewDoFs;

	for (size_t dof1 = 0; dof1 < locJ_tri.num_all_row_dof(0); ++dof1)
	{
		size_t global_index1 = rowInd.index(0, dof1);

	// if global_index1 is index of m_vertex: compute global index:
		if ( m_spInterfaceHandlerLocal->lies_onInterface(dof1) )
			global_index1 = numAllDoFs + m_spInterfaceHandlerLocal->get_index_for_global_index_Nitsche(global_index1);

		indexRow = DoFIndex(global_index1, rowInd.comp(0, dof1));

		for (size_t dof2 = 0; dof2 < locJ_tri.num_all_col_dof(0); ++dof2)
		{
			size_t global_index2 = rowInd.index(0, dof2);

		// if global_index2 is index of m_vertex: compute global index:
			if ( m_spInterfaceHandlerLocal->lies_onInterface(dof2) )
				global_index2 = numAllDoFs + m_spInterfaceHandlerLocal->get_index_for_global_index_Nitsche(global_index2);

			indexCol = DoFIndex(global_index2, colInd.comp(0, dof2));

		// finally add loc to glob:
			DoFRef(mat, indexRow, indexCol) += locJ_tri.value(0, dof1, 0, dof2);

		}
	}

///////////////////////////////////////////////////////////////
/// SECOND: add loc to glob for locJ_tri:
///    -> same as FIRST, but:
///			(1) lmat instead of locJ_tri
///		    (2) if ( !m_spInterfaceHandlerLocal->lies_onInterface(dof) )
///				instead of
///				if ( m_spInterfaceHandlerLocal->lies_onInterface(dof) )
///////////////////////////////////////////////////////////////

	for (size_t dof1 = 0; dof1 < lmat.num_all_row_dof(0); ++dof1)
	{
		size_t global_index1 = rowInd.index(0, dof1);

	// if global_index1 is index of m_vertex: compute global index:
		if ( !m_spInterfaceHandlerLocal->lies_onInterface(dof1) )
			global_index1 = numAllDoFs + m_spInterfaceHandlerLocal->get_index_for_global_index_Nitsche(global_index1);

		indexRow = DoFIndex(global_index1, rowInd.comp(0, dof1));

		for (size_t dof2 = 0; dof2 < lmat.num_all_col_dof(0); ++dof2)
		{
			size_t global_index2 = rowInd.index(0, dof2);

		// if global_index2 is index of m_vertex: compute global index:
			if ( !m_spInterfaceHandlerLocal->lies_onInterface(dof2) )
				global_index2 = numAllDoFs + m_spInterfaceHandlerLocal->get_index_for_global_index_Nitsche(global_index2);

			indexCol = DoFIndex(global_index2, colInd.comp(0, dof2));

		// finally add loc to glob:
			DoFRef(mat, indexRow, indexCol) += lmat.value(0, dof1, 0, dof2);

		}
	}


}

//get_real_index(dof):
// if dof NOT on interface: returns orig corner index
// if dof ON interface: returns index of m_vertex-array

template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
add_local_mat_to_global_interface(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	DoFIndex indexRow, indexCol;

/////////////////////////////////////////////////////////////////////////////////
/// FIRST: add the boundary contribution stored in 'locJ_tri' to global matrix:
/////////////////////////////////////////////////////////////////////////////////
    
    LocalMatrix locJ_tri = get_local_jacobian(ROID_TRIANGLE);
 	const bool shift_global_index_tri = m_spInterfaceHandlerLocal->get_index_shift_tri();

 	size_t numAllDoFs = m_numGridNodes;

 	if ( shift_global_index_tri )
 		numAllDoFs = m_numGridNodes + m_numNewDoFs;

// loop all entries for mapping
	for (size_t dof1 = 0; dof1 < locJ_tri.num_all_row_dof(0); ++dof1)
	{
    // get the real index: this data was stored in 'InterfaceHandlerDiffusion::m_vRealCornerID_'
    // case1: node lies on interface => real_dof = entry within map 'InterfaceHandlerDiffusion::m_MapNearVertices'
    // case2: node lies on an original mesh node: => real_dof = usual, local index of vertex
		const size_t dof1_real = m_spInterfaceHandlerLocal->real_index_tri(dof1);

	// case1: dof1_real is index of m_vertex (i.e. interface-DoF): compute global index:
		if ( m_spInterfaceHandlerLocal->lies_onInterface_tri(dof1) )
			indexRow = DoFIndex(numAllDoFs + dof1_real, 0);
		else // case2: use dof1_real as local index
			indexRow = DoFIndex(rowInd.index(0, dof1_real), rowInd.comp(0, dof1_real));

		for (size_t dof2 = 0; dof2 < locJ_tri.num_all_col_dof(0); ++dof2)
		{
			const size_t dof2_real = m_spInterfaceHandlerLocal->real_index_tri(dof2);

		// if dof1_real is index of m_vertex (i.e. interface-DoF): compute global index:
			if ( m_spInterfaceHandlerLocal->lies_onInterface_tri(dof2) )
				indexCol = DoFIndex(numAllDoFs + dof2_real, 0);
			else
				indexCol = DoFIndex(colInd.index(0, dof2_real), colInd.comp(0, dof2_real));

		// finally add loc to glob:
			DoFRef(mat, indexRow, indexCol) += locJ_tri.value(0, dof1, 0, dof2);
		}
	}

////////////////////////////////////////////////////////////////////////////////////
/// SECOND: add the boundary contribution stored in 'locJ_quad' to global matrix:
////////////////////////////////////////////////////////////////////////////////////
    
    LocalMatrix locJ_quad = get_local_jacobian(ROID_QUADRILATERAL);
  	const bool shift_global_index_quad = m_spInterfaceHandlerLocal->get_index_shift_quad();

// reset numAllDoFs!
  	 numAllDoFs = m_numGridNodes;

 	if ( shift_global_index_quad )
 		numAllDoFs = m_numGridNodes + m_numNewDoFs;

// loop all entries for mapping
	for (size_t dof1 = 0; dof1 < locJ_quad.num_all_row_dof(0); ++dof1)
	{
		size_t dof1_real = m_spInterfaceHandlerLocal->real_index_quad(dof1);

    // if dof1_real is index of m_vertex: compute global index:
        if ( m_spInterfaceHandlerLocal->lies_onInterface_quad(dof1)
         && !m_spInterfaceHandlerLocal->check_vertex_modus(ON_INTERFACE, dof1, -1) )
 			indexRow = DoFIndex(numAllDoFs + dof1_real, 0);
		else
			indexRow = DoFIndex(rowInd.index(0, dof1_real), rowInd.comp(0, dof1_real));

		for (size_t dof2 = 0; dof2 < locJ_quad.num_all_col_dof(0); ++dof2)
		{
			size_t dof2_real = m_spInterfaceHandlerLocal->real_index_quad(dof2);

        // if dof1_real is index of m_vertex: compute global index:
            if ( m_spInterfaceHandlerLocal->lies_onInterface_quad(dof2)
             && !m_spInterfaceHandlerLocal->check_vertex_modus(ON_INTERFACE, dof2, -1) )
 				indexCol = DoFIndex(numAllDoFs + dof2_real, 0);
			else
				indexCol = DoFIndex(colInd.index(0, dof2_real), colInd.comp(0, dof2_real));

		// finally add loc to glob:
			DoFRef(mat, indexRow, indexCol) += locJ_quad.value(0, dof1, 0, dof2);

		}
	}

}


template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
add_local_vec_to_global_Nitsche(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
{
	const LocalIndices& ind = lvec.get_indices();
	DoFIndex index;

	///////////////////////////////////////////////////////////////
	/// FIRST: add loc to glob for 'locD_tri':
	///////////////////////////////////////////////////////////////
    LocalVector locD_tri = get_local_defect(ROID_TRIANGLE);
	const bool shift_global_index_tri = m_spInterfaceHandlerLocal->get_index_shift_tri();

	if ( lvec.num_all_dof(0) != locD_tri.num_all_dof(0) )
		UG_THROW("in 'add_local_vec_to_global_Nitsche': non-consistent sizees!\n");

 	size_t numAllDoFs = m_numGridNodes;

 	if ( shift_global_index_tri )
 	{		numAllDoFs = m_numGridNodes + m_numNewDoFs;
 	}

	for (size_t dof = 0; dof < locD_tri.num_all_dof(0); ++dof)
	{
		size_t global_index = ind.index(0, dof);

	// if dof_real is index of m_vertex: compute global index:
		if ( m_spInterfaceHandlerLocal->lies_onInterface(dof) )
			global_index = numAllDoFs + m_spInterfaceHandlerLocal->get_index_for_global_index_Nitsche(global_index);

		index = DoFIndex(global_index, ind.comp(0, dof));

	// finally add loc to glob:
		DoFRef(vec, index) += locD_tri.value(0, dof);

	}


	///////////////////////////////////////////////////////////////
	/// SECOND: add loc to glob for 'lvec':
	///    -> same as FIRST, but:
	///			(1) lvec instead of locD_tri
	///		    (2) if ( !m_spInterfaceHandlerLocal->lies_onInterface(dof) )
	///				instead of
	///				if ( m_spInterfaceHandlerLocal->lies_onInterface(dof) )
	///////////////////////////////////////////////////////////////

	for (size_t dof = 0; dof < lvec.num_all_dof(0); ++dof)
	{
		size_t global_index = ind.index(0, dof);

	// if dof_real is index of m_vertex: compute global index:
		if ( !m_spInterfaceHandlerLocal->lies_onInterface(dof) )
			global_index = numAllDoFs + m_spInterfaceHandlerLocal->get_index_for_global_index_Nitsche(global_index);

		index = DoFIndex(global_index, ind.comp(0, dof));

	// finally add loc to glob:
		DoFRef(vec, index) += lvec.value(0, dof);

	}

}

template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
add_local_vec_to_global_interface(vector_type& vec, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd)
{
	DoFIndex index_print;

	const LocalIndices& ind = lvec.get_indices();
	DoFIndex index;

///////////////////////////////////////////////////////////////////////////////////
/// FIRST: add the boundary contribution stored in 'locD_tri' to global defect:
///////////////////////////////////////////////////////////////////////////////////

    LocalVector locD_tri = get_local_defect(ROID_TRIANGLE);
	const bool shift_global_index_tri = m_spInterfaceHandlerLocal->get_index_shift_tri();

 	size_t numAllDoFs = m_numGridNodes;

 	if ( shift_global_index_tri )
        numAllDoFs = m_numGridNodes + m_numNewDoFs;

//	UG_LOG("in DiffusionInterfaceMapper::add_loc_vec(): locD_tri = " << locD_tri << "\n");

	for (size_t dof = 0; dof < locD_tri.num_all_dof(0); ++dof)
	{
		const size_t dof_real = m_spInterfaceHandlerLocal->real_index_tri(dof);

	// if dof_real is index of m_vertex: compute global index:
        if ( m_spInterfaceHandlerLocal->lies_onInterface_tri(dof) && !m_spInterfaceHandlerLocal->check_vertex_modus(ON_INTERFACE, dof, 1) )
            index = DoFIndex(numAllDoFs + dof_real,0);
		else
			index = DoFIndex(ind.index(0, dof_real), ind.comp(0, dof_real));

	// finally add loc to glob:
		DoFRef(vec, index) += locD_tri.value(0, dof);

	}


/////////////////////////////////////////////////////////////////////////////////
/// SECOND: add the boundary contribution stored in 'locD_quad' to global defect:
/////////////////////////////////////////////////////////////////////////////////

    LocalVector locD_quad = get_local_defect(ROID_QUADRILATERAL);
	const bool shift_global_index_quad = m_spInterfaceHandlerLocal->get_index_shift_quad();

// reset numAllDoFs!
    numAllDoFs = m_numGridNodes;
    
 	if ( shift_global_index_quad )
 		numAllDoFs = m_numGridNodes + m_numNewDoFs;

//	UG_LOG("in DiffusionInterfaceMapper::add_loc_vec(): locD_quad = " << locD_quad << "\n");

	for (size_t dof = 0; dof < locD_quad.num_all_dof(0); ++dof)
	{
		size_t dof_real = m_spInterfaceHandlerLocal->real_index_quad(dof);

	// if dof_real is index of m_vertex: compute global index:
        if ( m_spInterfaceHandlerLocal->lies_onInterface_quad(dof)
         && !m_spInterfaceHandlerLocal->check_vertex_modus(ON_INTERFACE, dof, 1) )
            index = DoFIndex(numAllDoFs + dof_real,0);
  		else
			index = DoFIndex(ind.index(0, dof_real), ind.comp(0, dof_real));

	// finally add loc to glob:
		DoFRef(vec, index) += locD_quad.value(0, dof);
	}

    
}

} // namespace ConvectionDiffusionPlugin
} // end ug namespace

#endif /* DIFFUSION_INTERFACE_MAPPER_IMPL_H_ */
