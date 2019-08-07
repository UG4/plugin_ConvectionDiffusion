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

			if ( fabs(lmat.value(fct1, dof1, fct1, dof1)) < 0.000001 && BlockRef(mat(rowIndex, rowIndex), rowComp, rowComp) == 0.0 )
				BlockRef(mat(rowIndex, rowIndex), rowComp, rowComp) = 1.0;
		}

}


template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
adjust_mat_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
	DoFIndex index;

	for (size_t i = 0; i < m_numDoFs; ++i)
	{
		if ( BlockRef(mat(i, i), 0, 0) == 0.0 )
			BlockRef(mat(i, i), 0, 0) = 1.0;
	}
}


template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
modify_GlobalSol(vector_type& vecMod, const vector_type& vec, ConstSmartPtr<DoFDistribution> dd)
{
	size_t numDoFs = vecMod.size();
	const size_t numNewDoFs = numDoFs - m_numDoFs;

	UG_LOG("---------------modify_GlobalSol--------------\n");
	UG_LOG(" vecMod.size(): " << numDoFs << "\n");
 	UG_LOG("m_numDoFs: " << m_numDoFs << "\n");
	UG_LOG(" computed numNewDoFs: " << numNewDoFs << "\n");
	UG_LOG(" m_numNewDoFs: " << m_numNewDoFs << "\n");

	DoFIndex index;
    std::vector<double > verticesValue;
    verticesValue.clear();

	for (size_t i = 0; i < numNewDoFs; ++i)
	{
		index = DoFIndex(m_numDoFs + i,0);
	// if vrt on interface == DoF: get value from vec
	// if vrt on interface != DoF: written via 'compute_values()'
		double value;
		if ( 1 )
		{	value = DoFRef(vec, index); }
		else
		{
			value = 0.0; //m_spInterfaceHandlerLocal->compute_value(m_spInterfaceHandlerLocal->m_verticesPos[i]);
		}

		verticesValue.push_back(value);
	}

	// call InterfaceHandlerLocal-method:
	// no not doing it! Done locally in add_def: locU_tri and locU_quad:
	write_solution(verticesValue);

}

template<typename TDomain, typename TAlgebra>
void DiffusionInterfaceMapper<TDomain, TAlgebra>::
add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat, ConstSmartPtr<DoFDistribution> dd)
{
	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	DoFIndex indexRow, indexCol;

	bool print = false;

// resize global matrix dynamically:
	const size_t numDoFs = mat.num_rows();

	if ( print ) UG_LOG("*** vorher: vec.size(): " << numDoFs << "\n");

	size_t numAllDoFs = m_numDoFs + m_numNewDoFs;
	if ( m_scaleDoFs )
		numAllDoFs = m_numDoFs + 2*m_numNewDoFs;

	const int diffDoFs = numAllDoFs - numDoFs;

// resize global defect ONCE:
	if ( diffDoFs > 0 )
	{
		mat.resize_and_keep_values(numAllDoFs, numAllDoFs);
		if ( print )
		{
			UG_LOG("*** m_numDoFs: " << m_numDoFs << "\n");
			UG_LOG("*** m_numNewDoFs: " << m_numNewDoFs << "\n");
			UG_LOG("*** numAllDoFs: " << numAllDoFs << "\n");
			UG_LOG("*** nachher: mat.num_rows(): " << mat.num_rows() << "\n");
		}
	}
	else if (  diffDoFs == 0  )
	{   if ( print ) UG_LOG("no resizing!\n");}
	else if ( diffDoFs < 0 )
	{
		if ( print ) UG_LOG("diffDoFs = " << diffDoFs << "\n");
		UG_THROW("error in add_local_mat_to_global: diffDofs < 0\n");
	}


	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

	switch(modus)
	{
		case OUTSIDE_DOM:
			AddLocalMatrixToGlobal(mat, lmat);
			break;
		case INSIDE_DOM:
			AddLocalMatrixToGlobal(mat, lmat);
			//set_identity_mat_constraint(mat, lmat, dd);
			break;
		case CUT_BY_INTERFACE:
            //AddLocalMatrixToGlobal(mat, lmat);
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
	bool print = false;

	const size_t numDoFs = vec.size();

	if ( print ) UG_LOG("*** vorher: vec.size(): " << numDoFs << "\n");

	size_t numAllDoFs = m_numDoFs + m_numNewDoFs;
	if ( m_scaleDoFs )
		numAllDoFs = m_numDoFs + 2*m_numNewDoFs;

	const int diffDoFs = numAllDoFs - numDoFs;

// resize global defect ONCE:
	if ( diffDoFs > 0 )
	{
		vec.resize(numAllDoFs, true);
        vec.set(0.0);
/*
        bool bJac = m_spInterfaceHandlerLocal->get_jac_bool();
        if ( !bJac )
        {
           // vec.set(0.0);
            m_spInterfaceHandlerLocal->set_jac_bool(true);
        }
*/
		if ( print )
		{
			UG_LOG("*** m_numDoFs: " << m_numDoFs << "\n");
			UG_LOG("*** m_numNewDoFs: " << m_numNewDoFs << "\n");
			UG_LOG("*** numAllDoFs: " << numAllDoFs << "\n");
			UG_LOG("*** nachher: vec.size(): " << vec.size() << "\n");
		}
	}
	else if (  diffDoFs == 0  )
	{   if ( print ) UG_LOG("no resizing!\n");}
	else if ( diffDoFs < 0 )
	{
		if ( print ) UG_LOG("diffDoFs = " << diffDoFs << "\n");
		UG_THROW("error in add_local_vec_to_global: diffDofs < 0\n");
	}


	ElementModus modus = m_spInterfaceHandlerLocal->elementModus();

	switch(modus)
	{
		case OUTSIDE_DOM:
			AddLocalVector(vec, lvec);
			break;
		case INSIDE_DOM:
			AddLocalVector(vec, lvec);
			break;
		case CUT_BY_INTERFACE:
 //           AddLocalVector(vec, lvec);
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

	UG_LOG("------------------------ START DiffusionInterfaceMapper - Nitsche ------------------------\n");

	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	DoFIndex indexRow, indexCol;

	///////////////////////////////////////////////////////////////
	/// FIRST: add loc to glob for locJ_tri:
	///////////////////////////////////////////////////////////////
	const LocalMatrix& locJ_tri = get_local_jacobian_tri();
 	const bool shift_global_index_tri = m_spInterfaceHandlerLocal->get_index_shift_tri();

	if ( lmat.num_all_row_dof(0) != locJ_tri.num_all_row_dof(0) )
		UG_THROW("in 'add_local_mat_to_global_Nitsche': non-consistent sizees!\n");

 	size_t numAllDoFs = m_numDoFs;

 	if ( shift_global_index_tri )
 		numAllDoFs = m_numDoFs + m_numNewDoFs;

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

	indexRow = DoFIndex(3, 0);
	indexCol = DoFIndex(0, 0);
	number value = DoFRef(mat, indexRow, indexCol);

	if ( DoFRef(mat, indexRow, indexCol) > 0 )
		UG_LOG("waiting...\n");
	UG_LOG("value: " << value << "\n");


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
	// UG_LOG("------------------------ START DiffusionInterfaceMapper ------------------------\n");

	const LocalIndices& rowInd = lmat.get_row_indices();
	const LocalIndices& colInd = lmat.get_col_indices();

	DoFIndex indexRow, indexCol;

	///////////////////////////////////////////////////////////////
	/// FIRST: add loc to glob for locJ_tri:
	///////////////////////////////////////////////////////////////
	const LocalMatrix& locJ_tri = get_local_jacobian_tri();
 	const bool shift_global_index_tri = m_spInterfaceHandlerLocal->get_index_shift_tri();

 	size_t numAllDoFs = m_numDoFs;

 	if ( shift_global_index_tri )
 		numAllDoFs = m_numDoFs + m_numNewDoFs;

//	UG_LOG("in DiffusionInterfaceMapper::add_loc_mat(): locJ_tri = " << locJ_tri << "\n");


	for (size_t dof1 = 0; dof1 < locJ_tri.num_all_row_dof(0); ++dof1)
	{
		const size_t dof1_real = m_spInterfaceHandlerLocal->real_index_tri(dof1);

	// if dof1_real is index of m_vertex: compute global index:
		if ( m_spInterfaceHandlerLocal->lies_onInterface_tri(dof1) )
			indexRow = DoFIndex(numAllDoFs + dof1_real, 0);
		else
			indexRow = DoFIndex(rowInd.index(0, dof1_real), rowInd.comp(0, dof1_real));

		for (size_t dof2 = 0; dof2 < locJ_tri.num_all_col_dof(0); ++dof2)
		{
			const size_t dof2_real = m_spInterfaceHandlerLocal->real_index_tri(dof2);

		// if dof1_real is index of m_vertex: compute global index:
			if ( m_spInterfaceHandlerLocal->lies_onInterface_tri(dof2) )
				indexCol = DoFIndex(numAllDoFs + dof2_real, 0);
			else
				indexCol = DoFIndex(colInd.index(0, dof2_real), colInd.comp(0, dof2_real));

		// finally add loc to glob:
			DoFRef(mat, indexRow, indexCol) += locJ_tri.value(0, dof1, 0, dof2);

			if ( indexRow[0] == 0 || indexRow[0] == 30 ){
				UG_LOG("------------------------ tri: += " << locJ_tri.value(0, dof1, 0, dof2) << "\n");
				UG_LOG("(indexRow, indexCol) = " << indexRow[0] << "," << indexCol[0] << "\n");
			}
		}
	}

	///////////////////////////////////////////////////////////////
	/// SECOND: add loc to glob for locJ_tri:
	///////////////////////////////////////////////////////////////
	const LocalMatrix& locJ_quad = get_local_jacobian_quad();
  	const bool shift_global_index_quad = m_spInterfaceHandlerLocal->get_index_shift_quad();

  	// reset numAllDoFs!
  	 numAllDoFs = m_numDoFs;

 	if ( shift_global_index_quad )
 		numAllDoFs = m_numDoFs + m_numNewDoFs;

//	UG_LOG("in DiffusionInterfaceMapper::add_loc_mat(): locJ_quad = " << locJ_quad << "\n");

	for (size_t dof1 = 0; dof1 < locJ_quad.num_all_row_dof(0); ++dof1)
	{
		size_t dof1_real = m_spInterfaceHandlerLocal->real_index_quad(dof1);

		// if dof1_real is index of m_vertex: compute global index:
        if ( m_spInterfaceHandlerLocal->lies_onInterface_quad(dof1) && !m_spInterfaceHandlerLocal->check_vertex_modus(ON_INTERFACE, dof1, -1) )
 			indexRow = DoFIndex(numAllDoFs + dof1_real, 0);
		else
			indexRow = DoFIndex(rowInd.index(0, dof1_real), rowInd.comp(0, dof1_real));

		for (size_t dof2 = 0; dof2 < locJ_quad.num_all_col_dof(0); ++dof2)
		{
			size_t dof2_real = m_spInterfaceHandlerLocal->real_index_quad(dof2);

			// if dof1_real is index of m_vertex: compute global index:
            if ( m_spInterfaceHandlerLocal->lies_onInterface_quad(dof2) && !m_spInterfaceHandlerLocal->check_vertex_modus(ON_INTERFACE, dof2, -1) )
 				indexCol = DoFIndex(numAllDoFs + dof2_real, 0);
			else
				indexCol = DoFIndex(colInd.index(0, dof2_real), colInd.comp(0, dof2_real));

		// finally add loc to glob:
			DoFRef(mat, indexRow, indexCol) += locJ_quad.value(0, dof1, 0, dof2);

			if ( indexRow[0] == 0 || indexRow[0] == 30){
				UG_LOG("------------------------ quad: += " << locJ_quad.value(0, dof1, 0, dof2) << "\n");
				UG_LOG("(indexRow, indexCol) = " << indexRow[0] << "," << indexCol[0] << "\n");
			}
		}
	}

	//UG_LOG("------------------------ END DiffusionInterfaceMapper ------------------------\n");

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
	const LocalVector& locD_tri = get_local_defect_tri();
	const bool shift_global_index_tri = m_spInterfaceHandlerLocal->get_index_shift_tri();

	if ( lvec.num_all_dof(0) != locD_tri.num_all_dof(0) )
		UG_THROW("in 'add_local_vec_to_global_Nitsche': non-consistent sizees!\n");

 	size_t numAllDoFs = m_numDoFs;

 	if ( shift_global_index_tri )
 	{		numAllDoFs = m_numDoFs + m_numNewDoFs;
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

	bool print = false;
	DoFIndex index_print;

	const LocalIndices& ind = lvec.get_indices();
	DoFIndex index;

	///////////////////////////////////////////////////////////////
	/// FIRST: add loc to glob for locD_tri:
	///////////////////////////////////////////////////////////////
	const LocalVector& locD_tri = get_local_defect_tri();
	const bool shift_global_index_tri = m_spInterfaceHandlerLocal->get_index_shift_tri();

 	size_t numAllDoFs = m_numDoFs;

 	if ( shift_global_index_tri )
 	{		numAllDoFs = m_numDoFs + m_numNewDoFs;
 			UG_LOG("it is shifted!\n");
 	}
//	UG_LOG("in DiffusionInterfaceMapper::add_loc_vec(): locD_tri = " << locD_tri << "\n");

	for (size_t dof = 0; dof < locD_tri.num_all_dof(0); ++dof)
	{
		const size_t dof_real = m_spInterfaceHandlerLocal->real_index_tri(dof);

	 	if ( (numAllDoFs + dof_real) > 352 )
	 	{	 		print = true;index_print = DoFIndex(numAllDoFs + dof_real,0);}

	// if dof_real is index of m_vertex: compute global index:
        if ( m_spInterfaceHandlerLocal->lies_onInterface_tri(dof) && !m_spInterfaceHandlerLocal->check_vertex_modus(ON_INTERFACE, dof, 1) )
            index = DoFIndex(numAllDoFs + dof_real,0);
		else
			index = DoFIndex(ind.index(0, dof_real), ind.comp(0, dof_real));

	// finally add loc to glob:
		DoFRef(vec, index) += locD_tri.value(0, dof);

	}


	///////////////////////////////////////////////////////////////
	/// SECOND: add loc to glob for locU_quad:
	///////////////////////////////////////////////////////////////
	const LocalVector& locD_quad = get_local_defect_quad();
	const bool shift_global_index_quad = m_spInterfaceHandlerLocal->get_index_shift_quad();

// reset numAllDoFs!
    numAllDoFs = m_numDoFs;
    
 	if ( shift_global_index_quad )
 		numAllDoFs = m_numDoFs + m_numNewDoFs;

//	UG_LOG("in DiffusionInterfaceMapper::add_loc_vec(): locD_quad = " << locD_quad << "\n");

	for (size_t dof = 0; dof < locD_quad.num_all_dof(0); ++dof)
	{
		size_t dof_real = m_spInterfaceHandlerLocal->real_index_quad(dof);

	 	if ( (numAllDoFs + dof_real) > 352 )
	 	{	 		print = true;index_print = DoFIndex(numAllDoFs + dof_real,0);}

	// if dof_real is index of m_vertex: compute global index:
        if ( m_spInterfaceHandlerLocal->lies_onInterface_quad(dof) && !m_spInterfaceHandlerLocal->check_vertex_modus(ON_INTERFACE, dof, 1) )
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
