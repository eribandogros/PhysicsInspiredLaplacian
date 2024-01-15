//
//  grid_mesh.hpp
//  
//
//  Created by Emily Ribando-Gros
//

#ifndef grid_mesh_hpp
#define grid_mesh_hpp

#include <stdio.h>
#include <iostream>

#include "SparseMatrix.h"
#include "Vector3d.h"

#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

class grid_mesh{
public:
    double m_grid_length;
    int m_num_dimensions_x, m_num_dimensions_y, m_num_dimensions_z;
    double m_min_x, m_min_y, m_min_z, m_max_x, m_max_y, m_max_z;
    
    int m_num_grid_points;
    int m_estimated_num_edges;
    int m_estimated_num_faces;
    int m_estimated_num_cells;
    
    CSparseMatrix *m_d0;
    CSparseMatrix *m_d1;
    CSparseMatrix *projection1;
    
    bool *m_is_edge_counted, *m_is_counted_face;
    bool *m_is_cell_included;
    bool *m_is_boundary_edge, *m_is_boundary_vertex;
    void decode_face_index(int face_index,int &x,int &y,int &z,int &type);
    void check_cell_info();
    bool is_face_included(int i,int j,int k,int type,bool &is_interior, int &num_np);
    bool *m_is_interior_face; //indicate whether this face is an interior face or not
    double *m_faces_size;
    
    CVector3d return_grid_pos(int i,int j,int k,double grid_length,CVector3d start_pos);
    
    int return_point_index(int index_x,int index_y,int index_z);
    int return_edge_index(int index_x,int index_y,int index_z);
    int return_face_index(int index_x,int index_y,int index_z);
    int return_face_index_with_type(int index_x,int index_y,int index_z,int type);
    int return_cell_index(int index_x,int index_y,int index_z);
    
    void tag_faces_edges(int c_i,int c_j,int c_k);
    bool is_edge_counted(int i,int j,int k, int type,int *coeff,
                         int *adj_face_index,int &p1_index,int &p2_index);
    
    void init();
    void initial_mesh();
    void initial_mesh_cell();
    void initial_mesh_face();
    void initial_mesh_edge_vertex();
    
};

#endif /* grid_mesh_hpp */
