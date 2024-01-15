//
//  grid_mesh.cpp
//  
//
//  Created by Emily Ribando-Gros
//

#include "grid_mesh.hpp"

CVector3d grid_mesh::return_grid_pos(int i,int j,int k,double grid_length,CVector3d start_pos){
    CVector3d diff = CVector3d(i*grid_length,j*grid_length,k*grid_length);
    return (diff+start_pos);
}

int grid_mesh::return_point_index(int index_x,int index_y,int index_z){
    //out of range
    if(index_x<0 || index_y<0 || index_z<0)
        return -1;
    if(index_x>=m_num_dimensions_x || index_y>=m_num_dimensions_y || index_z>=m_num_dimensions_z)
        return -1;
    return (index_z*m_num_dimensions_y*m_num_dimensions_x + index_y*m_num_dimensions_x + index_x);
}

int grid_mesh::return_edge_index(int index_x,int index_y,int index_z){
    return 3*(index_z*m_num_dimensions_y*m_num_dimensions_x +
        index_y*m_num_dimensions_x + index_x);
}

int grid_mesh::return_face_index(int index_x,int index_y,int index_z){
    return 3*(index_z*m_num_dimensions_y*m_num_dimensions_x +
        index_y*m_num_dimensions_x + index_x);
}

int grid_mesh::return_face_index_with_type(int index_x,int index_y,int index_z,int type){
    if(index_x<0 || index_y<0 || index_z<0)
        return -1;
    if(index_x>=m_num_dimensions_x || index_y>=m_num_dimensions_y || index_z>=m_num_dimensions_z)
        return -1;

    //return -1 if dangling face
    if (type == 0 && index_x==m_num_dimensions_x-1)
        return -1;
    else if (type == 1 && index_y==m_num_dimensions_y-1)
        return -1;
    else if (type == 2 && index_z==m_num_dimensions_z-1)
        return -1;
    
    return (return_face_index(index_x,index_y,index_z) + type);
}

void grid_mesh::decode_face_index(int face_index,int &x,int &y,int &z,int &type){
    type = face_index%3;

    int tmp1 = (face_index-type)/3;
    z = tmp1/(m_num_dimensions_x*m_num_dimensions_y);

    int tmp2 = tmp1 - z*m_num_dimensions_x*m_num_dimensions_y;
    y = tmp2/m_num_dimensions_x;

    x = tmp2%m_num_dimensions_x;
}

int grid_mesh::return_cell_index(int index_x,int index_y,int index_z){
    if(index_x<0 || index_x>=m_num_dimensions_x-1)
        return -1;
    if(index_y<0 || index_y>=m_num_dimensions_y-1)
        return -1;
    if(index_z<0 || index_z>=m_num_dimensions_z-1)
        return -1;
    return (index_z*(m_num_dimensions_y-1)*(m_num_dimensions_x-1) +
        index_y*(m_num_dimensions_x-1) + index_x);
}

void grid_mesh::tag_faces_edges(int c_i,int c_j,int c_k){
    
    int adj_face_index[6];
    adj_face_index[0] = return_face_index_with_type(c_i,  c_j,  c_k,  0);
    adj_face_index[1] = return_face_index_with_type(c_i,  c_j,  c_k+1,0);
    adj_face_index[2] = return_face_index_with_type(c_i,  c_j,  c_k,  1);
    adj_face_index[3] = return_face_index_with_type(c_i,  c_j+1,c_k,  1);
    adj_face_index[4] = return_face_index_with_type(c_i,  c_j,  c_k,  2);
    adj_face_index[5] = return_face_index_with_type(c_i+1,c_j,  c_k,  2);

    int adj_edge_index[12];
    adj_edge_index[0] = return_edge_index(c_i,  c_j,  c_k  );
    adj_edge_index[1] = return_edge_index(c_i+1,c_j,  c_k  )+2;
    adj_edge_index[2] = return_edge_index(c_i,  c_j,  c_k+1);
    adj_edge_index[3] = return_edge_index(c_i,  c_j,  c_k  )+2;

    adj_edge_index[4] = return_edge_index(c_i,  c_j+1,c_k  );
    adj_edge_index[5] = return_edge_index(c_i+1,c_j+1,c_k  )+2;
    adj_edge_index[6] = return_edge_index(c_i,  c_j+1,c_k+1);
    adj_edge_index[7] = return_edge_index(c_i,  c_j+1,c_k  )+2;

    adj_edge_index[8]  = return_edge_index(c_i,  c_j,  c_k  )+1;
    adj_edge_index[9]  = return_edge_index(c_i+1,c_j,  c_k  )+1;
    adj_edge_index[10] = return_edge_index(c_i,  c_j,  c_k+1)+1;
    adj_edge_index[11] = return_edge_index(c_i+1,c_j,  c_k+1)+1;

    for(int i=0;i<6;i++)
        m_is_counted_face[adj_face_index[i]] = true;
}

bool grid_mesh::is_edge_counted(int i,int j,int k,int type,int *coeff,
                                int *adj_face_index,int &p1_index,int &p2_index){
    p1_index = return_point_index(i,j,k);
    
    if(type==0){
        if(i==m_num_dimensions_x-1)
            return false;
        p2_index = return_point_index(i+1,j,k);
    }
    if(type==1){
        if(j==m_num_dimensions_y-1)
            return false;
        p2_index = return_point_index(i,j+1,k);
    }
    if(type==2){
        if(k==m_num_dimensions_z-1)
            return false;
        p2_index = return_point_index(i,j,k+1);
    }

    if(type==0){
        adj_face_index[0] = return_face_index_with_type(i,j,  k,  2); coeff[0] = 1;
        adj_face_index[1] = return_face_index_with_type(i,j-1,k,  2); coeff[1] = -1;
        adj_face_index[2] = return_face_index_with_type(i,j,  k-1,1); coeff[2] = 1;
        adj_face_index[3] = return_face_index_with_type(i,j,  k,  1); coeff[3] = -1;
    }
    if(type==1){
        adj_face_index[0] = return_face_index_with_type(i-1,j,k,  0); coeff[0] = 1;
        adj_face_index[1] = return_face_index_with_type(i,  j,k,  0); coeff[1] = -1;
        adj_face_index[2] = return_face_index_with_type(i,  j,k,  2); coeff[2] = 1;
        adj_face_index[3] = return_face_index_with_type(i,  j,k-1,2); coeff[3] = -1;
    }
    if(type==2){
        adj_face_index[0] = return_face_index_with_type(i-1,j,k,1); coeff[0] = -1;
        adj_face_index[1] = return_face_index_with_type(i,  j,k,1); coeff[1] = 1;
        adj_face_index[2] = return_face_index_with_type(i,  j,k,0); coeff[2] = -1;
        adj_face_index[3] = return_face_index_with_type(i,j-1,k,0); coeff[3] = 1;
    }

    return true;
}

void grid_mesh::initial_mesh_cell(){
    m_is_counted_face = new bool[m_estimated_num_faces];
    m_is_edge_counted = new bool[m_estimated_num_edges];
    
    memset(m_is_counted_face,false,sizeof(bool)*m_estimated_num_faces);
    memset(m_is_edge_counted,false,sizeof(bool)*m_estimated_num_edges);
    
    for(int k=0;k<m_num_dimensions_z-1;k++)
        for(int j=0;j<m_num_dimensions_y-1;j++)
            for(int i=0;i<m_num_dimensions_x-1;i++){
                int cell_index = return_cell_index(i,j,k);
                
                tag_faces_edges(i,j,k);
                
                int adj_face_index[6];
                int coeff[6];
                adj_face_index[0] = return_face_index_with_type(i,  j,  k,  0); coeff[0] = -1;
                adj_face_index[1] = return_face_index_with_type(i,  j,  k+1,0); coeff[1] = 1;
                adj_face_index[2] = return_face_index_with_type(i,  j,  k,  1); coeff[2] = -1;
                adj_face_index[3] = return_face_index_with_type(i,  j+1,k,  1); coeff[3] = 1;
                adj_face_index[4] = return_face_index_with_type(i,  j,  k,  2); coeff[4] = -1;
                adj_face_index[5] = return_face_index_with_type(i+1,j,  k,  2); coeff[5] = 1;

            }
}

//check whether the face existed
bool grid_mesh::is_face_included(int i,int j,int k,int type,bool &is_interior, int &num_np){
    is_interior = false;
    CVector3d origin_pos = CVector3d(0.0,0.0,0.0);
    CVector3d pos[4];
    if(type==0){
        //centered point index (i+0.5, j+0.5, k)
        //the face is out of range
        if(i==m_num_dimensions_x-1 || j==m_num_dimensions_y-1)
            return false;
        pos[0] = return_grid_pos(i,  j,  k,  m_grid_length,origin_pos);
        pos[1] = return_grid_pos(i+1,j,  k,  m_grid_length,origin_pos);
        pos[2] = return_grid_pos(i+1,j+1,k,  m_grid_length,origin_pos);
        pos[3] = return_grid_pos(i,  j+1,k,  m_grid_length,origin_pos);
        
        if(k==0 || k==m_num_dimensions_z-1)
            is_interior = false;
        else{
            if(m_is_cell_included[return_cell_index(i,j,k)] && m_is_cell_included[return_cell_index(i,j,k-1)]){ // two cells must be included
                is_interior = true;
            }
            else
                is_interior = false;
        }
    }

    if(type==1){
        //centered point index (i+0.5, j, k+0.5)
        //the face is out of range
        if(i==m_num_dimensions_x-1 || k==m_num_dimensions_z-1)
            return false;
        
        pos[0] = return_grid_pos(i,  j,  k,  m_grid_length,origin_pos);
        pos[1] = return_grid_pos(i,  j,  k+1,m_grid_length,origin_pos);
        pos[2] = return_grid_pos(i+1,j,  k+1,m_grid_length,origin_pos);
        pos[3] = return_grid_pos(i+1,j,  k,  m_grid_length,origin_pos);

        if(j==0 || j==m_num_dimensions_y-1)
            is_interior = false;
        else{
            if(m_is_cell_included[return_cell_index(i,j,k)] && m_is_cell_included[return_cell_index(i,j-1,k)]){
                is_interior = true;
            }
            else
                is_interior = false;
        }
    }

    if(type==2){
        //centered point index (i, j+0.5, k+0.5)
        //the face is out of range
        if(j==m_num_dimensions_y-1 || k==m_num_dimensions_z-1)
            return false;
        
        pos[0] = return_grid_pos(i,  j,  k,  m_grid_length,origin_pos);
        pos[1] = return_grid_pos(i,  j+1,k,  m_grid_length,origin_pos);
        pos[2] = return_grid_pos(i,  j+1,k+1,m_grid_length,origin_pos);
        pos[3] = return_grid_pos(i,  j,  k+1,m_grid_length,origin_pos);

        if(i==0 || i==m_num_dimensions_x-1)
            is_interior = false;
        else{
            if(m_is_cell_included[return_cell_index(i,j,k)] && m_is_cell_included[return_cell_index(i-1,j,k)]){
                is_interior = true;
            }
            else
                is_interior = false;
        }
    }

    return true;
}

void grid_mesh::initial_mesh_face(){
    m_is_interior_face = new bool[m_estimated_num_faces]; //indicate whether this face is an interior face or not
    
    int m_num_interior_faces = 0;
    m_d1 = new CSparseMatrix(m_estimated_num_faces,m_estimated_num_edges);
    //initialization
    for(int i=0;i<m_estimated_num_faces;i++){
        m_is_interior_face[i] = false;
    }
}

void grid_mesh::initial_mesh_edge_vertex(){
    m_is_boundary_edge = new bool[m_estimated_num_edges];
    m_is_boundary_vertex = new bool[m_num_grid_points];
    m_d0 = new CSparseMatrix(m_estimated_num_edges,m_num_grid_points);
    projection1 = new CSparseMatrix(m_estimated_num_edges,m_estimated_num_edges);
    
    //initialization
    for(int i=0;i<m_estimated_num_edges;i++)
        m_is_boundary_edge[i] = false;
    for(int i=0;i<m_num_grid_points;i++)
        m_is_boundary_vertex[i] = false;
    
    //check for all the interior faces
    for(int i=0;i<m_estimated_num_faces;i++)
       if(m_is_counted_face[i] && !m_is_interior_face[i]){
            int x,y,z,type;
            decode_face_index(i,x,y,z,type);
            if(type==0){
                m_is_boundary_edge[return_edge_index(x,  y,  z  )  ] = true;
                m_is_boundary_edge[return_edge_index(x+1,y,  z  )+1] = true;
                m_is_boundary_edge[return_edge_index(x,  y+1,z  )  ] = true;
                m_is_boundary_edge[return_edge_index(x,  y,  z  )+1] = true;

                m_is_boundary_vertex[return_point_index(x,  y,  z)] = true;
                m_is_boundary_vertex[return_point_index(x+1,y,  z)] = true;
                m_is_boundary_vertex[return_point_index(x+1,y+1,z)] = true;
                m_is_boundary_vertex[return_point_index(x,  y+1,z)] = true;
            }

            if(type==1){
                m_is_boundary_edge[return_edge_index(x,  y,  z  )  ] = true;
                m_is_boundary_edge[return_edge_index(x+1,y,  z  )+2] = true;
                m_is_boundary_edge[return_edge_index(x,  y,  z+1)  ] = true;
                m_is_boundary_edge[return_edge_index(x,  y,  z  )+2] = true;

                m_is_boundary_vertex[return_point_index(x,  y,  z  )] = true;
                m_is_boundary_vertex[return_point_index(x+1,y,  z  )] = true;
                m_is_boundary_vertex[return_point_index(x+1,y,  z+1)] = true;
                m_is_boundary_vertex[return_point_index(x,  y,  z+1)] = true;
            }

            if(type==2){
                m_is_boundary_edge[return_edge_index(x,  y,  z  )+1] = true;
                m_is_boundary_edge[return_edge_index(x,  y+1,z  )+2] = true;
                m_is_boundary_edge[return_edge_index(x,  y,  z+1)+1] = true;
                m_is_boundary_edge[return_edge_index(x,  y,  z  )+2] = true;

                m_is_boundary_vertex[return_point_index(x,  y,  z  )] = true;
                m_is_boundary_vertex[return_point_index(x,  y,  z+1)] = true;
                m_is_boundary_vertex[return_point_index(x,  y+1,z+1)] = true;
                m_is_boundary_vertex[return_point_index(x,  y+1,z  )] = true;
            }
       }

    int count_edge = 0;
    for(int k=0;k<m_num_dimensions_z;k++)
        for(int j=0;j<m_num_dimensions_y;j++)
            for(int i=0;i<m_num_dimensions_x;i++){
                int general_edge_id = return_edge_index(i,j,k);
                //std::cout << general_edge_id << std::endl;
                
                for(int type=0;type<3;type++){
                    int coeff[4],adj_face_index[4];
                    int p1_index,p2_index;
                    
                    bool tmp = is_edge_counted(i,j,k,type,coeff,adj_face_index,p1_index,p2_index);

                    if (tmp){
                        m_is_edge_counted[general_edge_id+type] = true;
                        //std::cout << i << " " << j << " " << k <<std::endl;
                        m_d0->add1Value(general_edge_id+type,p1_index,1);
                        m_d0->add1Value(general_edge_id+type,p2_index,-1);
                        projection1->add1Value(count_edge++,general_edge_id+type,1.0);
                    }
                    else
                        m_is_edge_counted[general_edge_id+type] = false;
                    
                    //avoids dangling edges/faces
                    for(int l=0;l<4;l++)
                        if(adj_face_index[l]!=-1){
                            m_is_interior_face[adj_face_index[l]] = true; //need for matrix cleanup
                            if(m_is_edge_counted[general_edge_id+type] == true)
                                m_d1->add1Value(adj_face_index[l],general_edge_id+type,coeff[l]);
                        }
                }
            }
}

void grid_mesh::initial_mesh(){
    initial_mesh_cell();
    initial_mesh_face();
    initial_mesh_edge_vertex();
}

void grid_mesh::init(){
    
    m_estimated_num_cells = (m_num_dimensions_z-1)*(m_num_dimensions_y-1)*(m_num_dimensions_x-1);
    m_estimated_num_faces = 3*m_num_dimensions_z*m_num_dimensions_y*m_num_dimensions_x;
    m_estimated_num_edges = 3*m_num_dimensions_z*m_num_dimensions_y*m_num_dimensions_x;
    m_num_grid_points = m_num_dimensions_z*m_num_dimensions_y*m_num_dimensions_x;
    
    std::cout << "number of vertices: " << m_num_grid_points << ", edges: " << m_estimated_num_edges << ", faces: " << m_estimated_num_faces << std::endl;

    initial_mesh();

}
