//
//  main.cpp
//
//
//  Created by Emily Ribando-Gros
//

#include "main.hpp"
#include "grid_mesh.hpp"
#include "hamiltonian.hpp"

#include <string>
#include <iostream>
#include <cmath>

int main(int argc, char* argv[]){
    
    if(argc < 2){
        std::cout << "Instructions: " << std::endl;
        std::cout << "argv[1]: filename" << std::endl;
        std::cout << "argv[2]: scale option (1->0.6 or 2->0.7)" << std::endl;
        std::cout << "argv[3]: L0 or L1 based Hamiltonian (0 or 1)" << std::endl;
        std::cout << "argv[4]: potential (g->Gaussian or l->Lorentz)" << std::endl;
    }
    
    molecule m;
    m.m_filename = argv[1]; //Example: ../dataset/xyz-folder-train/train_0.npy
    
    m.read_molecule();
    m.init();
    
    //assumes min and maxes are always - and + resp.
    if (m.m_max_x < 0 || m.m_max_y < 0 || m.m_max_z < 0 || m.m_min_x > 0
        || m.m_min_y > 0 || m.m_min_z > 0){
        std::cout << "ERROR - bounding box" << std::endl;
    }
    
    std::cout << "input file: " << m.m_filename << std::endl;

    grid_mesh g;
    g.m_grid_length = 0.7;
    
    g.m_max_x = ceil(m.m_max_x)+1; g.m_max_y = ceil(m.m_max_y)+1; g.m_max_z = ceil(m.m_max_z)+1;
    g.m_min_x = floor(m.m_min_x)-1; g.m_min_y = floor(m.m_min_y)-1; g.m_min_z = floor(m.m_min_z)-1;
    
    double x_size = g.m_max_x+abs(g.m_min_x);
    g.m_num_dimensions_x = ceil(x_size/g.m_grid_length)+1;
    double y_size = g.m_max_y+abs(g.m_min_y);
    g.m_num_dimensions_y = ceil(y_size/g.m_grid_length)+1;
    double z_size = g.m_max_z+abs(g.m_min_z);
    g.m_num_dimensions_z = ceil(z_size/g.m_grid_length)+1;

    std::cout << "grid dimensions: " << g.m_num_dimensions_z << " " << g.m_num_dimensions_y
        << " " << g.m_num_dimensions_x << std::endl;

    g.init();
    
    //for element specific file naming
//    char delim = '/';
//    std::string trim = m.m_filename;
//    trim.erase(0, trim.find(delim)+1);
//    trim.erase(0, trim.find(delim)+1);
//    trim.erase(0, trim.find(delim)+1);
//    trim.erase(trim.find(delim), trim.length()-1);
//    trim.erase(0, trim.find_last_of('-')+1);
//    m.elements = trim;

    hamiltonian h(g, m);
    h.m_num_eigenvalues = 20;
    
    h.m_eta_option = atoi(argv[2]);
    h.order = atoi(argv[3]);
    h.m_potential_option = *argv[4];
    
    h.calculate_potential();
    
    h.matlabEIGS();
    
    return 0;
}
