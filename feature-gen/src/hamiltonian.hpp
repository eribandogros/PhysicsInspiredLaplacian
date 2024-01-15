//
//  hamiltonian.hpp
//  
//
//  Created by Emily Ribando-Gros
//

#ifndef hamiltonian_hpp
#define hamiltonian_hpp

#include <stdio.h>
#include <fstream>
#include <string>
#include <vector>

#include "grid_mesh.hpp"

class molecule{
public:
    std::string m_filename;
    std::string elements;
    int m_numAtoms;
    double m_avgdist = 0;
    double m_mindist = 5;

    double m_min_x = 0, m_min_y = 0, m_min_z = 0, m_max_x = 0, m_max_y = 0, m_max_z = 0;
    
    molecule(std::string f) : m_filename{f} {};
    molecule() = default;
    
    struct atom{
        CVector3d a;
        
        atom(CVector3d A) : a{A} {};
    };
    std::vector<atom> m_atoms;
    
    void read_molecule();
    void read_xyz();
    void read_npy();
    void simple_example();
    void init();
    
};

class hamiltonian{
public:
    grid_mesh m_g;
    molecule m_m;
    
    std::vector<double> m_V;
    std::vector<double> m_E;
    std::vector<double> m_features;
    int m_center_index = 0;
    int m_decay_rate = 2;
    
    double m_eta;
    int m_eta_option = 1;
    char m_potential_option = 'g';
    char m_dataset = 'c';
    int m_num_eigenvalues = 100;
    int order = 0;
    
    hamiltonian(grid_mesh G, molecule M) : m_g{G}, m_m{M} {};
    
    void calculate_potential();
    void set_eta();
    
    void matlabEIGS();
    void write_npy();
    void write_vtk_vf(int ev);
    
    //potentials
    double potential(double distance); //generic
    double power(double distance, double rate);
    double gaussian(double distance);
    double lorentz(double distance);
};

double distance(CVector3d p1, CVector3d p2);

#endif /* hamiltonian_hpp */
