//
//  hamiltonian.cpp
//  
//
//  Created by Emily Ribando-Gros
//

#include "hamiltonian.hpp"
#include "cnpy.h"

#define PI 3.1415926535897932384

void molecule::read_molecule(){
    //extract file extension
    std::string ext = m_filename.substr(m_filename.find_last_of(".")+1);
    
    if(ext == "xyz")
        read_xyz();
    else if(ext == "npy")
        read_npy();
}

void molecule::read_xyz(){
    std::string f = m_filename;
    std::ifstream file(f);

    std::string line;
    
    while(std::getline(file, line)){
        std::stringstream ss(line);
        double x, y, z;
        ss>>x>>y>>z;
        m_atoms.push_back(CVector3d(x,y,z));
    }
    file.close();
    m_numAtoms = m_atoms.size();
}

void molecule::read_npy(){
    cnpy::NpyArray arr = cnpy::npy_load(m_filename);
    double* loaded_data = arr.data<double>();
    
    for(int i = 0; i < arr.shape[0];){
        double x = loaded_data[i++];
        double y = loaded_data[i++];
        double z = loaded_data[i++];
        m_atoms.push_back(CVector3d(x,y,z));
    }
    m_numAtoms = m_atoms.size();
    
    //output the xyz
//    std::string trim_name = m_filename.substr(m_filename.find_first_of("/")+1);
//    trim_name = trim_name.substr(0, trim_name.find_last_of("."));
//    trim_name += ".xyz";
//    std::ofstream xyzout(trim_name);
//    for(int i = 0; i<m_numAtoms; i++){
//        xyzout << m_atoms[i].a.x() << " " << m_atoms[i].a.y() << " " << m_atoms[i].a.z() << std::endl;
//    }
//    xyzout.close();
}

void molecule::init(){
    
    for(int i = 0; i < m_numAtoms; i++){
        //min_x, min_y, min_z, max_x, max_y, max_z
        double x = m_atoms[i].a.x();
        if(x < m_min_x)
            m_min_x = x;
        if(m_max_x < x)
            m_max_x = x;
        
        double y = m_atoms[i].a.y();
        if(y < m_min_y)
            m_min_y = y;
        if(m_max_y < y)
            m_max_y = y;
        
        double z = m_atoms[i].a.z();
        if(z < m_min_z)
            m_min_z = z;
        if(m_max_z < z)
            m_max_z = z;
    }
}

//Euclidean distance between 2 3D points
double distance(CVector3d p1, CVector3d p2){
    float diff = p1.x()-p2.x();
    float dist = diff*diff;
    diff = p1.y()-p2.y();
    dist += diff*diff;
    diff = p1.z()-p2.z();
    dist += diff*diff;
    return sqrt(dist);
}

double hamiltonian::power(double distance, double rate){
    return pow(distance, rate);
}

double hamiltonian::gaussian(double distance){
    return std::exp(-distance*distance);
}

double hamiltonian::lorentz(double distance){
    return 1./(1.+distance*distance);
}

double hamiltonian::potential(double distance){
    if(m_potential_option == 'g')
        return -gaussian(distance);
    else if (m_potential_option == 'q')
        return power(distance, 2);
    else if (m_potential_option == 'c')
        return power(distance, 3);
    else if (m_potential_option == 'l')
        return lorentz(distance);
}

void hamiltonian::calculate_potential(){
    CVector3d start_pos = CVector3d(m_g.m_min_x, m_g.m_min_y, m_g.m_min_z);
    set_eta();
    m_E = std::vector<double>(m_g.m_estimated_num_edges, 0);
    int count = 0;
    int index = 0;
    for(int k=0;k<m_g.m_num_dimensions_z;k++)
            for(int j=0;j<m_g.m_num_dimensions_y;j++)
                for(int i=0;i<m_g.m_num_dimensions_x;i++){
                    
                    CVector3d v = m_g.return_grid_pos(i,j,k,m_g.m_grid_length,start_pos);
                    //V0
                    double Vsum;
                    Vsum = 0;
                    for(molecule::atom a : m_m.m_atoms){
                        Vsum += potential(distance(v, a.a)/m_eta); //0-form potential
                    }
                    m_V.push_back(Vsum);
                }
    
    if (order == 1){//1-form potential
        for(int k=0;k<m_g.m_num_dimensions_z;k++)
            for(int j=0;j<m_g.m_num_dimensions_y;j++)
                for(int i=0;i<m_g.m_num_dimensions_x;i++){
                    int general_edge_id = m_g.return_edge_index(i,j,k);
                    int coeff[4],adj_face_index[4];
                    int p1_index,p2_index;
                    for(int type=0;type<3;type++){
                        bool tmp = m_g.is_edge_counted(i,j,k,type,coeff,adj_face_index,p1_index,p2_index);
                        
                        if(general_edge_id < m_g.m_estimated_num_edges){
                            m_E[general_edge_id+type] = .5*(m_V[p1_index]+m_V[p2_index]);//.5* //add some HS constant
                        }
                    }
                }
    }
}

void hamiltonian::set_eta(){
    
    if(m_potential_option == 'g' || m_potential_option == 'l'){
        if (order == 1){
            switch(m_eta_option){
                case 1:
                    m_eta = 0.6;
                    break;
                case 2:
                    m_eta = 0.5;
                    break;
                case 3:
                    m_eta = 0.8;
            }
        }
        else{ //order 0
            switch(m_eta_option){
                case 1:
                    m_eta = 0.6;
                    break;
                case 2:
                    m_eta = 0.7;
                    break;
                case 3:
                    m_eta = 0.8;
            }
        }
    }
    else if (m_potential_option == 'q' || m_potential_option == 'c'){
        switch(m_eta_option){
            case 1:
                m_eta = 100;
                break;
            case 2:
                m_eta = 200;
        }
    }
    else {
        std::cout << "ERROR - potential" << std::endl;
    }
}

void hamiltonian::matlabEIGS(){
    using namespace matlab::engine;

    std::cout << "Start MATLAB" << std::endl;
    std::unique_ptr<MATLABEngine> mat = startMATLAB();
    matlab::data::ArrayFactory factory;
    matlab::data::ArrayDimensions ad;
    
    std::vector<double> gl = {m_g.m_grid_length};
    matlab::data::TypedArray<double> mgl = factory.createArray<std::vector<double>::iterator>({ 1, 1 }, gl.begin(), gl.end());
    
    mat->setVariable(u"grid_length", std::move(mgl));
    
    std::vector<int> numev = {m_num_eigenvalues};
    matlab::data::TypedArray<int> mnumev = factory.createArray<std::vector<int>::iterator>({ 1, 1 }, numev.begin(), numev.end());
    
    mat->setVariable(u"numev", std::move(mnumev));

    // D0
    int numEld0 = 0;
    CMatrixElement* theElem;
    for(int i = 0; i < m_g.m_d0->numRows; i++){
        for(theElem = m_g.m_d0->rowList[i]; theElem != NULL; theElem =
            theElem->rowNext){
                numEld0++;
        }
    }
    double rowd0[numEld0];
    double cold0[numEld0];
    double vald0[numEld0];
    int count = 0;
    for(int i = 0; i < m_g.m_d0->numRows; i++){
        for(theElem = m_g.m_d0->rowList[i]; theElem != NULL; theElem =
            theElem->rowNext){
            rowd0[count] = theElem->i+1;
            cold0[count] = theElem->j+1;
            vald0[count] = theElem->value;
            count++;
        }
    }
    
    ad.push_back(numEld0);
    ad.push_back(1);
    matlab::data::TypedArray<double> mrowd0 = factory.createArray(ad, rowd0, rowd0+numEld0);
    matlab::data::TypedArray<double> mcold0 = factory.createArray(ad, cold0, cold0+numEld0);
    matlab::data::TypedArray<double> mvald0 = factory.createArray(ad, vald0, vald0+numEld0);
    ad.clear();
    std::vector<double> vsized0 = {static_cast<double>(m_g.m_d0->numRows), static_cast<double>(m_g.m_d0->numCols)};
    matlab::data::TypedArray<double> msize0 = factory.createArray<std::vector<double>::iterator>({1, 2}, vsized0.begin(), vsized0.end());
    
    //D1
    int numEld1 = 0;
    for(int i = 0; i < m_g.m_d1->numRows; i++){
        for(theElem = m_g.m_d1->rowList[i]; theElem != NULL; theElem =
            theElem->rowNext){
                numEld1++;
        }
    }
    double rowd1[numEld1];
    double cold1[numEld1];
    double vald1[numEld1];
    count = 0;
    for(int i = 0; i < m_g.m_d1->numRows; i++){
        for(theElem = m_g.m_d1->rowList[i]; theElem != NULL; theElem = theElem->rowNext){
            rowd1[count] = theElem->i+1;
            cold1[count] = theElem->j+1;
            vald1[count] = theElem->value;
            count++;
        }
    }
    
    ad.push_back(numEld1);
    ad.push_back(1);
    matlab::data::TypedArray<double> mrowd1 = factory.createArray(ad, rowd1, rowd1+numEld1);
    matlab::data::TypedArray<double> mcold1 = factory.createArray(ad, cold1, cold1+numEld1);
    matlab::data::TypedArray<double> mvald1 = factory.createArray(ad, vald1, vald1+numEld1);
    ad.clear();
    std::vector<double> vsized1 = {static_cast<double>(m_g.m_d1->numRows), static_cast<double>(m_g.m_d1->numCols)};
    matlab::data::TypedArray<double> msized1 = factory.createArray<std::vector<double>::iterator>({ 1, 2 }, vsized1.begin(), vsized1.end());
    
    //P1
    int numElp1 = 0;
    for(int i = 0; i < m_g.projection1->numRows; i++){
        for(theElem = m_g.projection1->rowList[i]; theElem != NULL; theElem =
            theElem->rowNext){
            numElp1++;
        }
    }
    double rowp1[numElp1];
    double colp1[numElp1];
    double valp1[numElp1];
    count = 0;
    for(int i = 0; i < m_g.projection1->numRows; i++){
        for(theElem = m_g.projection1->rowList[i]; theElem != NULL; theElem = theElem->rowNext){
            rowp1[count] = theElem->i+1;
            colp1[count] = theElem->j+1;
            valp1[count] = theElem->value;
            count++;
        }
    }
    
    ad.push_back(numElp1);
    ad.push_back(1);
    matlab::data::TypedArray<double> mrowp1 = factory.createArray(ad, rowp1, rowp1+numElp1);
    matlab::data::TypedArray<double> mcolp1 = factory.createArray(ad, colp1, colp1+numElp1);
    matlab::data::TypedArray<double> mvalp1 = factory.createArray(ad, valp1, valp1+numElp1);
    ad.clear();
    std::vector<double> vsizep1 = {static_cast<double>(m_g.projection1->numRows), static_cast<double>(m_g.projection1->numCols)};
    matlab::data::TypedArray<double> msizep1 = factory.createArray<std::vector<double>::iterator>({ 1, 2 }, vsizep1.begin(), vsizep1.end());
    
    mat->setVariable(u"row0", std::move(mrowd0));
    mat->setVariable(u"col0", std::move(mcold0));
    mat->setVariable(u"val0", std::move(mvald0));
    mat->setVariable(u"size0", std::move(msize0));
    
    mat->setVariable(u"rowd1", std::move(mrowd1));
    mat->setVariable(u"cold1", std::move(mcold1));
    mat->setVariable(u"vald1", std::move(mvald1));
    mat->setVariable(u"sized1", std::move(msized1));
    
    mat->setVariable(u"rowp1", std::move(mrowp1));
    mat->setVariable(u"colp1", std::move(mcolp1));
    mat->setVariable(u"valp1", std::move(mvalp1));
    mat->setVariable(u"sizep1", std::move(msizep1));
    
    ad.push_back(m_g.m_num_grid_points);
    ad.push_back(1);
    matlab::data::TypedArray<double> mV = factory.createArray<double>(ad, m_V.data(), m_V.data()+m_V.size());
    ad.clear();
    
    mat->setVariable(u"V", std::move(mV));
    
    ad.push_back(m_g.m_estimated_num_edges);
    ad.push_back(1);
    matlab::data::TypedArray<double> mE = factory.createArray<double>(ad, m_E.data(), m_E.data()+m_E.size());
    ad.clear();
    
    mat->setVariable(u"E", std::move(mE));

    //L0
    mat->eval(u"D0=sparse(row0,col0,val0,size0(1),size0(2));");
    mat->eval(u"S1=speye(size0(1))*grid_length;");
    mat->eval(u"L0=D0'*S1*D0;");
    mat->eval(u"P0=spdiags(V,0,size0(2),size0(2));"); //potential
    
    if (order == 0){
        mat->eval(u"H=L0+P0;");
        mat->eval(u"[V,D]=eigs(H*1/(grid_length*grid_length*grid_length),speye(size0(2)),numev,'smallestreal');");
    }
    
    //L1
    else if (order == 1){
        
        mat->eval(u"S0_inv=speye(size0(2))*1/(grid_length*grid_length*grid_length);");
        mat->eval(u"S2=speye(sized1(1))*1/grid_length;");
        mat->eval(u"P1=sparse(rowp1,colp1,valp1,sizep1(1),sizep1(2));"); //projection
        mat->eval(u"D1=sparse(rowd1,cold1,vald1,sized1(1),sized1(2));");
        
        mat->eval(u"S1P=P1*S1*P1';");
        mat->eval(u"D0P=P1*D0;");
        mat->eval(u"D1P=D1*P1';");

        mat->eval(u"L1_right=S1P*D0P*S0_inv*D0P'*S1P;"); // s1 D0 s0^-1 D0^T s1
        mat->eval(u"L1_left=D1P'*S2*D1P;"); // D1^T S2 D1
        mat->eval(u"L1=L1_left+L1_right;"); //Hodge L1n
        mat->eval(u"P=P1*(spdiags(E,0,sized1(2),sized1(2)))*P1';"); //potential
        mat->eval(u"H=L1+P;");//+speye(sized1(2))*1e-10 //(H+H')/2

        mat->eval(u"H=(H+H')/2.;");
        mat->eval(u"H=H+speye(size(H,1))*1e-10;");
        mat->eval(u"[V,D]=eigs(H,P1*S1*P1'+speye(size(S1,1))*1e-10,numev,'sr');");
        //mat->eval(u"diag(D)");

    }

    //write eigenvalues to file
    matlab::data::TypedArray<double> d = mat->getVariable(u"D");
    std::vector<double> eval;
    eval.reserve(d.getDimensions()[0]);
    for(int i=0; i<d.getDimensions()[0]; i++){
        eval.push_back(d[i][i]);
        //m_features.push_back(d[i][i]);
    }
    
    std::string trim_name = m_m.m_filename.substr(m_m.m_filename.find_last_of("/")+1);
    trim_name = trim_name.substr(0, trim_name.find_last_of("."));
    // other naming options...
//    std::string fname = "features/features_" + trim_name + "_" +  std::to_string(m_eta_option) + ".txt";
//    std::string fname = "features/L"+ std::to_string(order) + "/";
    std::string fname = "../features";
//    fname += m_dataset;
//    fname += "/" + std::to_string(m_eta_option) + "/" + trim_name + ".txt";
//    fname += "/" + m_m.elements + "/" + trim_name + ".txt";
    fname += "/" + trim_name + ".txt";

    std::ofstream eigval(fname);
    for(int i=0; i<d.getDimensions()[0]; i++)
        eigval << eval[i] << std::endl;
    eigval.close();
    

    //write eigenvectors to file
//    matlab::data::TypedArray<double> v = mat->getVariable(u"V");
//    std::vector<double> evec;
//
//    evec.reserve(v.getDimensions()[0]);
//    for(int j=0; j<; j++)
//        for(int i=0; i<v.getDimensions()[j]; i++)
//            evec.push_back(v[i][j]);

//    for(int j=0; j<10; j++){
//        std::string fnamevec = "evec_" + trim_name + std::to_string(j) + ".txt";
//        std::ofstream eigvec(fnamevec);
//        for(int i=0; i<v.getDimensions()[j]; i++)
//            eigvec << v[i][j] << std::endl;
//        eigvec.close();
//    }
//
    //write grid point positions
//    CVector3d start_pos = CVector3d(m_g.m_min_x, m_g.m_min_y, m_g.m_min_z);
//    std::ofstream file("grid_points");
//    for(int k=0;k<m_g.m_num_dimensions_z;k++)
//            for(int j=0;j<m_g.m_num_dimensions_y;j++)
//                for(int i=0;i<m_g.m_num_dimensions_x;i++){
//                    CVector3d grid_pos = m_g.return_grid_pos(i,j,k,m_g.m_grid_length,start_pos);
//
//                    file << grid_pos.x() << " " << grid_pos.y() << " " << grid_pos.z() << std::endl;
//                }
//    file.close();

//    write_vtk_vf(5); // the eigenvector field associated with the 5th eigenvalue

}

void hamiltonian::write_vtk_vf(int ev){
    
    std::ofstream out("mesh_vf"+std::to_string(ev)+".vtk");

    out<<"# vtk DataFile Version 2.0"<<std::endl;
    out<<"Structured field"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET STRUCTURED_GRID"<<std::endl;
    out<<"DIMENSIONS " << m_g.m_num_dimensions_x << " " << m_g.m_num_dimensions_y << " " << m_g.m_num_dimensions_z <<std::endl;
    out<<"POINTS " << m_g.m_num_grid_points << " double" <<std::endl;
    
    std::string filename = "grid_points";

    std::ifstream file(filename);
    for(int i=0; i<m_g.m_num_grid_points;i++){
        double x, y, z;
        file >> x >> y >> z;
        out << x << " " << y << " " << z << std::endl;
    }
    file.close();
    
    std::string trim_name = m_m.m_filename.substr(m_m.m_filename.find_first_of("/")+1);
    trim_name = trim_name.substr(0, trim_name.find_last_of("."));
    
    std::string fnamevec0 = "evec_" + trim_name + std::to_string(ev) + ".txt";
    std::ifstream file0(fnamevec0);
    out<<"POINT_DATA "<<m_g.m_num_grid_points<<std::endl;
    out<<"SCALARS U" << ev << " double 1"<<std::endl;
    out<<"LOOKUP_TABLE default"<<std::endl;
    for(int i=0; i<m_g.m_num_grid_points;i++){
        double v;
        file0 >> v;
        out << std::fixed << v << std::endl;
    }
    file0.close();

    out.close();
}

