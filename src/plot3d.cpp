#include "plot3d.hpp"

using namespace std;

PLOT3D::PLOT3D(){
    flip_normal = false;
}

PLOT3D::~PLOT3D(){

}

void PLOT3D :: set_surface_filename(std::string name){
    surface_filename = name;
}

void PLOT3D :: set_surface(std::shared_ptr<Surface> surf){
    surface = surf;
}

void PLOT3D :: flip_normals(bool val){
    flip_normal = val;
}

void PLOT3D :: read_surface(std::string name){

    set_surface_filename(name);

// open mesh file to read
    ifstream mesh(surface_filename);

// assert if mesh file is open
    assert(mesh.is_open());

// read number of blocks
    mesh >> blocks;

    assert(blocks = 1);

    int KMAX;
// read number of nodes in each direction
    mesh >> IMAX >> JMAX >> KMAX;

//assert number of nodes in z direction is equal to 1
    assert(KMAX = 1);

    cout << "Mesh " << surface->get_name_surface() <<" has " << (IMAX-1)*(JMAX-1) << " panels." << endl;
    
    surface->set_IMAX_JMAX(IMAX, JMAX);
    surface->nodes.resize(IMAX*JMAX);

//read_coordinates
    int index ;
    for(int dim = 0; dim < 3; dim++){
        for(int i = 0; i < IMAX; i++){
            for(int j = 0; j < JMAX; j++){
                index = i*JMAX + j;
                // convert 2d index to 1d
                mesh >> surface->nodes[index][dim];

            } /* end j loop */
        } /* end i loop */
    } /* end dim loop */

//populate panels with its node numbers
    for(int j = 0; j < JMAX - 1; j++){
        for(int i = 0; i < IMAX - 1; i++){

            // plot3d being structured grid, all panels will have 4 nodes
            vector<int> new_panel(4);

            if(flip_normal){
                new_panel[0] = j*IMAX + i ;
                new_panel[1] = j*IMAX + (i+1) ;
                new_panel[2] = (j+1)*IMAX + (i+1) ;
                new_panel[3] = (j+1)*IMAX + i ;
            }else{
                new_panel[0] = j*IMAX + (i+1) ;
                new_panel[1] = j*IMAX + i ;
                new_panel[2] = (j+1)*IMAX + i ;
                new_panel[3] = (j+1)*IMAX + (i+1) ;
            }

            //cout << new_panel[0] << "\t" << new_panel[1] << "\t" << new_panel[2] << "\t" << new_panel[3] << endl;

            surface->panels.push_back(new_panel);

        } /* end i loop */
    }/* end j loop */

// build topology
    build_topology();
}

void PLOT3D :: build_topology(){

    //compute neighbours
    int imax = IMAX - 1;
    int jmax = JMAX - 1;

    for(int j = 0; j < jmax; j++){
        for(int i = 0; i < imax; i++){

            vector<int> new_neighbours;
            vector<int> new_neighbours2;

            if(i == 0){
                new_neighbours.push_back(j * imax + (i+1));
                
                new_neighbours2.push_back(j * imax + (i+1)); // XD added
                new_neighbours2.push_back(j * imax + (i+2));
            }else if(i == imax - 1){
                new_neighbours.push_back(j * imax + (i-1));
                
                new_neighbours2.push_back(j * imax + (i-1)); // XD added
                new_neighbours2.push_back(j * imax + (i-2));
            }else{
                new_neighbours.push_back(j * imax + (i+1));
                new_neighbours.push_back(j * imax + (i-1));
                
                new_neighbours2.push_back(j * imax + (i+1)); // XD added
                new_neighbours2.push_back(j * imax + (i-1));
            }

            if(j == 0){
                new_neighbours.push_back((j+1) * imax + i);
                
                new_neighbours2.push_back((j+1) * imax + i); // XD added
                new_neighbours2.push_back((j+2) * imax + i);
            }else if(j == jmax - 1){
                new_neighbours.push_back((j-1) * imax + i);
                
                new_neighbours2.push_back((j-1) * imax + i); // XD added
                new_neighbours2.push_back((j-2) * imax + i);
            }else{
                new_neighbours.push_back((j+1) * imax + i);
                new_neighbours.push_back((j-1) * imax + i);
                
                new_neighbours2.push_back((j+1) * imax + i);
                new_neighbours2.push_back((j-1) * imax + i);
            }
//            for (int k=0;k<new_neighbours2.size();k++)
//                cout << new_neighbours2[k] << " " ;
//            cout << endl;

            surface->panel_neighbours.push_back(new_neighbours);
            surface->neighbours_two.push_back(new_neighbours2); // XD added
    // order in neighbours_two: front / back    / right / left
    //      when i=0, TE lower: front / front+1 / right / left
    //   when i=imax, TE upper: back  / back-1  / right / left
    //     when j=0, first row: front / back    / right / right+1
    //   when j=jmax, last row: front / back    / left  / left+1

        } /* end i loop */
    } /* end j loop */

    vector3d beam_node,beam_node_collocation,leading_node,trailing_node;
    int index ;
    //compute leading and trailing edge nodes
    for(int j = 0; j < JMAX; j++){
        index = j*IMAX + (IMAX+1)/2 - 1;
        surface->leading_edge_nodes.push_back(index);
        
        leading_node = surface->nodes[index];
        
        index = j*IMAX ;
        surface->trailing_edge_nodes.push_back(index);
        
        trailing_node = surface->nodes[index];
        
// beam_node is located at the border of the panels
        beam_node = leading_node + (trailing_node-leading_node)*0.25 ;// definition of beam node at chord /4
        surface->beam_nodes.push_back(beam_node);
// beam_node_collocation is located at the middle of the panels
        if (j>0) {
            beam_node_collocation = (surface->beam_nodes[j-1] + surface->beam_nodes[j])*0.5 ;
            surface->beam_nodes_collocation.push_back(beam_node_collocation);
        }
        
    } /* end j loop */

    //compute trailing edge panels
    imax = IMAX - 1;
    jmax = JMAX - 1;
    
    surface->is_TE_panel.clear();
    surface->is_TE_panel.resize(imax*jmax,0);
    for(int j = 0; j < jmax; j++){
        if(flip_normal){
            surface->upper_TE_panels.push_back(j*imax);
            surface->lower_TE_panels.push_back(imax*(j+1) - 1);
            surface->is_TE_panel[j*imax]        =1;   
            surface->is_TE_panel[imax*(j+1) - 1]=-1;
        }
        else {
            surface->lower_TE_panels.push_back(j*imax);
            surface->upper_TE_panels.push_back(imax*(j+1) - 1);
            surface->is_TE_panel[j*imax]        =-1;
            surface->is_TE_panel[imax*(j+1) - 1]=1;
        }
    }

//    for(int l = 0; l < surface->upper_TE_panels.size(); l++){
//        cout << surface->upper_TE_panels[l] << endl;
//    }

}

void PLOT3D :: set_domain(std::shared_ptr<Domain> dom){
    domain = dom;
}

void PLOT3D :: read_domain(std::string name){

    domain_filename = name;

    // open mesh file to read
    ifstream mesh(domain_filename);

    // assert if mesh file is open
    assert(mesh.is_open());

    // read number of blocks
    mesh >> blocks;

    assert(blocks = 1);

    int KMAX;
    // read number of nodes in each direction
    mesh >> IMAX >> JMAX >> KMAX;

    domain->set_domain_ranges(IMAX,JMAX,KMAX);

    domain->nodes.resize(IMAX*JMAX*KMAX);

    //read_coordinates (read such that fastest index is x, then y, then z so that it can produce correct vtk output)
    for(int dim = 0; dim < 3; dim++){
        int index = 0;
        for(int k = 0; k < KMAX; k++){
            for(int j = 0; j < JMAX; j++){
                for(int i = 0; i < IMAX; i++){

                    //cout << index << endl;
                    mesh >> domain->nodes[index][dim];
                    index++;

                } /* end i loop */
            } /* end j loop */
        } /* end k loop */
    } /* end dim loop */

//    for(size_t n = 0; n < domain->nodes.size(); n++)
//        cout << domain->nodes[n] << endl;


}
