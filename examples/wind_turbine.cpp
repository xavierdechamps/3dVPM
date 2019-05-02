#include <iostream>

#include "plot3d.hpp"
#include "vtk_writer.hpp"
#include "solver.hpp"
#include "wake.hpp"
#include "domain.hpp"

#include <vector>
using namespace std;

// turn off assertions
#define NDEBUG

int main(int argc, char** args)
{
    Parameters::unsteady_problem = true;
//    Parameters :: switch_off_root_vortex = true;
//    Parameters :: switch_off_tip_vortex = true;
    Parameters :: trailing_edge_wake_shed_factor = 0.25;
    Parameters::use_vortex_core_model=true;
    
    double time_step = 0.015;     // seconds
    double fluid_density = 1.225; // kg / m3
//    std::string filename = "NREL_5MW_59x151_3dVPMbis.txt";
    std::string filename = "blade_4.775_modified.x";
//    std::string filename = "NREL_Phase_VI_Blade_101_81.txt";
    int num_blades = 2;
    int num_time_steps = 100 ;
	
// read mesh file
    PLOT3D mesh;
	
// Initialize a vector of pointers to the different wet surfaces of the wind turbine (blade, hub)
	vector<shared_ptr<Surface>> surface ;
    double shift_value = 360. / num_blades ;
    vector3d angular_shifting_blade(0,shift_value,0); angular_shifting_blade = angular_shifting_blade * M_PI / 180.0 ;
    
    for (int i=0; i<num_blades; i++) {
        surface.push_back( shared_ptr<Surface>(new Surface()) );
// Set the names of the surfaces
        surface[i]->set_name_surface("blade_"+std::to_string(i+1));
// Associate a surface and the mesh reader, then read the mesh for each surface of the vector.
        mesh.set_surface(surface[i]);
        mesh.read_surface(filename);
// Rotate blades to build the whole wind turbine
        if (i>0) 
            surface[i]->rotate_surface(angular_shifting_blade*i,true);
	}
		
// set free stream and angular velocities NREL Phase VI
// 71.63RPM = 1.193833333 rotations per second
    vector3d free_stream_velocity(0,7,0); // free stream velocity
    vector3d surface_angular_velocity(0,71.63,0); // angular velocity
    
// set free stream and angular velocities NREL 5MW
// 9.6RPM = 0.16 rotations per second
//    vector3d free_stream_velocity(0,8,0); // free stream velocity
//    vector3d surface_angular_velocity(0,9.6,0); // angular velocity

// Initialize a vector of pointers to the different wakes of the wind turbine (blade, hub)
    vector<shared_ptr<Wake>> wake ;

	for (int i=0;i<surface.size();i++) {
		surface[i]->set_angular_velocity(surface_angular_velocity,false); 
		surface[i]->compute_panel_components();
		
	    wake.push_back( shared_ptr<Wake>(new Wake()) ); // add a new wake to the vector
		wake[i]->set_name_surface("wake_"+std::to_string(i+1));
	    wake[i]->add_lifting_surface(surface[i]); // link this new wake with the corresponding wet surface
		wake[i]->initialize(free_stream_velocity,time_step); // initialize this wake
	}

// Initialize the VTK writer
    shared_ptr<vtk_writer> writer(new vtk_writer());

// Initialize an object of class Solver
	Solver solver(argc,args);
	
	solver.add_surface(surface);
	solver.add_wake(wake);
    solver.add_logger(writer);
    solver.set_free_stream_velocity(free_stream_velocity);
    solver.set_reference_velocity(free_stream_velocity);
    solver.set_fluid_density(fluid_density);

    vector3d angular_displacement = surface_angular_velocity * (2 * M_PI / 60.0) * time_step;

//write coordinates of the collocation points to data variable
    vector<vector<double>> data (4,vector<double>(surface[0]->n_panels()));
    for(int p = 0; p < surface[0]->n_panels(); p++){
        const vector3d& point = surface[0]->get_collocation_point(p);
        data[0][p] = point[0];
        data[1][p] = point[1];
        data[2][p] = point[2];
    }
    vector<double> CP_surf(surface[0]->n_panels());
    std::string file_name;

    for(int i = 0; i < num_time_steps; i++){
        solver.solve(time_step,i);
        solver.convect_wake(time_step, i*time_step);
		for (int j=0;j<surface.size();j++) {
			surface[j]->rotate_surface(angular_displacement,true);
			wake[j]->shed_wake(free_stream_velocity,time_step);
		}
		solver.finalize_iteration(i);
		
        solver.write_output(i);
        solver.write_global_forces(i, time_step);
        
// test motion of beam
//        solver.move_beam_test(time_step,i);
        
        
        int IMAX,JMAX ;
        surface[0]->get_IMAX_JMAX(IMAX, JMAX) ;
        IMAX-- ; JMAX--; // IMAX/JMAX = number of panels along chord and span directions
        for(int p = 0; p < surface[0]->n_panels(); p++){
            data[3][p] = solver.get_pressure_coefficient(0,p);
        }
        
        file_name = "Output/pressure_data_iteration_"+std::to_string(i)+".txt";
        ofstream file(file_name);
        file << "x,y,z,cp" << endl;
        
        for (int jj=0; jj<JMAX; jj++) {
            for (int ii=0; ii<IMAX; ii++) {
                int p = jj * (IMAX) + ii ;
                file << scientific << data[0][p] << " " << data[1][p] << " " << data[2][p] << " " << data[3][p] << endl;
            }
            file << endl << endl ;
        }
        
        /*
        for(int p = 0; p < surface[0]->n_panels(); p++){
            file << scientific << data[0][p] << " " << data[1][p] << " " << data[2][p] << " " << data[3][p] << endl;
        }*/
        file.close();
        
    }
    
    cout << "Exiting program." << endl;

    return 0;
}
