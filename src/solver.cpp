#include "solver.hpp"

using namespace std;

Solver::Solver(const int argC,char** argS) : argc(argC), args(argS) {
    free_stream_velocity = 0;
    reference_velocity = 0;
    density = 0;
    total_n_panels = 0;
}

Solver::~Solver(){    
    if (this->ofile_global_forces.is_open())
        this->ofile_global_forces.close();
    
    if (RHS)                      VecDestroy(&RHS);
    if (solution)                 VecDestroy(&solution);
    if (doublet_influence_matrix) MatDestroy(&doublet_influence_matrix);
    if (ksp_doublet)              KSPDestroy(&ksp_doublet);
    
    if (Jacobian_doublets) MatDestroy(&Jacobian_doublets);
    
// PetscFinalize() causes random crashes, let's simply initialize Petsc once and never finalize
//    PetscFinalize();
}

void Solver :: add_surface(const std::vector<std::shared_ptr<Surface>> &surf) {
    
    this->surface = surf;
    for (int i=0;i<surface.size();i++) {
        assert(surface[i]->n_panels()>0); // the added surface must contain some panels
        total_n_panels += surface[i]->n_panels();
//        cout << "In Solver::add_surface with "<<surface[i]->get_name_surface()<<endl;
    }
}

void Solver :: add_wake(std::vector<std::shared_ptr<Wake>> &surf){
    this->wake = surf;
//    for (int i=0;i<wake.size();i++) 
//        cout << "In Solver::add_wake with "<<wake[i]->get_name_surface()<<endl;
}

void Solver :: add_logger(const std::shared_ptr<vtk_writer> writer){
    log = writer;
}

void Solver :: add_logger(const std::shared_ptr<matlab_writer> writer){
    mlog = writer;
}

void Solver :: set_free_stream_velocity(const vector3d& vel){
    free_stream_velocity = vel;
}

void Solver :: set_reference_velocity(const vector3d& vel){
    reference_velocity = vel;
}

void Solver :: set_fluid_density(const double value){
    density = value;
}

void Solver :: solve(const double dt, int iteration){

    cout << "ITERATION : " << iteration + 1 << endl;
        
    std::vector<int> ind_start_surface(surface.size(),0);
    for (int i=1;i<surface.size();i++)
        ind_start_surface[i] = ind_start_surface[i-1] + surface[i]->n_panels();    

    std::vector<int> ind_start_wake(wake.size(),0);
    for (int i=1;i<wake.size();i++) 
        ind_start_wake[i] = ind_start_wake[i-1] + wake[i]->get_size_doublet_strength(); 
    
// compute source strengths sigma = -u*n
    get_source_array(ind_start_surface);
    
// compute source and doublet coefficients and build influence matrix
    get_source_doublet_matrices(ind_start_surface);

// Compute influence coefficient of the wake doublets
    get_wake_influence_matrix(ind_start_surface,ind_start_wake);
        
// Initialize Petsc variables
    if(iteration == 0)
        initialize_petsc_variables();
    
// Here the A matrix is not modified by the Morino Kutta condition yet (A = doublet coefficient influence matrix)
//     and the matrix W is already built (W = influence coefficient of the wake doublets)
//     and Petsc is initialized
//        -> compute the Jacobian J_mu = - (inv)A * W
    get_Jacobian_wake_on_body(ind_start_surface,ind_start_wake);
    
// Apply Morino linear Kutta condition
    apply_Kutta_Morino(ind_start_surface);
    
// Setup matrices A and B of the system A*X = B
    setup_linear_system();
    
// solve linear system & set unknown doublet strengths
    solve_linear_system();
    
// Get doublets of wake at trailing edge mu_wake = mu_body_upper - mu_body_lower
// Used as a first estimation in the iterative Kutta condition
    get_doublet_wake_TE(ind_start_surface);
    
// Iterate over the doublets to reach zero pressure difference at trailing edge
    iterator_conventional_Kutta(iteration,dt,ind_start_surface) ;
    
// compute surface velocity / surface potential / Pressure coefficient
    get_velocity_potential_pressure(iteration,dt,ind_start_surface);
    
// compute wake strength
    get_new_wake_doublet(iteration,ind_start_surface);

// Compute the global forces and the forces along the beam
    int n_mat = 0;
    for (int i=0;i<surface.size();i++)
        n_mat += surface[i]->beam_nodes_collocation.size() ;
    this->beam_collocation_forces.clear();
    this->beam_collocation_forces.resize(n_mat);
    
    body_forces = 0. ;
    n_mat = 0;
    for (int i=0;i<surface.size();i++){
        body_forces += compute_body_forces(i);

        for (int j=0; j<surface[i]->beam_nodes_collocation.size(); j++) {
            this->beam_collocation_forces[n_mat] = compute_local_beam_forces(i,j) ;
//            cout << surface[i]->beam_nodes_collocation[j] << " " << this->beam_collocation_forces[n_mat] << endl;
            n_mat += 1;
        }
    }


/*
// compute body forces and force coefficients
    body_forces = compute_body_forces();
    body_force_coefficients = compute_body_force_coefficients();
*/
}

double Solver::compute_source_strength(const int surf,const int panel) const {
    const vector3d& node = surface[surf]->get_collocation_point(panel);
    vector3d vel = surface[surf]->get_kinematic_velocity(node) - free_stream_velocity;
    return (vel.dot(surface[surf]->get_panel_normal(panel)));
}

void Solver :: initialize_petsc_variables(){

    assert(total_n_panels > 0);

    // initialize PETSc
    PetscInitialize(&argc,&args,(char*)0,NULL);

    // create PETSc Vec RHS and solution, and PETSc Mat doublet_influence_matrix
    petsc_vec_create(RHS                     ,total_n_panels );
    petsc_vec_create(solution                ,total_n_panels );
    petsc_mat_create(doublet_influence_matrix,total_n_panels,total_n_panels);

    // create KSP solver
    KSPCreate(PETSC_COMM_WORLD,&ksp_doublet);
}

void Solver :: setup_linear_system(){

    // clear previous data, if any.
    VecSet(RHS,0.0);
    VecSet(solution,0.0);
    MatZeroEntries(doublet_influence_matrix);

    int col[total_n_panels];
    
    // setup column number for matrix insertion
    for(int i = 0; i < total_n_panels ; i++)
        col[i] = i;

    // setup RHS = source_coefficient * source_strength
    double *_RHS;
    VecGetArray(RHS,&_RHS);
    
//    PetscLogDouble v1,v2,elapsed_time;
//    PetscTime(&v1);
        
// copy doublet_influence to doublet_influence_matrix //
// WARNING Petsc is not thread-safe -> if there is a run-time error complaining about some corrupted memory, 
//                                     please comment the following line (OpenMP) -> it will be slower but safer
//#pragma omp parallel for
    for(int i = 0;i < total_n_panels; i++){
        double *pt1,*pt2;
        pt1 = &*source_influence[i].begin();
        pt2 = &*source_strength.begin();
        _RHS[i] = cblas_ddot(total_n_panels,pt1,1,pt2,1); // if MKL_CBLAS is not available: comment this line and uncomment the line in the next for loop
//        for(int j = 0; j < total_n_panels; j++){
//            _RHS[i] += source_influence[i][j] * source_strength[j];
//        }
        
        double val[total_n_panels];
        std::copy( doublet_influence[i].begin(), doublet_influence[i].end(), val );
        MatSetValues(doublet_influence_matrix,1,&i,total_n_panels,col,val,INSERT_VALUES); // not thread-safe -> does not work safely with OpenMP
    }
    MatAssemblyBegin(doublet_influence_matrix,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(doublet_influence_matrix,MAT_FINAL_ASSEMBLY);

// RHS -= wake_doublet_coefficient * wake_doublet_strength
    if(wake_doublet_strength.size() > 0){
#pragma omp parallel for
        for(int i = 0; i < total_n_panels; i++){
            for(int j = 0; j < (int)wake_doublet_strength.size(); j++)
                _RHS[i] -= wake_doublet_influence[i][j] * wake_doublet_strength[j] ;
        }
    }

// copy doublet_influence to doublet_influence_matrix (using a faster method to multiply matrix by vector)

/*
    std::vector<double>::const_iterator iit1,iit2,iit3;
    for(int i = 0;i < total_n_panels; i++){
        double val[total_n_panels];
        iit1 = doublet_influence[i].begin();
        iit2 = source_influence[i].begin();
        iit3 = source_strength.begin();
        _RHS[i] = 0. ;
        for(int j = 0; j < total_n_panels; j++){
            val[j] = (*iit1++);
            _RHS[i] += (*iit2++) * (*iit3++);
        }
        MatSetValues(doublet_influence_matrix,1,&i,total_n_panels,col,val,INSERT_VALUES);
    }
    
    MatAssemblyBegin(doublet_influence_matrix,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(doublet_influence_matrix,MAT_FINAL_ASSEMBLY);

// RHS -= wake_doublet_coefficient * wake_doublet_strength (using a faster method to multiply matrix by vector)
    if(wake_doublet_strength.size() > 0){
        double ytemp;
        for(int i = 0; i < total_n_panels; i++){
            iit1 = wake_doublet_influence[i].begin();
            iit2 = wake_doublet_strength.begin();
            ytemp=0;
            for(int j = 0; j < (int)wake_doublet_strength.size(); j++)
                ytemp += (*iit1++) * (*iit2++);
            _RHS[i] -= ytemp;
        }
    
    }
*/    

//    PetscTime(&v2);
//    cout << v2 - v1 << " seconds to build the Petsc matrix and the RHS" << endl;
    
    VecRestoreArray(RHS,&_RHS);
    
}

void Solver :: solve_linear_system(){

    int itn;
    PC pc;

    KSPSetOperators(ksp_doublet,doublet_influence_matrix,doublet_influence_matrix);
    KSPGetPC(ksp_doublet,&pc);
    PCSetType(pc,PCLU);
    KSPSetType(ksp_doublet,KSPPREONLY);
    KSPSetTolerances(ksp_doublet,1e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(ksp_doublet);
    KSPSetUp(ksp_doublet);
    KSPSolve(ksp_doublet,RHS,solution);
    KSPGetIterationNumber(ksp_doublet,&itn);
//    PetscPrintf(PETSC_COMM_WORLD,"Solution converged in %d iterations.\n",itn+1);

    PetscReal *_SOL;

    VecGetArray(solution,&_SOL);

    /* set doublet strength of panels */
    doublet_strength.clear();
    doublet_strength.resize(total_n_panels);
    for(int p = 0; p < total_n_panels; p++){
        doublet_strength[p] = _SOL[p];
    }

    VecRestoreArray(solution,&_SOL);
}

void Solver :: iterator_conventional_Kutta(const int& iteration, const double &dt, const std::vector<int> ind_start_surface) {
    assert(total_wake_panels_on_edge>0);
    
    Mat Jacobian_global_petsc = PETSC_NULL;
    petsc_mat_create(Jacobian_global_petsc,total_wake_panels_on_edge,total_wake_panels_on_edge);
    
    Vec delta_pressure_TE_petsc = PETSC_NULL;
    petsc_vec_create(delta_pressure_TE_petsc,total_wake_panels_on_edge );
    
    Vec local_sol = PETSC_NULL;
    petsc_vec_create(local_sol,total_wake_panels_on_edge );
    
    Vec delta_mu_k_petsc = PETSC_NULL;
    petsc_vec_create(delta_mu_k_petsc,total_n_panels );
    
    KSP ksp_jacobian = PETSC_NULL;
    KSPCreate(PETSC_COMM_WORLD,&ksp_jacobian);
    
    vector<double> delta_mu_k(total_n_panels,0.);
    
    double maxval = 0.;
    int current_iteration = 0;
        
// Iteration loop to converge towards zero pressure difference at trailing edge
    do {
// compute surface velocity / surface potential / Pressure coefficient
        get_velocity_potential_pressure(iteration,dt,ind_start_surface);
    
// Get the Jacobian (Delta_p_beta - Delta_p^k) / (mu_wake_beta - mu_wake^k)
// The difference in pressure Delta_p^k is saved in the vector delta_pressure_TE for use afterwards to update mu_wake^k+1
        get_pressure_difference_TE(ind_start_surface) ;
        
// Compute the Jacobian only in the first inner iteration and for the 3 first outer time iterations.
// The jacobian is then assumed to remain constant for the next outer time iterations -> conventional C in 
// A fast method to realize the pressure Kutta condition in boundary element method for lifting bodies, Wand Youjiang, Ocean Engineering 2017
//        if (current_iteration==0)
//        if (current_iteration==0 && iteration<4) 
            get_global_Jacobian_pressure_Kutta_conventional(iteration,dt,ind_start_surface);
        
        maxval = 0.;
        for (int TE=0; TE<delta_pressure_TE.size(); TE++){
            if ( maxval < abs(delta_pressure_TE[TE]) )
                maxval = abs(delta_pressure_TE[TE]) ;
        }
        cout << maxval << endl;
        
// Update mu_wake_TE = mu_wake_k+1 = mu_wake_k - inv(J) * Delta_p^k
        VecSet(delta_pressure_TE_petsc,0.0);
        MatZeroEntries(Jacobian_global_petsc);
// Copy the entries of Delta_p^k and Jacobian J^k in the Mat and Vec of Petsc
        for (PetscInt i=0; i<total_wake_panels_on_edge; i++){
            int col[total_wake_panels_on_edge];
            double val[total_wake_panels_on_edge];
            for(int j = 0; j < total_wake_panels_on_edge ; j++){
                col[j] = j;
                val[j] = Jacobian_global[i][j];
            }
            double delta_2_add = delta_pressure_TE[i];
            VecSetValues(delta_pressure_TE_petsc,1,&i,&delta_2_add,INSERT_VALUES);
            MatSetValues(Jacobian_global_petsc,1,&i,total_wake_panels_on_edge,col,val,INSERT_VALUES);
        }
        MatAssemblyBegin(Jacobian_global_petsc,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(Jacobian_global_petsc,MAT_FINAL_ASSEMBLY);
// Solve the system of equations 
        PC pc;
        KSPSetOperators(ksp_jacobian,Jacobian_global_petsc,Jacobian_global_petsc);
        KSPGetPC(ksp_jacobian,&pc);
        PCSetType(pc,PCLU);
        KSPSetType(ksp_jacobian,KSPPREONLY);
        KSPSetTolerances(ksp_jacobian,1e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
        KSPSetFromOptions(ksp_jacobian);
        KSPSetUp(ksp_jacobian);
        
        KSPSolve(ksp_jacobian,delta_pressure_TE_petsc,local_sol);
        
        VecSet(delta_pressure_TE_petsc,0.0); // used as tmp vector to update mu_body
// update mu_wake_TE = mu_wake_k+1 = mu_wake_k + value (where value is the opposite 
//                                        of the solution of the system of eqt)
        PetscReal *_SOL;
        VecGetArray(local_sol,&_SOL);
        for (int j=0; j<total_wake_panels_on_edge; j++){
            double value = -_SOL[j] ;          
            mu_wake_TE[j] += value;
            VecSetValues(delta_pressure_TE_petsc,1,&j,&value,INSERT_VALUES); // delta_pressure_TE_petsc holds mu_wake^(k+1) - mu_wake^(k)
       }
        VecRestoreArray(local_sol,&_SOL);
        
// update doublet_strength = mu on body at iteration k+1
        MatMult(Jacobian_doublets,delta_pressure_TE_petsc,delta_mu_k_petsc);
        VecGetArray(delta_mu_k_petsc,&_SOL);
        for (int j=0; j<total_n_panels; j++){
            doublet_strength[j] += _SOL[j];
       }
        VecRestoreArray(delta_mu_k_petsc,&_SOL);
        
        current_iteration ++;
    }
    while (maxval > Parameters::max_error_pressure_TE && current_iteration <= Parameters::max_number_iterations_pressure_TE) ;
    
    if (current_iteration>Parameters::max_number_iterations_pressure_TE)
        cout << "WARNING The pressure difference at the trailing edge is bigger than the tolerance"<<endl;
    else
        cout << "Number of iterations for pressure condition: "<<current_iteration<<endl;
    
    MatDestroy(&Jacobian_global_petsc);
    VecDestroy(&delta_pressure_TE_petsc);
    VecDestroy(&local_sol);
    VecDestroy(&delta_mu_k_petsc);
    KSPDestroy(&ksp_jacobian);
    
}

void Solver :: get_velocity_potential_pressure(const int iteration, const double dt, const std::vector<int> ind_start_surface ) {
// compute surface velocity / surface potential / Pressure coefficient
// for the current doublet distribution on every surfaces

    if(iteration == 0){
        surface_velocity.clear();
        surface_velocity.resize(total_n_panels);
        
        surface_potential.clear();
        surface_potential.resize(total_n_panels);
        
        pressure_coefficient.clear();
        pressure_coefficient.resize(total_n_panels);
    }
#pragma omp parallel
    for (int i=0;i<surface.size();i++){
        
//    surface[i]->panel_tangents.clear() ;
//    surface[i]->panel_tangents.resize(surface[i]->n_panels());
        
#pragma omp for
        for(int p = 0; p < surface[i]->n_panels(); p++){
            int n_mat = ind_start_surface[i] + p ;
            surface_velocity[n_mat]     = compute_surface_velocity(i,p) ;
//            surface_velocity[n_mat]     = compute_surface_velocity_FD(i,p,ind_start_surface[i]) ;
            
            surface_potential[n_mat]    = compute_surface_potential(n_mat);
            pressure_coefficient[n_mat] = compute_pressure_coefficient(i,p,n_mat,iteration,dt);
        }
    }
}

void Solver :: get_pressure_for_perturbation(const int iteration, const double dt, const std::vector<int> ind_start_surface, const std::vector<double> &doublets, std::vector<double> &pressure  ) {
// compute the pressure coefficient for all surfaces
// for the doublet distribution given in parameter

    pressure.clear();
    pressure.resize(total_n_panels);
    
    std::vector<vector3d> velocity (total_n_panels);
    std::vector<double> potential(total_n_panels,0.);
    
#pragma omp parallel
    for (int i=0; i<surface.size(); i++){
#pragma omp for
        for(int p=0; p<surface[i]->n_panels(); p++){
            int n_mat = ind_start_surface[i] + p ;
            velocity[n_mat]  = compute_surface_velocity(i,p,doublets) ;
//            velocity[n_mat]  = compute_surface_velocity_FD(i,p,n_mat,doublets) ;
            potential[n_mat] = - doublets[n_mat];
            pressure[n_mat]  = compute_pressure_coefficient(i,p,n_mat,iteration,dt,potential[n_mat],velocity[n_mat]);
        }
    }
}

void Solver :: get_pressure_for_perturbation(const int iteration, const double dt, const int isurf, const int ipanel, const std::vector<int> ind_start_surface, const std::vector<double> &doublets, double &pressure  ) {
// compute the pressure coefficient at a given panel n_mat = ind_start_surface[isurf] + ipanel 
// for the doublet distribution given in parameter

    int n_mat         = ind_start_surface[isurf] + ipanel ;
    vector3d velocity = compute_surface_velocity(isurf,ipanel,doublets) ;
//    vector3d velocity = compute_surface_velocity_FD(isurf,ipanel,n_mat,doublets) ;
    double potential  = - doublets[n_mat];
    pressure          = compute_pressure_coefficient(isurf,ipanel,n_mat,iteration,dt,potential,velocity);
    
}

vector3d Solver :: compute_surface_velocity(const int isurf,const int panel) const {

    vector3d local_velocity = surface[isurf]->get_kinematic_velocity(surface[isurf]->get_collocation_point(panel)) - free_stream_velocity;

    vector3d local_velocity_transformed = surface[isurf]->transform_vector_panel(panel,local_velocity);

    /* computes tangential velocity using Least-Square approach
     * refer : CFD : Principles and Applications, by J. Blazek (pg. 162)
     *
     * Least square problem solved with Lapack's divide and conquor routine : dgelsd_
     * Example implementation of dgelsd_ :
     * https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgelsd_ex.c.htm
     * More info on LLSQ : http://www.netlib.org/lapack/lug/node27.html
     */

    const vector<int>& neighbour_panels = surface[isurf]->panel_neighbours[panel];
    int neighbour_size = (int)neighbour_panels.size();
    assert(neighbour_size > 0);

    int dim = 2;    // neglect variation in z
    double rhs[neighbour_size];
    double mat[dim*neighbour_size];

    int ind_start = 0;
    for (int i = 0; i < isurf; i++)
        ind_start += surface[i]->n_panels();
    
    // setup RHS
    for(int i = 0; i < neighbour_size; i++)
        rhs[i] = doublet_strength[ind_start+neighbour_panels[i]] - doublet_strength[ind_start+panel];

    // setup matrix (in column major layout)
    for(int i = 0; i < neighbour_size; i++){

        // transform CP of neighbouring node to panel's coordinates
        vector3d neighbour_node = surface[isurf]->transform_point_panel(panel,surface[isurf]->get_collocation_point(neighbour_panels[i]));

        for(int j = 0; j < dim; j++){
            mat[j*neighbour_size+i] = neighbour_node[j];
        }
    }

    /* solve least square problem */
    /* Local variables to dgelsd_ */
    int m = neighbour_size, n = dim, nrhs = 1, lda = m, ldb = max(m,n), info, lwork, rank;
    double rcond = -1.0;
    double wkopt;
    /* iwork dimension should be at least 3*min(m,n)*nlvl + 11*min(m,n),
     * where nlvl = max( 0, int( log_2( min(m,n)/(smlsiz+1) ) )+1 )
     * and smlsiz = 25 */
    int iwork[11*min(m,n)];
    double s[m];

    /* Query and allocate the optimal workspace */
    lwork = -1;
    dgelsd_( &m, &n, &nrhs, mat, &lda, rhs, &ldb, s, &rcond, &rank, &wkopt, &lwork, iwork, &info );
    lwork = (int)wkopt;
    double work[lwork];

    /* Solve the equations A*X = B */
    dgelsd_( &m, &n, &nrhs, mat, &lda, rhs, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info );

    // check if solution returned success
    assert(info == 0);

    // notice negative sign on rhs terms
    // also notice third component is kept zero
    vector3d total_velocity = vector3d(-rhs[0],-rhs[1],0) - vector3d(local_velocity_transformed[0],local_velocity_transformed[1],0);

    // transform back to global coordinates
    return surface[isurf]->transform_vector_panel_inverse(panel,total_velocity);
}

vector3d Solver :: compute_surface_velocity(const int isurf,const int panel,const std::vector<double> &doublets) const {

    assert(doublets.size()>0);
    
    vector3d local_velocity = surface[isurf]->get_kinematic_velocity(surface[isurf]->get_collocation_point(panel)) - free_stream_velocity;

    vector3d local_velocity_transformed = surface[isurf]->transform_vector_panel(panel,local_velocity);

    const vector<int>& neighbour_panels = surface[isurf]->panel_neighbours[panel];
    int neighbour_size = (int)neighbour_panels.size();
    assert(neighbour_size > 0);

    int dim = 2;    // neglect variation in z
    double rhs[neighbour_size];
    double mat[dim*neighbour_size];

    int ind_start = 0;
    for (int i = 0; i < isurf; i++)
        ind_start += surface[i]->n_panels();
    
    // setup RHS
    for(int i = 0; i < neighbour_size; i++)
        rhs[i] = doublets[ind_start+neighbour_panels[i]] - doublets[ind_start+panel];

    // setup matrix (in column major layout)
    for(int i = 0; i < neighbour_size; i++){

        // transform CP of neighbouring node to panel's coordinates
        vector3d neighbour_node = surface[isurf]->transform_point_panel(panel,surface[isurf]->get_collocation_point(neighbour_panels[i]));

        for(int j = 0; j < dim; j++){
            mat[j*neighbour_size+i] = neighbour_node[j];
        }
    }

    /* solve least square problem */
    /* Local variables to dgelsd_ */
    int m = neighbour_size, n = dim, nrhs = 1, lda = m, ldb = max(m,n), info, lwork, rank;
    double rcond = -1.0;
    double wkopt;
    /* iwork dimension should be at least 3*min(m,n)*nlvl + 11*min(m,n),
     * where nlvl = max( 0, int( log_2( min(m,n)/(smlsiz+1) ) )+1 )
     * and smlsiz = 25 */
    int iwork[11*min(m,n)];
    double s[m];

    /* Query and allocate the optimal workspace */
    lwork = -1;
    dgelsd_( &m, &n, &nrhs, mat, &lda, rhs, &ldb, s, &rcond, &rank, &wkopt, &lwork, iwork, &info );
    lwork = (int)wkopt;
    double work[lwork];

    /* Solve the equations A*X = B */
    dgelsd_( &m, &n, &nrhs, mat, &lda, rhs, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info );

    // check if solution returned success
    assert(info == 0);

    // notice negative sign on rhs terms
    // also notice third component is kept zero
    vector3d total_velocity = vector3d(-rhs[0],-rhs[1],0) - vector3d(local_velocity_transformed[0],local_velocity_transformed[1],0);

    // transform back to global coordinates
    return surface[isurf]->transform_vector_panel_inverse(panel,total_velocity);
}

vector3d Solver :: compute_surface_velocity_FD(const int isurf, const int panel, const int n_mat) const {

// get the local coordinates of the 4 neighbour panels
    vector<vector3d> dist_neighbours ;
    vector<int> local_neighbours;
    surface[isurf]->compute_neighbour_distances(panel, dist_neighbours) ;
    surface[isurf]->get_local_neighbours(panel,local_neighbours);
    
// test if panel at trailing edge / tip / root
    int imax,jmax;
    surface[isurf]->get_IMAX_JMAX(imax, jmax) ;
    imax-- ; jmax--; // number of panels along chord and span respectively
    
    int ipanel,jpanel;
    jpanel = panel/imax; // local indices of the current panel
    ipanel = panel%imax;
//    printf (" ( %d,%d,%d )\n",panel,ipanel,jpanel);
/*    
    printf (" ( %d,%d,%d ) ",panel,ipanel,jpanel);
    for (int i=0;i<4;i++)
        cout << "(" << dist_neighbours[i] << ") "; 
    cout << endl;*/
    
    double deriv1=0.,deriv2=0.;
    
    int ind0 = n_mat + panel; // global index of the panel in the unknown array
    int ind1 = n_mat + local_neighbours[0];
    int ind2 = n_mat + local_neighbours[1];
    double ds1 = fabs( (dist_neighbours[0])[0] ); // distance along chord to the 2 neighbours along chord
    double ds2 = fabs( (dist_neighbours[1])[0] );
    
    if (ipanel==0)  // trailing edge lower
        deriv1 = ( ds1*ds1*doublet_strength[ind2] - ds2*ds2*doublet_strength[ind1] + (ds2*ds2-ds1*ds1)*doublet_strength[ind0] ) / (ds1*ds2*(ds1-ds2));
    else if (ipanel==(imax-1))  // trailing edge upper
        deriv1 = ( ds1*ds1*doublet_strength[ind2] - ds2*ds2*doublet_strength[ind1] + (ds2*ds2-ds1*ds1)*doublet_strength[ind0] ) / (ds1*ds2*(ds2-ds1));
    else  // standard derivative along chordwise direction
        deriv1 = ( ds2*ds2*doublet_strength[ind1] - ds1*ds1*doublet_strength[ind2] - (ds2*ds2-ds1*ds1)*doublet_strength[ind0] ) / (ds1*ds2*(ds1+ds2));
//    cout << ipanel<<" "<<deriv1<<endl;
    
    ind1 = n_mat + local_neighbours[2]; // global index of neighbours panels along span
    ind2 = n_mat + local_neighbours[3];
    ds1    = fabs( (dist_neighbours[2])[1] ); // distance along span to the 2 neighbours along span
    ds2    = fabs( (dist_neighbours[3])[1] );
    
// TEST TEST TEST
// dirTgt is a vector that joins the collocations points of the two neighbours panels along the spanwise direction
// This vector is not perpendicular to the longitudinal vector of the current panel, so the derivative along the spanwise direction must be corrected
// by the projection of the vector dirTgt on the spanwise vector (vec_transv) and on the longitudinal vector (vec_long)
    vector3d dirTgt = surface[isurf]->get_collocation_point(local_neighbours[2]) - surface[isurf]->get_collocation_point(local_neighbours[3]) ;
    dirTgt.normalize() ;
    if (jpanel==0)
        dirTgt.flip();
    vector3d vec_transv = surface[isurf]->get_panel_transverse(panel);
    vector3d vec_long = surface[isurf]->get_panel_longitudinal(panel);
    double cosPsi = vec_transv.dot(dirTgt);
    double sinPsi = vec_long.dot(dirTgt);
    
//    surface[isurf]->panel_tangents[panel] = dirTgt;
    
    if (jpanel==0)  
        deriv2 = ( ds1*ds1*doublet_strength[ind2] - ds2*ds2*doublet_strength[ind1] + (ds2*ds2-ds1*ds1)*doublet_strength[ind0] ) / (ds1*ds2*(ds1-ds2));
    else if (jpanel==(jmax-1))
        deriv2 = ( ds1*ds1*doublet_strength[ind2] - ds2*ds2*doublet_strength[ind1] + (ds2*ds2-ds1*ds1)*doublet_strength[ind0] ) / (ds1*ds2*(ds2-ds1));
    else  // standard derivative along span direction
        deriv2 = ( ds2*ds2*doublet_strength[ind1] - ds1*ds1*doublet_strength[ind2] - (ds2*ds2-ds1*ds1)*doublet_strength[ind0] ) / (ds1*ds2*(ds1+ds2));
//    cout << ipanel << " " << jpanel<<" "<<deriv2<<endl;
    
// Correction of the spanwise derivative because the tangent vector is not perpendicular to the longitudinal vector
    deriv2 = ( deriv2 - deriv1*sinPsi ) / cosPsi ; // TEST TEST TEST 
    
// deriv1 and deriv2 are the components of the surface perturbation velocity   
    vector3d local_velocity = surface[isurf]->get_kinematic_velocity(surface[isurf]->get_collocation_point(panel)) - free_stream_velocity;
    vector3d local_velocity_transformed = surface[isurf]->transform_vector_panel(panel,local_velocity);
    vector3d total_velocity = vector3d(-deriv1,-deriv2,0.) - vector3d(local_velocity_transformed[0],local_velocity_transformed[1],0);

    // transform back to global coordinates
    return surface[isurf]->transform_vector_panel_inverse(panel,total_velocity);
}

vector3d Solver :: compute_surface_velocity_FD(const int isurf, const int panel, const int n_mat,const std::vector<double> &doublets) const {

    assert(doublets.size()>0);
    
// get the local coordinates of the 4 neighbour panels
    vector<vector3d> dist_neighbours ;
    vector<int> local_neighbours;
    surface[isurf]->compute_neighbour_distances(panel, dist_neighbours) ;
    surface[isurf]->get_local_neighbours(panel,local_neighbours);
    
    int imax,jmax;
    surface[isurf]->get_IMAX_JMAX(imax, jmax) ;
    imax-- ; jmax--; // number of panels along chord and span respectively
    
    int ipanel,jpanel;
    jpanel = panel/imax; // local indices of the current panel
    ipanel = panel%imax;
    
    double deriv1=0.,deriv2=0.;
    
    int ind0 = n_mat + panel; // global index of the panel in the unknown array
    int ind1 = n_mat + local_neighbours[0];
    int ind2 = n_mat + local_neighbours[1];
    double ds1 = fabs( (dist_neighbours[0])[0] ); // distance along chord to the 2 neighbours along chord
    double ds2 = fabs( (dist_neighbours[1])[0] );
    
    if (ipanel==0)  // trailing edge lower
        deriv1 = ( ds1*ds1*doublets[ind2] - ds2*ds2*doublets[ind1] + (ds2*ds2-ds1*ds1)*doublets[ind0] ) / (ds1*ds2*(ds1-ds2));
    else if (ipanel==(imax-1))  // trailing edge upper
        deriv1 = ( ds1*ds1*doublets[ind2] - ds2*ds2*doublets[ind1] + (ds2*ds2-ds1*ds1)*doublets[ind0] ) / (ds1*ds2*(ds2-ds1));
    else  // standard derivative along chordwise direction
        deriv1 = ( ds2*ds2*doublets[ind1] - ds1*ds1*doublets[ind2] - (ds2*ds2-ds1*ds1)*doublets[ind0] ) / (ds1*ds2*(ds1+ds2));
//    cout << ipanel<<" "<<deriv1<<endl;
    
    ind1 = n_mat + local_neighbours[2]; // global index of neighbours panels along span
    ind2 = n_mat + local_neighbours[3];
    ds1    = fabs( (dist_neighbours[2])[1] ); // distance along span to the 2 neighbours along span
    ds2    = fabs( (dist_neighbours[3])[1] );
    
// TEST TEST TEST
    vector3d dirTgt = surface[isurf]->get_collocation_point(local_neighbours[2]) - surface[isurf]->get_collocation_point(local_neighbours[3]) ;
    dirTgt.normalize() ;
    if (jpanel==0)
        dirTgt.flip();
    vector3d vec_transv = surface[isurf]->get_panel_transverse(panel);
    vector3d vec_long = surface[isurf]->get_panel_longitudinal(panel);
    double cosPsi = vec_transv.dot(dirTgt);
    double sinPsi = vec_long.dot(dirTgt);
    
    if (jpanel==0)  
        deriv2 = ( ds1*ds1*doublets[ind2] - ds2*ds2*doublets[ind1] + (ds2*ds2-ds1*ds1)*doublets[ind0] ) / (ds1*ds2*(ds1-ds2));
    else if (jpanel==(jmax-1))
        deriv2 = ( ds1*ds1*doublets[ind2] - ds2*ds2*doublets[ind1] + (ds2*ds2-ds1*ds1)*doublets[ind0] ) / (ds1*ds2*(ds2-ds1));
    else  // standard derivative along span direction
        deriv2 = ( ds2*ds2*doublets[ind1] - ds1*ds1*doublets[ind2] - (ds2*ds2-ds1*ds1)*doublets[ind0] ) / (ds1*ds2*(ds1+ds2));
//    cout << ipanel << " " << jpanel<<" "<<deriv2<<endl;
    
    deriv2 = ( deriv2 - deriv1*sinPsi ) / cosPsi ; // TEST TEST TEST 
    
// deriv1 and deriv2 are the components of the surface perturbation velocity   
    vector3d local_velocity = surface[isurf]->get_kinematic_velocity(surface[isurf]->get_collocation_point(panel)) - free_stream_velocity;
    vector3d local_velocity_transformed = surface[isurf]->transform_vector_panel(panel,local_velocity);
    vector3d total_velocity = vector3d(-deriv1,-deriv2,0.) - vector3d(local_velocity_transformed[0],local_velocity_transformed[1],0);

    // transform back to global coordinates
    return surface[isurf]->transform_vector_panel_inverse(panel,total_velocity);
}

double Solver :: compute_pressure_coefficient(const int isurf,const int panel,const int index, const int iteration, const double dt) const {

//compute dphidt
    double dphidt = 0;
    if(iteration > 0 && Parameters::unsteady_problem){
        assert(dt > 0);
        
        if (iteration>1) // first derivative of phi, second order of precision
            dphidt = (1.5*surface_potential[index] - 2*surface_potential_old[index] + 0.5*surface_potential_old_old[index]) / dt ;
        else // first derivative of phi, first order of precision
            dphidt = (surface_potential[index] - surface_potential_old[index]) / dt ;
    }

    vector3d ref_vel = free_stream_velocity - surface[isurf]->get_kinematic_velocity(surface[isurf]->get_collocation_point(panel));

    assert(ref_vel.squared_norm() != 0);

    // compute pressure_coefficient
    double Cp = 1.0 - (surface_velocity[index].squared_norm() + 2.0 * dphidt) / ref_vel.squared_norm();

    return Cp;
}

double Solver :: compute_pressure_coefficient(const int isurf,const int panel,const int index, const int iteration, const double dt, const double potential, const vector3d velocity ) const {

//compute dphidt
    double dphidt = 0;
    if(iteration > 0 && Parameters::unsteady_problem){
        assert(dt > 0);
        
        if (iteration>1) // first derivative of phi, second order of precision
            dphidt = (1.5*potential - 2*surface_potential_old[index] + 0.5*surface_potential_old_old[index]) / dt ;
        else // first derivative of phi, first order of precision
            dphidt = (potential - surface_potential_old[index]) / dt ;
    }

    vector3d ref_vel = free_stream_velocity - surface[isurf]->get_kinematic_velocity(surface[isurf]->get_collocation_point(panel));

    assert(ref_vel.squared_norm() != 0);

    // compute pressure_coefficient
    double Cp = 1.0 - (velocity.squared_norm() + 2.0 * dphidt) / ref_vel.squared_norm();

    return Cp;
}

double Solver :: compute_surface_potential(const int& panel) const {
    double potential = - doublet_strength[panel];

    // add contribution from free stream velocity and body velocity (phi_infinity = U*x+V*y+W*z)
    //vector3d local_velocity = surface->get_kinematic_velocity(surface->get_collocation_point(panel)) - free_stream_velocity;
    //potential += local_velocity.dot(surface->get_collocation_point(panel));

    return potential;
}

double Solver :: compute_surface_potential_infinity(const int& panel, const int &isurf) const {
// Contribution from free stream velocity and body velocity (phi_infinity = U*x+V*y+W*z)
    vector3d local_velocity = surface[isurf]->get_kinematic_velocity(surface[isurf]->get_collocation_point(panel)) - free_stream_velocity;
    double potential = local_velocity.dot(surface[isurf]->get_collocation_point(panel));

    return potential;
}

vector3d Solver :: compute_body_forces(const int isurf) const {
    assert(density > 0);
    
    vector3d Force(0,0,0);
    int ind_start = 0;
    for (int i=0;i<isurf;i++)
        ind_start += surface[i]->n_panels();
    
    double dynamic_pressure = 0. ;
    vector3d ref_vel ;
    // compute force
    for(int p = 0; p < surface[isurf]->n_panels(); p++){
        ref_vel = free_stream_velocity - surface[isurf]->get_kinematic_velocity(surface[isurf]->get_collocation_point(p));
        dynamic_pressure = 0.5 * density * ref_vel.squared_norm();
        Force -= surface[isurf]->get_panel_normal(p) * (dynamic_pressure * pressure_coefficient[ind_start+p] * surface[isurf]->get_panel_area(p));
    }

    return Force;
}

vector3d Solver :: compute_local_beam_forces(const int isurf,const int ibeam) const {
    assert(density > 0);
    vector3d Force(0,0,0);
    
    int ind_start = 0;
    for (int i=0;i<isurf;i++)
        ind_start += surface[i]->n_panels(); // access to the start index for the current blade isurf
    
    int IMAX,JMAX;
    surface[isurf]->get_IMAX_JMAX(IMAX, JMAX) ;
    int ind_start2 = ibeam*(IMAX-1) ; // access to the panel at the lower TE of the current section ibeam
    
    vector3d ref_vel ;
    double dynamic_pressure = 0. ;
    for(int p = ind_start2; p < ind_start2+IMAX-1 ; p++ ) { // all the panels for the current blade isurf at the given section ibeam
        ref_vel = free_stream_velocity - surface[isurf]->get_kinematic_velocity(surface[isurf]->get_collocation_point(p));
        dynamic_pressure = 0.5 * density * ref_vel.squared_norm();
        Force -= surface[isurf]->get_panel_normal(p) * (dynamic_pressure * pressure_coefficient[ind_start+p] * surface[isurf]->get_panel_area(p));
    }
    
    return Force;
}

vector3d Solver :: compute_total_velocity(const int isurf,const vector3d& x, const double cur_time) const {

    assert(source_strength.size()  > 0);
    assert(doublet_strength.size() > 0);

    vector3d velocity(0,0,0);

    // compute velocity due to surface panels
    int ind_start = 0;
    for (int i=0;i<isurf;i++)
        ind_start += surface[i]->n_panels();
        
    for(int sp = 0; sp < surface[isurf]->n_panels(); sp++){
        velocity += surface[isurf]->compute_source_panel_unit_velocity(sp,x)  * source_strength[ind_start+sp];
        velocity += surface[isurf]->compute_doublet_panel_unit_velocity(sp,x,cur_time) * doublet_strength[ind_start+sp];
    }
    // compute velocity due to shedded wake panels
    ind_start = 0;
    for (int i=0;i<isurf;i++)
        ind_start += wake[i]->get_size_doublet_strength();
    
    for(int wp = 0; wp < wake[isurf]->get_size_doublet_strength(); wp++){
        velocity -= wake[isurf]->compute_doublet_panel_unit_velocity(wp,x,cur_time) * wake_doublet_strength[ind_start+wp];
    }
    // add free stream velocity
    velocity += free_stream_velocity;

    return velocity;
}

void Solver :: convect_wake(const double& dt, const double cur_time){
// this function only convects nodes which are already present in the wake (except trailing edge nodes)
// must call shed_wake() function to add new wake panel row

    if(Parameters::static_wake == false){
        int nodes_to_convect;
        for (int i=0; i<wake.size();i++) {
        
            assert(wake[i]->n_panels() > 0);

// nodes_to_convect does not include nodes on trailing edge.
// The nodes on the TE are the last being added, so they are indexed between wake->n_nodes() - surface->n_trailing_edge_nodes()   AND   wake->n_nodes()-1
            nodes_to_convect = wake[i]->n_nodes() - surface[i]->n_trailing_edge_nodes();
            assert(nodes_to_convect > 0);

            // compute velocity at the wake nodes which needs to be convected
            vector<vector3d> wake_velocity(nodes_to_convect);
#pragma omp parallel
{
#pragma omp for
            for(int wn = 0; wn < nodes_to_convect; wn++)
                wake_velocity[wn] = compute_total_velocity(i,wake[i]->nodes[wn],cur_time);
#pragma omp barrier
            // move the nodes with wake velocity
#pragma omp for
            for(int wn = 0; wn < nodes_to_convect; wn++)
                wake[i]->nodes[wn] += wake_velocity[wn] * dt ;

}
        }
    }

}

void Solver :: write_output(const int& iteration) const {

    // folder name
    const std::string foldername = "Output";

    // create directory if does not exist
    struct stat dir_err;
    if(stat(foldername.c_str(),&dir_err) == -1)
#ifdef COMP_WINDOWS // compile on Windows
        mkdir(foldername.c_str());
#else // compile on Linux
        mkdir(foldername.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);   
#endif

    // set surface and wake names
    std::string file_name ;
    std::string surf_name ;
    // write data
//    cout << "writing output...";
    
    ///////////////////////////////////////////////////////////////////////
    std::vector<double> data ;
    std::vector<double>::const_iterator iit1 ;
    std::vector<double>::const_iterator iit2 ;
    std::vector<vector3d> data2 ;
    std::vector<vector3d>::const_iterator iit3 ;
    std::vector<vector3d>::const_iterator iit4 ;
    
    int ibegin=0,iend=0;
// write the data for the surfaces
    for (int isurf=0;isurf<surface.size();isurf++) { 
        iend += surface[isurf]->n_panels(); // update this iterator to point to the last element to be copied
        if ((surface[isurf]->get_name_surface()).empty())
            surf_name = "surface_" + to_string(isurf+1);
        else
            surf_name = surface[isurf]->get_name_surface() ;
        
        file_name = foldername+"/"+surf_name+"_at_iteration_"+to_string(iteration);

        // normals
        data2.clear();
        data2.resize(surface[isurf]->n_panels());
        for (int j=0;j<surface[isurf]->n_panels();j++) {
            data2[j]=surface[isurf]->get_panel_normal(j);
        }
        log->write_surface_data(file_name,surface[isurf],data2,"normals",true);
        data2.clear();
        
        // longitudinals
        data2.clear();
        data2.resize(surface[isurf]->n_panels());
        for (int j=0;j<surface[isurf]->n_panels();j++) {
            data2[j]=surface[isurf]->get_panel_longitudinal(j);
        }
        log->write_surface_data(file_name,surface[isurf],data2,"longitudinals",false);
        data2.clear();

        // transverse
        data2.clear();
        data2.resize(surface[isurf]->n_panels());
        for (int j=0;j<surface[isurf]->n_panels();j++) {
            data2[j]=surface[isurf]->get_panel_transverse(j);
        }
        log->write_surface_data(file_name,surface[isurf],data2,"transverse",false);
        data2.clear();

        // tangents
/*        data2.clear();
        data2.resize(surface[isurf]->n_panels());
        for (int j=0;j<surface[isurf]->n_panels();j++) {
            data2[j]=surface[isurf]->panel_tangents[j];
        }
        log->write_surface_data(file_name,surface[isurf],data2,"tangent",false);
        data2.clear();*/
        
    // Source strength 
        iit1 = source_strength.begin() + ibegin ;
        iit2 = source_strength.begin() + iend ;
        data.assign(iit1,iit2);
        log->write_surface_data(file_name,surface[isurf],data,"Sigma",false);
        data.clear();
        
    // Doublet strength
        iit1 = doublet_strength.begin() + ibegin ;
        iit2 = doublet_strength.begin() + iend ;
        data.assign(iit1,iit2);
        log->write_surface_data(file_name,surface[isurf],data,"Mu",false);
        data.clear();
           
    // Velocity vector
        iit3 = surface_velocity.begin() + ibegin ;
        iit4 = surface_velocity.begin() + iend ;
        data2.assign(iit3,iit4);
        log->write_surface_data(file_name,surface[isurf],data2,"V",false);
        data2.clear();
        
    // Pressure coefficient
        iit1 = pressure_coefficient.begin() + ibegin ;
        iit2 = pressure_coefficient.begin() + iend ;
        data.assign(iit1,iit2);
        log->write_surface_data(file_name,surface[isurf],data,"CP",false);
        data.clear();
        
        ibegin = iend; // update this iterator to point to the first element to be copied
    }
    
    if(wake_doublet_strength.size() > 0) {
        ibegin=0,iend=0;
        for (int isurf=0;isurf<wake.size();isurf++) {
            iend += wake[isurf]->get_size_doublet_strength();
            
            if ((surface[isurf]->get_name_surface()).empty())
                surf_name = "wake_" + to_string(isurf+1);
            else
                surf_name = wake[isurf]->get_name_surface() ;
            file_name = foldername+"/"+surf_name+"_at_iteration_"+to_string(iteration);
            
            iit1 = wake_doublet_strength.begin() + ibegin ;
            iit2 = wake_doublet_strength.begin() + iend ;
            data.assign(iit1,iit2);
            data.insert(data.end(),surface[isurf]->n_trailing_edge_panels(),0.);
            log->write_surface_data(file_name,wake[isurf],data,"Mu",true);
            
            // normals
            data2.clear();
            data2.resize(wake[isurf]->n_panels());
            for (int j=0;j<wake[isurf]->n_panels();j++) {
                data2[j]=wake[isurf]->get_panel_normal(j);
            }
            log->write_surface_data(file_name,wake[isurf],data2,"normals",false);
            data2.clear();
            
            ibegin = iend;
        }
    }

//    cout << "Done." << endl;

}

void Solver :: finalize_iteration(const int iteration){

// If the problem is steady
    if(!Parameters::unsteady_problem){
        // clear doublet strength of the wake
        wake_doublet_strength.clear();
    }

// If the problem is unsteady
    else if (Parameters::unsteady_problem) {
        // save surface potential in previous variable
        if (iteration>0)
            surface_potential_old_old = surface_potential_old;
        
        surface_potential_old = surface_potential;
        
    }
}

vector3d Solver:: get_body_forces() const{
    return body_forces;
}

vector3d Solver :: get_body_force_coefficients() const{
    return body_force_coefficients;
}

double Solver :: get_pressure_coefficient(const int surf,const int panel) const{
    assert(surf>=0);
    assert(surf<surface.size());
    
    int ind_start = 0;
    for (int i=0;i<surf;i++)
        ind_start += surface[i]->n_panels(); 
    return pressure_coefficient[ind_start+panel];
}

std::vector<double> Solver :: get_pressure_coefficient_surface(const int surf) const{
    assert(surf>=0);
    assert(surf<surface.size());
    
    std::vector<double>::const_iterator it1 = pressure_coefficient.begin();
    int ind_start = 0;
    for (int i=0;i<surf;i++)
        ind_start += surface[i]->n_panels(); 
    it1 += ind_start;
    
    std::vector<double>::const_iterator it2;
    if ((surf+1)==surface.size())
        it2 = pressure_coefficient.end();
    else
        it2 = it1 + surface[surf]->n_panels();
    
    std::vector<double> CP_2_send (it1,it2) ;
    return CP_2_send;
}

void Solver :: write_global_forces(const int& iteration, const double &deltaT) {
    
    const std::string file_name = "Output/Global_forces.dat";
    
    if (!this->ofile_global_forces.is_open()) {
        this->ofile_global_forces.open(file_name);
        this->ofile_global_forces << "# Iteration   Time    Fx   Fy   Fz\n";
/*        
        struct stat buffer;
        if (stat (file_name.c_str(), &buffer) == 0){
            cout << "file exists"<<endl;
            this->ofile_global_forces.open(file_name, std::ios::app);
        }
        else {
            cout << "file does not exist"<<endl;
            this->ofile_global_forces.open(file_name);
        }
*/        
    }
    
    this->ofile_global_forces << iteration << " " << iteration*deltaT << " " << body_forces << endl;
}

void Solver :: get_source_array(const std::vector<int> ind_start_surface) {
    source_strength.clear();
    source_strength.resize(total_n_panels);

    for (int i=0;i<surface.size();i++){
#pragma omp parallel for
        for(int p = 0; p < surface[i]->n_panels(); p++){
            int n_mat = ind_start_surface[i] + p;
            source_strength[n_mat] = compute_source_strength(i,p);
        }
    }
#pragma omp barrier
}

void Solver :: get_source_doublet_matrices(const std::vector<int> ind_start_surface) {
    source_influence.clear();
    source_influence.resize(total_n_panels,vector<double>(total_n_panels));
    doublet_influence.clear();
    doublet_influence.resize(total_n_panels,vector<double>(total_n_panels));
    
#pragma omp parallel
    for (int i=0;i<surface.size();i++){
        for (int j=0;j<surface.size();j++){
#pragma omp for
            for(int n = 0; n < surface[i]->n_panels() ; n++){
                int n_mat = ind_start_surface[i] + n; // WARNING n_mat must stay private in this openMP loop
                for(int p = 0; p < surface[j]->n_panels(); p++){ // for each panel of surface J = influence on current collocation point
                    pair<double,double> influence = surface[j]->compute_source_doublet_panel_influence(p,surface[i]->get_collocation_point(n));
                    
                    int p_mat = ind_start_surface[j] + p; // WARNING p_mat must stay private in this openMP loop
                    if(p_mat == n_mat)
                        influence.second = -0.5;
                    source_influence[n_mat][p_mat]  = influence.first;
                    doublet_influence[n_mat][p_mat] = influence.second;
                }
            }
        }
    }
#pragma omp barrier
}

void Solver :: get_wake_influence_matrix(const std::vector<int> ind_start_surface,const std::vector<int> ind_start_wake) {
    if(wake_doublet_strength.size() > 0){        
// Reset and resize the vector wake_doublet_influence to the actual number of doublets in all the wakes
        wake_doublet_influence.clear();
        int total_n_doublet_wake = 0;
        for (int i=0;i<wake.size();i++)
            total_n_doublet_wake += wake[i]->get_size_doublet_strength();
        wake_doublet_influence.resize(total_n_panels,vector<double>(total_n_doublet_wake));
        
#pragma omp parallel
        for (int i=0;i<surface.size();i++){
#pragma omp for
            for(int n = 0; n < surface[i]->n_panels() ; n++){ // for each panel of surface I = current collocation point
                int n_mat = ind_start_surface[i] + n ;
                for (int j=0;j<wake.size();j++){
                    for(int p = 0; p < wake[j]->get_size_doublet_strength(); p++){ // for each doublet of wake J = influence on current collocation point                        
                        int p_mat = ind_start_wake[j] + p ;
                        wake_doublet_influence[n_mat][p_mat]  = -wake[j]->compute_doublet_panel_influence(p,surface[i]->get_collocation_point(n)) ;
                    }
                }
            }
        }
    }
#pragma omp barrier
}

void Solver :: get_new_wake_doublet(const int iteration,const std::vector<int> ind_start_surface) {
    
    std::vector<double> current_wake; // vector for the current wake, added to the end of wake_doublet_strength
    std::vector<double> tmp_wake = wake_doublet_strength ; // contains a copy of the data
    int cur_pos = 0;
    int wake_panel_start, wake_panel_end ;
    
    wake_doublet_strength.clear();
    
    vector<int> ind_start_wake_TE(surface.size(),0);
    for (int i=1;i<surface.size();i++)
        ind_start_wake_TE[i] = ind_start_wake_TE[i-1] + surface[i]->n_trailing_edge_panels();
    
    for (int j=0;j<wake.size();j++){
        current_wake.clear();
// get the old values for the current wake
        if (iteration>0) 
            current_wake.assign(tmp_wake.begin()+cur_pos,tmp_wake.begin()+cur_pos+wake[j]->get_size_doublet_strength());
        cur_pos += wake[j]->get_size_doublet_strength() ;
        
        if(Parameters::unsteady_problem){
            wake_panel_start = wake[j]->n_panels() - surface[j]->n_trailing_edge_panels();
            wake_panel_end = wake[j]->n_panels();
        }
        else{
            wake_panel_start = 0;
            wake_panel_end = wake[j]->n_panels();
        }
        
//////////////////////////////////////////// OLD METHOD WHEN NEW WAKE COMPUTED BY MORINO
/*
        int n_mat = ind_start_surface[j];
// Compute the new values for the current wake and push them back in the current vector
        for(int wp = wake_panel_start; wp < wake_panel_end; wp++){
            assert(doublet_strength.size() > 0);
            
            int TE_panel_counter = wp % surface[j]->n_trailing_edge_panels();
            int upper_panel = surface[j]->upper_TE_panels[TE_panel_counter];
            int lower_panel = surface[j]->lower_TE_panels[TE_panel_counter];
            
// The following standard Kutta approach is valid for steady problems and for airfoils with zero thickness at the trailing edge
            double wake_strength = doublet_strength[n_mat+upper_panel] - doublet_strength[n_mat+lower_panel];

// To correct the problem with thick trailing edge, one can add the contribution of the potential at infinity phi_inf = U_glob . coord_glob
//            double phi_inf_lower = compute_surface_potential_infinity(lower_panel, j) ;
//            double phi_inf_upper = compute_surface_potential_infinity(upper_panel, j) ;
//            wake_strength += phi_inf_upper - phi_inf_lower ;
            vector3d rte = surface[j]->get_collocation_point(upper_panel) - surface[j]->get_collocation_point(lower_panel) ;
            wake_strength += free_stream_velocity.dot(rte);

// One can compute the new wake value as the average between the one from Kutta and the old advected one
            if (iteration>0) {
                int previous_wake = current_wake.size() - surface[j]->n_trailing_edge_panels() ;
                wake_strength = 0.5*(wake_strength + current_wake[previous_wake]);
            }
            
            current_wake.push_back(wake_strength);
            wake[j]->push_back_doublet_strength(wake_strength);
        }
*/

//////////////////////////////////////////// NEW METHOD WHEN NEW WAKE COMPUTED THROUGH ITERATIVE METHOD
        for(int wp = 0; wp < surface[j]->n_trailing_edge_panels(); wp++){
            double wake_strength = mu_wake_TE[wp + ind_start_wake_TE[j] ] ;
            
            if ( (Parameters :: switch_off_root_vortex) && (wp==(surface[j]->n_trailing_edge_panels()-1) ) )
                wake_strength = 0. ; 
            if ( (Parameters :: switch_off_tip_vortex) && (wp==0) ) 
                wake_strength = 0. ; 
            
            current_wake.push_back(wake_strength);
            wake[j]->push_back_doublet_strength(wake_strength);
        }
        
////////////////////////////////////////////////////////////////////////////////        
// remove the first wake cells which have oscillatory values
/*        if (iteration==10 |iteration==11 |iteration==12 | (iteration>15 & iteration%10==0)) { 
            for(int wp = wake_panel_start; wp < wake_panel_end; wp++){ // remove the panel values
                wake[j]->pop_front_double_strength() ;
                current_wake.erase(current_wake.begin());
            }
            for(int wp = wake_panel_start; wp < wake_panel_end+1; wp++){ // remove the corresponding wake nodes
                wake[j]->nodes.erase( wake[j]->nodes.begin() ) ;
            }
            wake[j]->build_topology();
            wake[j]->compute_panel_components();
        }*/
////////////////////////////////////////////////////////////////////////////////
        
// push back the new current wake in the global wake vector
        wake_doublet_strength.insert(wake_doublet_strength.end(), current_wake.begin(), current_wake.end());
        
    }
    
}

void Solver :: apply_Kutta_Morino(const std::vector<int> ind_start_surface) {
// compute wake panel range to consider in Kutta-condition, the last row of panels that has been added during the last time step
// FOR NOW IT IS ASSUMED THAT THE KUTTA CONDITION LINKS THE WAKE AND ITS GENERATING SURFACE BUT NOT THE OTHER SURFACES
    int wake_panel_start, wake_panel_end,TE_panel_counter;
    double influence2;
#pragma omp parallel
    for (int i=0;i<surface.size();i++){ // for each wake 
        if(Parameters::unsteady_problem){
            wake_panel_start = wake[i]->n_panels() - surface[i]->n_trailing_edge_panels();
            wake_panel_end = wake[i]->n_panels();
        }
        else{
            wake_panel_start = 0;
            wake_panel_end = wake[i]->n_panels();
        }
        
        int p_mat = ind_start_surface[i] ;        
#pragma omp for
        for(int sp = 0; sp < surface[i]->n_panels(); sp++){
            int n_mat = ind_start_surface[i] + sp ;
            for(int wp = wake_panel_start; wp < wake_panel_end; wp++){

                TE_panel_counter = wp % surface[i]->n_trailing_edge_panels();
                int upper_panel = surface[i]->upper_TE_panels[TE_panel_counter];
                int lower_panel = surface[i]->lower_TE_panels[TE_panel_counter];

// remember to use negative sign when computing doublet coeff of the wake (as normal is in opposite direction)
                influence2 = - wake[i]->compute_doublet_panel_influence(wp,surface[i]->get_collocation_point(sp));

                doublet_influence[n_mat][p_mat+upper_panel] += influence2;
                doublet_influence[n_mat][p_mat+lower_panel] -= influence2;
            }
        }
    }
#pragma omp barrier
}

void Solver :: get_Jacobian_wake_on_body(const std::vector<int> ind_start_surface,const std::vector<int> ind_start_wake) {
// Petsc must be initialized before this step  
    PetscBool mybool;
    PetscInitialized(&mybool);
    assert (mybool) ;
    
    total_wake_panels_on_edge = 0;
    for (int i=0; i<surface.size();i++)
        total_wake_panels_on_edge += surface[i]->n_trailing_edge_panels();
    
    assert(total_wake_panels_on_edge > 0);
    assert(total_n_panels > 0);
    
// matrix J_mu = - inv(A) * W
    petsc_mat_create(Jacobian_doublets,total_n_panels,total_wake_panels_on_edge);
    
// matrix that holds the influence coefficients of the body, without the Morino Kutta condition
    Mat local_mat_A = PETSC_NULL ;
    petsc_mat_create(local_mat_A,total_n_panels,total_n_panels);
    
// matrix that holds the influence coefficients of the first wake row on the body
    Mat local_mat_W = PETSC_NULL ;
    petsc_mat_create(local_mat_W,total_n_panels,total_wake_panels_on_edge);
    
// setup column number for matrix insertion
    int col[total_n_panels];
    for(int i = 0; i < total_n_panels ; i++)
        col[i] = i;
    
// copy doublet_influence to local_mat_A
    for(int i = 0;i < total_n_panels; i++){
        double val[total_n_panels];
        std::copy( doublet_influence[i].begin(), doublet_influence[i].end(), val );
        MatSetValues(local_mat_A,1,&i,total_n_panels,col,val,INSERT_VALUES);
    }
    MatAssemblyBegin(local_mat_A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(local_mat_A,MAT_FINAL_ASSEMBLY);
    
//    WriteMat(local_mat_A,"local_mat_A");
    
// copy doublet_influence to local_mat_A
    int ind_start_first_wake=0;
    for(int i = 0;i < surface.size(); i++){
        int col[surface[i]->n_trailing_edge_panels()];
        for(int j = 0; j < surface[i]->n_trailing_edge_panels() ; j++)
            col[j] = ind_start_first_wake + j;
        
        double val[surface[i]->n_trailing_edge_panels()];
// The first row of wake, near the trailing edge, is the last one in wake[i]->panels
        for(int n = 0; n < surface[i]->n_panels() ; n++){
            int n_mat = ind_start_surface[i] + n;            
            for(int p = 0; p < surface[i]->n_trailing_edge_panels(); p++){
                int p_mat = wake[i]->get_size_doublet_strength() + p ;
// The minus sign in the RHS is taken into account here
                val[p] = wake[i]->compute_doublet_panel_influence(p_mat,surface[i]->get_collocation_point(n)) ;
            }
            MatSetValues(local_mat_W,1,&n_mat,surface[i]->n_trailing_edge_panels(),col,val,INSERT_VALUES);
        }
        ind_start_first_wake += surface[i]->n_trailing_edge_panels() ;
    }
    MatAssemblyBegin(local_mat_W,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(local_mat_W,MAT_FINAL_ASSEMBLY);
//    WriteMat(local_mat_W,"local_mat_W");
    
// At this point, matrices A and W are build. The Jacobian J_mu can be 
// computed through the solving of A * J_mu = - W
// WARNING the minus sign in the RHS is already taken into account in W
    KSP ksp_jacobian = PETSC_NULL;
    // create KSP solver
    KSPCreate(PETSC_COMM_WORLD,&ksp_jacobian);
    
    int itn;
    PC pc;
    KSPSetOperators(ksp_jacobian,local_mat_A,local_mat_A);
    KSPGetPC(ksp_jacobian,&pc);
    PCSetType(pc,PCLU);
    KSPSetType(ksp_jacobian,KSPPREONLY);
    KSPSetTolerances(ksp_jacobian,1e-16,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(ksp_jacobian);
    KSPSetUp(ksp_jacobian);
    
    Vec local_RHS,local_sol;
    petsc_vec_create(local_RHS,total_n_panels );
    petsc_vec_create(local_sol,total_n_panels );
// Solve the system of equations with multiple RHS
    for (PetscInt i=0; i<total_wake_panels_on_edge; i++){
// Get the current RHS
        MatGetColumnVector(local_mat_W,local_RHS,i);
        
// Solve the system of equation for the current RHS
        KSPSolve(ksp_jacobian,local_RHS,local_sol);
        
// Transfer the solution from Vec to standard c++ array
        double val[total_n_panels];
        PetscReal *_SOL;
        VecGetArray(local_sol,&_SOL);
        for (int j=0; j<total_n_panels; j++)
            val[j] = _SOL[j];
        VecRestoreArray(local_sol,&_SOL);
// Insert the solution in the corresponding column of the Jacobian 
        MatSetValues(Jacobian_doublets,total_n_panels,col,1,&i,val,INSERT_VALUES);
        
//        WriteVec(local_sol,"local_sol_");
    }
    MatAssemblyBegin(Jacobian_doublets,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Jacobian_doublets,MAT_FINAL_ASSEMBLY);
    
//    WriteMat(Jacobian_doublets,"local_Jacobian");
    
// Free the local Petsc matrices
    MatDestroy(&local_mat_A);
    MatDestroy(&local_mat_W);
    VecDestroy(&local_RHS);
    VecDestroy(&local_sol);
    KSPDestroy(&ksp_jacobian);
}

void Solver :: get_doublet_wake_TE(const std::vector<int> ind_start_surface){
// Get doublets of wake at trailing edge mu_wake = mu_body_upper - mu_body_lower (Morino)
    assert(total_wake_panels_on_edge>0);
    
//    mu_wake_TE.clear();
//    mu_wake_TE.resize(total_wake_panels_on_edge);
    
    std::vector<double> current_wake; // vector for the current wake, added to the end of wake_doublet_strength
    int wake_panel_start, wake_panel_end ;
    
    mu_wake_TE.clear();
    
    for (int j=0;j<wake.size();j++){
        current_wake.clear();
        
        wake_panel_start = wake[j]->n_panels() - surface[j]->n_trailing_edge_panels();
        wake_panel_end = wake[j]->n_panels();
        
        int n_mat = ind_start_surface[j];
// Compute the new values for the current wake and push them back in the current vector
        for(int wp = wake_panel_start; wp < wake_panel_end; wp++){
            assert(doublet_strength.size() > 0);
            
            int TE_panel_counter = wp % surface[j]->n_trailing_edge_panels();
            int upper_panel = surface[j]->upper_TE_panels[TE_panel_counter];
            int lower_panel = surface[j]->lower_TE_panels[TE_panel_counter];
            
            double wake_strength = doublet_strength[n_mat+upper_panel] - doublet_strength[n_mat+lower_panel];
            current_wake.push_back(wake_strength);
        }
        
// push back the new current wake in the global wake vector
        mu_wake_TE.insert(mu_wake_TE.end(), current_wake.begin(), current_wake.end());
    }
    
}

void Solver :: get_pressure_difference_TE(const std::vector<int> ind_start_surface) {
    
    vector<int> ind_start_wake_TE(surface.size(),0);
    for (int i=1;i<surface.size();i++)
        ind_start_wake_TE[i] = ind_start_wake_TE[i-1] + surface[i]->n_trailing_edge_panels();
    
    delta_pressure_TE.clear();
    delta_pressure_TE.resize(total_wake_panels_on_edge);
    
    for (int cur_wake=0; cur_wake<wake.size(); cur_wake++) { // loop over each wake
        for (int i=0; i<surface[cur_wake]->n_trailing_edge_panels(); i++) {// row-loop over each body panel at edge
            int n_mat = ind_start_wake_TE[cur_wake] + i ;
            
            int upper_panel = surface[cur_wake]->upper_TE_panels[i];
            int lower_panel = surface[cur_wake]->lower_TE_panels[i];
                
            vector3d ref_vel_upper = free_stream_velocity - surface[cur_wake]->get_kinematic_velocity(surface[cur_wake]->get_collocation_point(upper_panel));
            vector3d ref_vel_lower = free_stream_velocity - surface[cur_wake]->get_kinematic_velocity(surface[cur_wake]->get_collocation_point(lower_panel));
            assert(ref_vel_upper.squared_norm() != 0);
            assert(ref_vel_lower.squared_norm() != 0);
            
            double p_upper = pressure_coefficient[ind_start_surface[cur_wake]+upper_panel]*2*density*ref_vel_upper.squared_norm();
            double p_lower = pressure_coefficient[ind_start_surface[cur_wake]+lower_panel]*2*density*ref_vel_lower.squared_norm();
            
            delta_pressure_TE[n_mat] = p_upper - p_lower ;
        }
    }
    
}

void Solver :: get_global_Jacobian_pressure_Kutta_conventional(const int& iteration, const double &dt, const std::vector<int> ind_start_surface){
    
    assert(total_wake_panels_on_edge > 0);
    assert(total_n_panels > 0);
    assert(mu_wake_TE.size() > 0);
    
    vector<int> ind_start_wake_TE(surface.size(),0);
    for (int i=1;i<surface.size();i++)
        ind_start_wake_TE[i] = ind_start_wake_TE[i-1] + surface[i]->n_trailing_edge_panels();
    
    vector<double> mu_wake_beta  = mu_wake_TE; // contains the doublets of all wake panels at trailing edge
    vector<double> mu_body_beta  (total_n_panels,0.);
    vector<double> pressure_body_beta (total_n_panels,0.);
    
    Vec Delta_mu_wake = PETSC_NULL;
    Vec Delta_RHS     = PETSC_NULL;
    petsc_vec_create(Delta_mu_wake,total_wake_panels_on_edge);
    petsc_vec_create(Delta_RHS,total_n_panels);
    
    Jacobian_global.clear();
    Jacobian_global.resize(total_wake_panels_on_edge,vector<double>(total_wake_panels_on_edge));
    
    double beta = 0.01 ;
    double factor = 1+beta;
        
    for (int cur_wake=0; cur_wake<wake.size(); cur_wake++) { // loop over each wake
        for (int i=0; i<surface[cur_wake]->n_trailing_edge_panels(); i++) {// row-loop over each body panel at edge
        
            int n_mat = ind_start_wake_TE[cur_wake] + i ;
            int upper_panel = surface[cur_wake]->upper_TE_panels[i];
            int lower_panel = surface[cur_wake]->lower_TE_panels[i];
            vector3d ref_vel_upper = free_stream_velocity - surface[cur_wake]->get_kinematic_velocity(surface[cur_wake]->get_collocation_point(upper_panel));
            vector3d ref_vel_lower = free_stream_velocity - surface[cur_wake]->get_kinematic_velocity(surface[cur_wake]->get_collocation_point(lower_panel));
            assert(ref_vel_upper.squared_norm() != 0);
            assert(ref_vel_lower.squared_norm() != 0);
           
            for (int j=0; j<surface[cur_wake]->n_trailing_edge_panels(); j++){ //column-loop over each wake panel at edge
                int p_mat = ind_start_wake_TE[cur_wake] + j ;
                
// Set a perturbation on mu_wake_beta 
                mu_wake_beta[p_mat] *= factor;
                
// Compute the perturbed distribution of the doublets on the body mu_body_beta = mu_body^k + J_mu * (mu_wake_beta - mu_wake_k)
                double delta_mu_value = -beta*mu_wake_TE[p_mat];
                VecSetValues(Delta_mu_wake,1,&p_mat,&delta_mu_value,INSERT_VALUES);
                MatMult(Jacobian_doublets,Delta_mu_wake,Delta_RHS);
                PetscReal *_SOL;
                VecGetArray(Delta_RHS,&_SOL);
                for (int jj=0; jj<total_n_panels; jj++) 
                    mu_body_beta[jj] = _SOL[jj] + doublet_strength[jj];
                VecRestoreArray(Delta_RHS,&_SOL);
                
// Compute delta p_beta corresponding to the perturbed mu_body_beta
// this takes a lot of time -> better if the Jacobian is computed once during the first iteration
//                get_pressure_for_perturbation(iteration, dt, ind_start_surface, mu_body_beta, pressure_body_beta);
//                double p_upper_beta = pressure_body_beta[ind_start_surface[cur_wake]+upper_panel]*2*density*ref_vel_upper.squared_norm();
//                double p_lower_beta = pressure_body_beta[ind_start_surface[cur_wake]+lower_panel]*2*density*ref_vel_lower.squared_norm();

// Now compute only the 2 required pressure instead of the pressure on all the panels (much more expensive)
                double p_upper_beta,p_lower_beta;
                get_pressure_for_perturbation(iteration, dt, cur_wake, upper_panel, ind_start_surface, mu_body_beta, p_upper_beta ) ;
                get_pressure_for_perturbation(iteration, dt, cur_wake, lower_panel, ind_start_surface, mu_body_beta, p_lower_beta ) ;
                p_upper_beta *= 2*density*ref_vel_upper.squared_norm();
                p_lower_beta *= 2*density*ref_vel_upper.squared_norm();
                
// Compute Jacobian_global(i,j)^k
                Jacobian_global[n_mat][p_mat] = ( p_upper_beta - p_lower_beta - delta_pressure_TE[n_mat] ) / delta_mu_value;
                
// Cancel the perturbation on mu_wake_beta for the next j
                mu_wake_beta[p_mat] = mu_wake_TE[p_mat];
                VecSet(Delta_mu_wake,0.0);
            }
        }        
    }
    VecDestroy(&Delta_mu_wake);
    VecDestroy(&Delta_RHS);
}

void Solver :: get_global_Jacobian_pressure_Kutta_analytical(const int iteration,const double &deltaT){
// A fast method to realize the pressure Kutta condition in boundary element method for lifting bodies
// by Youjiang Wang, Moustafa Abdel-Maksoud and Baowei Song, Ocean Engineering 130, pp398-406, 2017
    
    assert(total_wake_panels_on_edge > 0);
    assert(total_n_panels > 0);
    
    Jacobian_global.clear();
    Jacobian_global.resize(total_wake_panels_on_edge,vector<double>(total_wake_panels_on_edge));
    
    
    double toadd = 0;
    double dpdvsq = -0.5 * density; // eqt (18) : d(p)/d(v^2) = -rho / 2
    for (int cur_wake=0; cur_wake<wake.size(); cur_wake++) { // loop over each wake
// d(p)/d( d(phi)/dt ) * d( d(phi)/dt )/ d(mu_body)
        if (iteration>1) // second order backward discretization for d(phi)/dt
            toadd = 1.5 * density / deltaT; // combine eqts (19) and (20)
        else // first order backward discretization for d(phi)/dt
            toadd = density / deltaT; // combine eqts (19) and (20)
            
        for (int i=0; i<surface[cur_wake]->n_trailing_edge_panels(); i++) {// row-loop over each body panel at edge
// Get index of upper/lower panels at the current wake position
            int upper_panel = surface[cur_wake]->upper_TE_panels[i];
            int lower_panel = surface[cur_wake]->lower_TE_panels[i];
            
            
            for (int j=0; j<surface[cur_wake]->n_trailing_edge_panels(); j++){ //column-loop over each wake panel at edge
                // loop over the neighbours of the current body panel i
                // first compute contribution of upper, then lower panels, then Jacobian value at i,j
            }
        }
    }
}

// for test only
void Solver :: move_beam_test(const double dt, int iteration){
/*
// oscillatory rigid flap motion
    int beam_size;
    vector3d node1,node2,vector_perp,vector_beam;
    double amplitude = 2; // amplitude of oscillation in degrees
    amplitude *= M_PI / 180.0 ;
    double omega = 23.565325441667691928903819677371; // omega = half the rotation frequency
    double sinus = sin(iteration*dt*omega);
    double motion = amplitude*sinus;
    
    vector3d vectorY = vector3d(0.,motion,0.);
    
    for (int i=0;i<surface.size();i++){
        beam_size = surface[i]->beam_nodes.size();
        node1 = surface[i]->beam_nodes[beam_size-1];
        node2 = surface[i]->beam_nodes[beam_size-2];
        vector_beam = node2-node1;
        vector_beam.normalize();
        
        vector_perp = vector_beam.cross(vectorY);
        
        surface[i]->rotate_surface(vector_perp,true);
    }
*/


// oscillatory elastic flap motion
    int beam_size;
    vector3d node1,node2,vector_perp,vector_beam,vector_perp_frac;
    double amplitude = 2; // amplitude of oscillation in degrees
    amplitude *= M_PI / 180.0 ;
    double omega = 2.3876666666666666666666666666667  * 2. * M_PI ; // omega = double of the rotation frequency // NREL Phase VI
//    double omega = 0.32   * 2. * M_PI ; // omega = double of the rotation frequency // NREL 5MW
    double sinus = sin(iteration*dt*omega);
    double motion = amplitude*sinus;
    double rmin,rmax,rcurrent,rfrac;
    vector3d vectorY = vector3d(0.,motion,0.);
    vector3d center = vector3d(0.,0.,0.);
    
    for (int i=0;i<surface.size();i++){
        beam_size = surface[i]->beam_nodes.size();
        node1 = surface[i]->beam_nodes[beam_size-1];
        node2 = surface[i]->beam_nodes[beam_size-2];
        vector_beam = node2-node1;
        vector_beam.normalize();
        vector_perp = vector_beam.cross(vectorY);
        
        rmin = 1000; 
        rmax = 0;
        for (int j=0;j<beam_size;j++) { // get the min/max radii of the beam 
            node1 = surface[i]->beam_nodes[j];
            rcurrent = node1.norm();
            if (rcurrent < rmin)
                rmin = rcurrent;
            if (rcurrent > rmax)
                rmax = rcurrent;
        }
        
        for (int j=0;j<beam_size;j++) {
            node1 = surface[i]->beam_nodes[j];
            rcurrent = node1.norm();
            
            rfrac = (rcurrent-rmin)/(rmax-rmin);
            vector_perp_frac = vector_perp*rfrac ;
            
            surface[i]->rotate_section_about_arbitrary_point(center, vector_perp_frac, j, true);
        }
        
        surface[i]->compute_panel_components();
        
    }


// oscillation around beam -> change the angle of attack
/*
    int beam_size;
    double amplitude = 3.1415926535897932384626433832795 / 24 ;   // 7.5 degrees
    double omega = 23.565325441667691928903819677371; // RPS // omega = half the rotation frequency
    double sinus = sin(iteration*dt*omega);
    double motion = amplitude*sinus;
    
    vector3d node1,node2,vector_perp,vector_beam,center;
    for (int i=0;i<surface.size();i++){
        beam_size = surface[i]->beam_nodes.size();
        node1 = surface[i]->beam_nodes[beam_size-1];
        node2 = surface[i]->beam_nodes[beam_size-2];
        vector_beam = node2-node1;
        vector_beam.normalize();
        
        vector_beam *= motion ;
        
        for (int j=0;j<beam_size;j++) {
            center = surface[i]->beam_nodes[j];
            surface[i]->rotate_section_about_arbitrary_point(center, vector_beam, j, true);
        }
        surface[i]->compute_panel_components();
    }
*/
}

// global functions

// create petsc vec
void petsc_vec_create(Vec& vec, int size){

    VecCreate(PETSC_COMM_WORLD,&vec);
    VecSetSizes(vec,PETSC_DECIDE,size);
    VecSetFromOptions(vec);
    VecSet(vec,0.0);
}

// create petsc mat
void petsc_mat_create(Mat& mat, const int rows, const int cols){

    MatCreateSeqDense(PETSC_COMM_WORLD,rows,cols,NULL,&mat);
    MatSetFromOptions(mat);
    MatSetUp(mat);
    MatZeroEntries(mat);
}

// write petsc mat to a file readable by MATLAB
void WriteMat(Mat& mat,char const *name){

    PetscViewer viewer;
    char filename[64] = "Mat_"; char pfix[12] = ".m";
    strcat(filename,name); strcat(filename,pfix);    
    PetscObjectSetName((PetscObject)mat,name);
    
// ASCII Version
//    PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
//    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

// Binary Version // can be open in matlab with the function 
//    \petsc-3.10.2\share\petsc\matlab\PetscBinaryRead.m
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    
    MatView(mat,viewer);
    PetscViewerDestroy(&viewer);
}


// write petsc mat to a file readable by MATLAB
void WriteVec(Vec& vec,char const *name){

    PetscViewer viewer;
    char filename[64] = "Mat_"; char pfix[12] = ".m";
    strcat(filename,name); strcat(filename,pfix);    
    PetscObjectSetName((PetscObject)vec,name);
    
// ASCII Version
//    PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer);
//    PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

// Binary Version // can be open in matlab with the function 
//    \petsc-3.10.2\share\petsc\matlab\PetscBinaryRead.m
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);
    
    VecView(vec,viewer);
    PetscViewerDestroy(&viewer);
}