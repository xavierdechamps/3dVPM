#include "parameters.hpp"

double Parameters :: inversion_tolerance = 1e-12;

double Parameters :: farfield_factor = 10.0;

double Parameters :: trailing_edge_wake_shed_factor = 0.25;

bool Parameters :: unsteady_problem = false;

double Parameters :: static_wake_length = 1.0;

bool Parameters :: static_wake = false;

bool Parameters :: use_vortex_core_model = false;

int Parameters :: max_number_iterations_pressure_TE = 100;
double Parameters :: max_error_pressure_TE = 0.1 ;

bool Parameters :: switch_off_root_vortex = false;
bool Parameters :: switch_off_tip_vortex = false;