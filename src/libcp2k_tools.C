#include "libcp2k_tools.H"

void intit_transport_parameters_from_cp2k(cp2k_transport_parameters& cp2k_transport_params, transport_parameters* transport_params)
{

    transport_params->n_occ                      = cp2k_transport_params.n_occ;
    transport_params->n_atoms                    = cp2k_transport_params.n_atoms;
    transport_params->evoltfactor                = cp2k_transport_params.evoltfactor;
    transport_params->method                     = cp2k_transport_params.method;
    transport_params->bandwidth                  = cp2k_transport_params.bandwidth;
    transport_params->n_cells                    = cp2k_transport_params.n_cells;
    transport_params->n_abscissae                = cp2k_transport_params.n_abscissae;
    transport_params->n_kpoint                   = cp2k_transport_params.n_kpoint;
    transport_params->num_interval               = cp2k_transport_params.num_interval;
    transport_params->num_contacts               = cp2k_transport_params.num_contacts;
    transport_params->ndof                       = cp2k_transport_params.ndof;
    transport_params->tasks_per_point            = cp2k_transport_params.tasks_per_point;
    transport_params->cores_per_node             = cp2k_transport_params.cores_per_node;
    transport_params->colzero_threshold          = cp2k_transport_params.colzero_threshold;
    transport_params->eps_limit                  = cp2k_transport_params.eps_limit;
    transport_params->eps_decay                  = cp2k_transport_params.eps_decay;
    transport_params->eps_singularity_curvatures = cp2k_transport_params.eps_singularity_curvatures;
    transport_params->eps_mu                     = cp2k_transport_params.eps_mu;
    transport_params->eps_eigval_degen           = cp2k_transport_params.eps_eigval_degen;
    transport_params->energy_interval            = cp2k_transport_params.energy_interval;
    transport_params->min_interval               = cp2k_transport_params.min_interval;
    transport_params->temperature                = cp2k_transport_params.temperature;

}
