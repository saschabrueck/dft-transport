#include "libcp2k_tools.H"

/**
 *   \brief initializes a transport_parameters type from a cp2k_transport_parameters type
 *   \author Mohammad Hossein Bani-Hashemian
 */
void intit_transport_parameters_from_cp2k(cp2k_transport_parameters& cp2k_transport_params, transport_parameters* transport_params)
{

    transport_params->n_occ                      = cp2k_transport_params.n_occ;
    transport_params->n_atoms                    = cp2k_transport_params.n_atoms;
    transport_params->energy_diff                = cp2k_transport_params.energy_diff;
    transport_params->evoltfactor                = cp2k_transport_params.evoltfactor;
    transport_params->method                     = cp2k_transport_params.method;
    transport_params->injection_method           = cp2k_transport_params.injection_method;
    transport_params->linear_solver              = cp2k_transport_params.linear_solver;
    transport_params->n_abscissae                = cp2k_transport_params.n_abscissae;
    transport_params->n_kpoint                   = cp2k_transport_params.n_kpoint;
    transport_params->num_interval               = cp2k_transport_params.num_interval;
    transport_params->num_contacts               = cp2k_transport_params.num_contacts;
    transport_params->tasks_per_point            = cp2k_transport_params.tasks_per_point;
    transport_params->colzero_threshold          = cp2k_transport_params.colzero_threshold;
    transport_params->eps_limit                  = cp2k_transport_params.eps_limit;
    transport_params->eps_decay                  = cp2k_transport_params.eps_decay;
    transport_params->eps_singularity_curvatures = cp2k_transport_params.eps_singularity_curvatures;
    transport_params->eps_mu                     = cp2k_transport_params.eps_mu;
    transport_params->eps_eigval_degen           = cp2k_transport_params.eps_eigval_degen;
    transport_params->energy_interval            = cp2k_transport_params.energy_interval;
    transport_params->min_interval               = cp2k_transport_params.min_interval;
    transport_params->temperature                = cp2k_transport_params.temperature;

    c_icopy(2, cp2k_transport_params.cutout, 1, transport_params->cutout, 1);
    transport_params->contacts_data = new int[cp2k_transport_params.num_contacts*5];
    c_icopy(cp2k_transport_params.num_contacts*5, cp2k_transport_params.contacts_data, 1, transport_params->contacts_data, 1);

}
