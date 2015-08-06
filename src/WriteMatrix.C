#include "WriteMatrix.H"

void write_matrix(TCSR<double> *Overlap,TCSR<double> *KohnSham,transport_parameters *transport_params) {
    TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,0,MPI_COMM_WORLD);
    TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,0,MPI_COMM_WORLD);
    int iam;MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    if (!iam) {
//KohnShamCollect->remove_thr(1.0E-6);
//OverlapCollect->remove_thr(1.0E-6);
        KohnShamCollect->change_findx(1);
        OverlapCollect->change_findx(1);
        KohnShamCollect->write_CSR("KohnSham");
        OverlapCollect->write_CSR("Overlap");
        KohnShamCollect->removepbc(transport_params->bandwidth,OverlapCollect->size_tot/transport_params->n_cells);
        OverlapCollect->removepbc(transport_params->bandwidth,OverlapCollect->size_tot/transport_params->n_cells);
        KohnShamCollect->write_CSR_bin("H_4.bin");
        OverlapCollect->write_CSR_bin("S_4.bin");
        ofstream paramoutfile("TransportParams");
        paramoutfile << transport_params->method << endl;
        paramoutfile << transport_params->bandwidth << endl;
        paramoutfile << transport_params->n_cells << endl;
        paramoutfile << transport_params->n_occ << endl;
        paramoutfile << transport_params->n_atoms << endl;
        paramoutfile << transport_params->n_abscissae << endl;
        paramoutfile << transport_params->n_kpoint << endl;
        paramoutfile << transport_params->num_interval << endl;
        paramoutfile << transport_params->num_contacts << endl;
        paramoutfile << transport_params->ndof << endl;
        paramoutfile << transport_params->tasks_per_point << endl;
        paramoutfile << transport_params->cores_per_node << endl;
        paramoutfile << transport_params->evoltfactor << endl;
        paramoutfile << transport_params->colzero_threshold << endl;
        paramoutfile << transport_params->eps_limit << endl;
        paramoutfile << transport_params->eps_decay << endl;
        paramoutfile << transport_params->eps_singularity_curvatures << endl;
        paramoutfile << transport_params->eps_mu << endl;
        paramoutfile << transport_params->eps_eigval_degen << endl;
        paramoutfile << transport_params->energy_interval << endl;
        paramoutfile << transport_params->min_interval << endl;
        paramoutfile << transport_params->temperature << endl;
        paramoutfile.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
}
