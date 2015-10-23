#include "WriteMatrix.H"

void write_matrix(TCSR<double> *Overlap,TCSR<double> *KohnSham,transport_parameters *transport_params) {
    int iam;MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    KohnSham->change_findx(1);
    Overlap->change_findx(1);
    if (transport_params->ndof) {
        int bigblocksize=Overlap->size_tot/3;
        int cutblocksize=transport_params->ndof;
        for (int i=0;i<3;i++) {
            TCSR<double> *CutG = new TCSR<double>(KohnSham,bigblocksize,cutblocksize,i*bigblocksize,cutblocksize);
            TCSR<double> *Cut = new TCSR<double>(CutG,0,MPI_COMM_WORLD);
            delete CutG;
            if (!iam) {
                stringstream mysstream;
                mysstream << "H_" << i+3 << ".bin";
                Cut->write_CSR_bin(mysstream.str().c_str());
            }
            delete Cut;
        }
        for (int i=0;i<3;i++) {
            TCSR<double> *CutG = new TCSR<double>(Overlap,bigblocksize,cutblocksize,i*bigblocksize,cutblocksize);
            TCSR<double> *Cut = new TCSR<double>(CutG,0,MPI_COMM_WORLD);
            delete CutG;
            if (!iam) {
                stringstream mysstream;
                mysstream << "S_" << i+3 << ".bin";
                Cut->write_CSR_bin(mysstream.str().c_str());
            }
            delete Cut;
        }
    } else {
        TCSR<double> *KohnShamCollect = new TCSR<double>(KohnSham,0,MPI_COMM_WORLD);
        if (!iam) {
            KohnShamCollect->write_CSR("KohnSham");
            KohnShamCollect->removepbc(transport_params->bandwidth,KohnShamCollect->size_tot/transport_params->n_cells);
            KohnShamCollect->write_CSR_bin("H_4.bin");
        }
        delete KohnShamCollect;
        TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,0,MPI_COMM_WORLD);
        if (!iam) {
            OverlapCollect->write_CSR("Overlap");
            OverlapCollect->removepbc(transport_params->bandwidth,OverlapCollect->size_tot/transport_params->n_cells);
            OverlapCollect->write_CSR_bin("S_4.bin");
        }
        delete OverlapCollect;
        if (!iam) {
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
    }
    MPI_Barrier(MPI_COMM_WORLD);
    exit(0);
}
