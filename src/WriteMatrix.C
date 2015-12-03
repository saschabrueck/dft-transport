#include "WriteMatrix.H"

void write_matrix(TCSR<double> *Overlap,TCSR<double> *KohnSham,int cutblocksize,int bandwidth,int ndof) {
    int iam;MPI_Comm_rank(MPI_COMM_WORLD,&iam);
    KohnSham->change_findx(1);
    Overlap->change_findx(1);
    if (cutblocksize) {
        int bigblocksize=Overlap->size_tot/3;
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
            KohnShamCollect->removepbc(bandwidth,ndof);
            KohnShamCollect->write_CSR_bin("H_4.bin");
        }
        delete KohnShamCollect;
        TCSR<double> *OverlapCollect = new TCSR<double>(Overlap,0,MPI_COMM_WORLD);
        if (!iam) {
            OverlapCollect->write_CSR("Overlap");
            OverlapCollect->removepbc(bandwidth,ndof);
            OverlapCollect->write_CSR_bin("S_4.bin");
        }
        delete OverlapCollect;
    }
}
