void makeCovMat()
{
    // TFile* fout = new TFile("cov.root","RECREATE");
    // int ndim = 20;
    // TMatrixDSym cov_matrix(ndim);
    // for(int r = 0; r < ndim; ++r)
    // {
    //     for(int c = 0; c < ndim; ++c)
    //     {
    //         cov_matrix[r][c] = r==c? 1 : 0;
    //     }
    // }
    // cov_matrix.Write("cov_matrix");
    // fout->Close();

    TFile* fout = new TFile("Lalpha_beta_prior.root","RECREATE");
    int ndim = 2;
    TMatrixDSym cov_matrix(ndim);
    cov_matrix[0][0] = 1.e12;
    cov_matrix[1][1] = 0.015*0.015;
    cov_matrix.Write("cov_matrix");

    TMatrixDSym cov_matrix2(1);
    cov_matrix2[0][0] = 100*100;
    cov_matrix2.Write("cov_matrix_Lalpha");

    fout->Close();

}