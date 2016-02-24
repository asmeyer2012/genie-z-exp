{
#include <TMatrixD.h>

  int nP = 4;
  TFile *fMat = new TFile("tmat.out.root","recreate");

  // covariance matrix
  double tDat[nP][nP] =
    {{ 0.0174208,   0.00609007, -0.113303,  0.0970835},
     { 0.00609007,  0.347222,   -0.369704, -0.397993},
     {-0.113303,   -0.369704,    1.074980, -0.336495},
     { 0.0970835,  -0.397993,   -0.336495,  2.240440}};

  // flatten matrix into 1-D array
  double tFlt[nP* nP] = {0.};
  for (int i=0;i<nP;i++) {
  for (int j=0;j<nP;j++) {
    if (tFlt[j*nP+i] != tDat[i][j] && i > j)
    {
      std::cout<< "Matrix needs to be symmetric!" << std::endl;
      std::cout<< i<<","<<j<<": "<<tFlt[j*nP+i] <<","<< tDat[i][j]<<std::endl;
      exit(1);
    }
    tFlt[i*nP+j] = tDat[i][j];
  }}
  
  // make into TMatrixD object and write to file
  TMatrixD tMat(nP,nP);
  tMat.SetMatrixArray(tFlt);
  tMat->Write("tMat",TObject::kOverwrite);

  fMat->Close();
  exit(0);
}
