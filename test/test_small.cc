// This test verifies that reading matrix market data is occuring properly.
// We also write the data just to confirm that there are no failures on writing.
#include <lasCSR.h>
#include <lasSparse.h>
#include <lasSparskit.h>
#include <fstream>
#include <iostream>
typedef las::csrMat skMat;
skMat * readAndPrint(char * fileName)
{
  std::ifstream in(fileName);
  if (!in.is_open())
  {
    std::cerr << "Could not open small.mtx for reading!" << std::endl;
    std::abort();
  }
  las::Mat * readMat = las::readSparskitMat(in, las::PrintType::mmarket);
  in.close();
  skMat * m = getSparskitMatrix(readMat);
  las::CSR * csr = m->getCSR();
  int nnz = csr->getNumNonzero();
  std::cout << "nnz: " << nnz << std::endl;
  std::cout << "num rows/cols" << csr->getNumRows() << " " << csr->getNumCols()
            << std::endl;
  int * rows = csr->getRows();
  int * cols = csr->getCols();
  double * vals = m->getVals();
  std::cout << "Rows: ";
  for (int i = 0; i < csr->getNumRows() + 1; ++i)
  {
    std::cout << rows[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Cols: ";
  for (int i = 0; i < nnz; ++i)
  {
    std::cout << cols[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "Vals: ";
  for (int i = 0; i < nnz; ++i)
  {
    std::cout << vals[i] << " ";
  }
  std::cout << std::endl;
  las::printSparskitMat(std::cout, readMat, las::PrintType::full, false);
  std::cout << "print mmarket" << std::endl;
  std::ofstream out("out.dat");
  las::printSparskitMat(std::cout, readMat, las::PrintType::mmarket, false);
  out.close();
  std::cout << "done" << std::endl;
  las::LasCreateMat * mb = las::getMatBuilder<las::sparskit>(0);
  // mb->destroy(readMat);
  // las::destroySparsity<las::CSR *>((las::Sparsity *)csr);
  return m;
}
int main()
{
  int rows_true[4] = {1, 3, 5, 6};
  int cols_true[5] = {1, 3, 2, 3, 3};
  double vals_true[5] = {1, 5, 3, 7, 6};
  skMat * m = readAndPrint("./small.mtx");
  las::CSR * csr = m->getCSR();
  int * rows = csr->getRows();
  int * cols = csr->getCols();
  double * vals = m->getVals();
  for (int i = 0; i < 4; ++i)
  {
    assert(rows_true[i] == rows[i]);
  }
  for (int i = 0; i < 5; ++i)
  {
    assert(cols_true[i] == cols[i]);
    assert(fabs(vals_true[i] - vals[i]) < 1E-10);
  }
  int rows_true2[5] = {1, 3, 4, 6, 7};
  int cols_true2[6] = {1, 3, 2, 1, 3, 4};
  double vals_true2[6] = {1, 5, 3, 5, 6, 1.23};
  skMat * m2 =
      readAndPrint("./small_symm_matrix.mtx");
  csr = m2->getCSR();
  rows = csr->getRows();
  cols = csr->getCols();
  vals = m2->getVals();
  for (int i = 0; i < 5; ++i)
  {
    assert(rows_true2[i] == rows[i]);
  }
  for (int i = 0; i < 6; ++i)
  {
    assert(cols_true2[i] == cols[i]);
    assert(fabs(vals_true2[i] - vals[i]) < 1E-10);
  }
  // release the matrix builder when everything is done to avoid
  // memory leaks
  las::LasCreateMat * mb = las::getMatBuilder<las::sparskit>(0);
  delete mb;
  mb = NULL;
}
