#include <lasSparse.h>
#include <lasCSR.h>
#include <lasSparskit.h>
#include <iostream>
#include <fstream>
typedef las::csrMat skMat;
int main() {
  std::ifstream in("/fasttmp/mersoj/develop/las/test/small.mtx");
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
  std::cout<<"nnz: "<<nnz<<std::endl;
  std::cout<<"num rows/cols"<<csr->getNumRows()<<" "<<csr->getNumCols()<<std::endl;
  int * rows = csr->getRows();
  int * cols = csr->getCols();
  double * vals = m->getVals();
  std::cout<<"Rows: ";
  for(int i=0; i<csr->getNumRows()+1; ++i) {
    std::cout<<rows[i]<<" ";
  }
  std::cout<<std::endl;
  std::cout<<"Cols: ";
  for(int i=0; i<nnz; ++i) {
    std::cout<<cols[i]<<" ";
  }
  std::cout<<std::endl;
  std::cout<<"Vals: ";
  for(int i=0; i<nnz; ++i) {
    std::cout<<vals[i]<<" ";
  }
  std::cout<<std::endl;
  las::printSparskitMat(std::cout, readMat, las::PrintType::full, false);
  std::cout<<"print mmarket"<<std::endl;
  las::printSparskitMat(std::cout, readMat, las::PrintType::mmarket, false);
  std::cout<<"done"<<std::endl;
}
