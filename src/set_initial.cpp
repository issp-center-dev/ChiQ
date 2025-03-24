#include <iostream>
#include <fstream>
#include "set_initial.hpp"

/**
 * Get matrix dimensions from multiple input files
 *
 * @param sFileName_chi0_loc Filename for local bare susceptibility
 * @param sFileName_chi_loc Filename for local full susceptibility  
 * @param sFileName_ch0_q Filename for q-dependent bare susceptibility
 * @return Vector containing matrix dimensions (block size, inner dimension)
 */
std::vector<int> setInitial::GetMatrixInfo
        (
                std::string sFileName_chi0_loc,
                std::string sFileName_chi_loc,
                std::string sFileName_ch0_q
        ) {
  std::vector<int> matInfo = {0, 0};
  std::vector<int> tmp_matInfo = {0, 0};
  // Get dimensions from each file and take minimum
  GetMatrixInfoCore(sFileName_chi0_loc, matInfo);
  GetMatrixInfoCore(sFileName_chi_loc, tmp_matInfo);
  for (int i = 0; i < matInfo.size(); i++) {
    if (tmp_matInfo[i] < matInfo[i]) matInfo[i] = tmp_matInfo[i];
  }
  GetMatrixInfoCore(sFileName_ch0_q, tmp_matInfo);
  for (int i = 0; i < matInfo.size(); i++) {
    if (tmp_matInfo[i] < matInfo[i]) matInfo[i] = tmp_matInfo[i];
  }
  std::vector<int> iMatSize(matInfo[0], matInfo[1]);
  return iMatSize;
}

/**
 * Get matrix dimensions from two input files
 * 
 * @param sFileName_chi0_loc Filename for local bare susceptibility
 * @param sFileName_ch0_q Filename for q-dependent bare susceptibility
 * @return Vector containing matrix dimensions (block size, inner dimension)
 */
std::vector<int> setInitial::GetMatrixInfo
        (
                std::string sFileName_chi0_loc,
                std::string sFileName_ch0_q
        ) {
  std::vector<int> matInfo = {0, 0};
  std::vector<int> tmp_matInfo = {0, 0};
  // Get dimensions from each file and take minimum
  GetMatrixInfoCore(sFileName_chi0_loc, matInfo);
  GetMatrixInfoCore(sFileName_ch0_q, tmp_matInfo);
  for (int i = 0; i < matInfo.size(); i++) {
    if (tmp_matInfo[i] < matInfo[i]) matInfo[i] = tmp_matInfo[i];
  }
  std::vector<int> iMatSize(matInfo[0], matInfo[1]);
  return iMatSize;
}

/**
 * Core function to read matrix dimensions from a single file
 *
 * @param sFileName Input filename
 * @param matrix_info Vector to store dimensions (must be size 2)
 * @return true if successful, false on error
 */
bool setInitial::GetMatrixInfoCore
        (
                std::string sFileName,
                std::vector<int> &matrix_info
        ) {
  int iBlockDim, imatDim;
  if (matrix_info.size() != 2) {
    std::cout << "Error: the size of matrix_info vector is incorrect." << std::endl;
    return false;
  }
  std::ifstream fin(sFileName);
  if (fin.fail()) {
    std::cout << "Error: Read file:" << sFileName << std::endl;
    return false;
  }
  // Read block dimension and matrix dimension
  fin >> iBlockDim >> imatDim;
  matrix_info[0] = iBlockDim;
  matrix_info[1] = imatDim;
  fin.close();
  return true;
}

/**
 * Get Matsubara frequency size from multiple files
 *
 * @param sFileNames Vector of input filenames
 * @return Number of Matsubara frequencies (must be consistent across files)
 */
int setInitial::GetMatsubaraSize
        (
                std::vector<std::string> const& sFileNames
        ) {
  std::vector<int> nMatsubaras;
  // Get Matsubara size from each file
  for (const auto &sFileName : sFileNames){
    nMatsubaras.push_back(GetMatsubaraSizeCore(sFileName));
  }

  // Check if all components are equivalent
  int nMatsubara = nMatsubaras[0];
  for (const auto &n : nMatsubaras){
    if (n != nMatsubara){
      std::cerr << "Number of fermionic Matsubara frequencies are inconsistent" << std::endl;
      exit(1);
    }
  }
  return nMatsubara;

  // Find minimum
  // int min_nMatsubara = *std::min_element(nMatsubaras.begin(), nMatsubaras.end());
  // return min_nMatsubara;
}

/**
 * Core function to read Matsubara frequency size from a single file
 *
 * @param sFileName Input filename
 * @return Number of Matsubara frequencies
 */
int setInitial::GetMatsubaraSizeCore
        (
                std::string sFileName
        ) {
  std::ifstream fin(sFileName);
  if (fin.fail()) {
    std::cout << "Error: Read file:" << sFileName << std::endl;
    return false;
  }
  int nMatsubara;
  fin >> nMatsubara;
  fin.close();
  return nMatsubara;
}

/**
 * Read matrix dimensions from a file containing block structure
 *
 * @param sFileName Input filename
 * @return Vector containing inner dimensions for each block
 */
std::vector<int> setInitial::GetMatrixInfo(std::string sFileName) {
  std::ifstream fin(sFileName);
  if (fin.fail()) {
    std::cout << "Error: Read file:" << sFileName << std::endl;
    exit(1);
  }
  // Read number of blocks
  int iBlockDim;
  fin >> iBlockDim;
  // Read inner dimension for each block
  std::vector<int> imatDim(iBlockDim);
  for(int i=0; i<iBlockDim; i++){
    fin >> imatDim[i];
  }
  fin.close();
  return imatDim;
}
