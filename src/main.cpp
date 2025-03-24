/**
 * @file main.cpp
 * @brief Main program for BSE (Bethe-Salpeter Equation) calculations
 *
 * This program performs various types of BSE calculations including:
 * - Chi0 (bare susceptibility): Calculates the non-interacting susceptibility
 * - RPA (Random Phase Approximation): Includes bubble diagrams with bare vertex
 * - BSE (Full Bethe-Salpeter Equation): Full vertex corrections included
 * - RRPA (Renormalized RPA): RPA with renormalized vertex
 * - SCL (Strong-Coupling Limit): https://www.physics.okayama-u.ac.jp/~otsuki/pdf/PhysRevB.99.165134.pdf
 *
 * The program reads input parameters and susceptibility data from files,
 * performs the requested calculation type, and outputs results to files.
 * It supports both local (momentum-integrated) and q-dependent calculations.
 */

#include "bse.hpp"
#include <string>
#include <sstream>
#include <complex>
#include <algorithm>
#include <map>
#include <iomanip>
#include "set_initial.hpp"
#include <fstream>

/**
 * @brief Debug utility to dump contents of a block matrix
 * @param block_mat Block matrix to dump
 * @param matsize Vector containing matrix block sizes
 */
void DumpMatrix(bse::block_matrix<std::complex<double> > block_mat, std::vector<int> matsize) {
  // Iterate through all block matrix elements
  for (int i = 0; i < matsize.size(); i++) {
    for (int j = 0; j < matsize.size(); j++) {
      // Only print if block exists
      if (block_mat.exists(i, j)) {
        auto bm_view = block_mat(i, j);
        // Print block indices and first element's real/imag parts
        std::cout << i << ", " << j << ", " << bm_view(0, 0).real() << ", " << bm_view(0, 0).imag() << std::endl;
      }
    }
  }
}

/**
 * @brief Structure containing input parameters and filenames
 */
struct Init {
  std::string FileNameX0_loc;    ///< Local bare susceptibility X0
  std::string FileNameX_loc;     ///< Local full susceptibility X
  std::string FileNameX0_q;      ///< q-dependent bare susceptibility X0
  std::string FileNameChi0_loc;  ///< Local bare susceptibility Chi0
  std::string FileNameChi_loc;   ///< Local full susceptibility Chi
  std::string FileNameMatrixInfo;///< Matrix size information
  /* For RPA */
  std::string FileNameGamma0;    ///< Bare vertex for RPA
  /* For SCL */
  std::string FileNamePhi;       ///< Full vertex
  std::string FileNamePhi_Sum;   ///< Summed full vertex
  std::string CalcType;          ///< Calculation type (chi0/rpa/bse/rrpa/scl)
  double beta;                   ///< Inverse temperature
  int suffix;                    ///< File suffix number
};

/**
 * @brief Helper function to set string parameter from input map
 * @param _info Map containing input parameters
 * @param _KW Keyword to look for
 * @param _fileName String parameter to set
 * @return true if successful, false if keyword not found
 */
bool SetString(std::map<std::string, std::string> _info, std::string _KW, std::string &_fileName) {
  // Find keyword in map
  std::map<std::string, std::string>::iterator itr;
  itr = _info.find(_KW);
  if (itr != _info.end()) {
    // Set string value if found
    _fileName = itr->second;
  } else {
    return false;
  }
  return true;
}

/**
 * @brief Helper function to set boolean flag from input map
 * @param _info Map containing input parameters  
 * @param _KW Keyword to look for
 * @param flag Boolean flag to set
 * @return true if successful, false if keyword not found or invalid value
 */
bool SetFlag(std::map<std::string, std::string> _info, std::string _KW, bool &flag) {
  std::map<std::string, std::string>::iterator itr;
  itr = _info.find(_KW);

  // Find keyword in map
  itr = _info.find(_KW);
  if (itr != _info.end()) {
    // Convert value to uppercase
    std::transform(itr->second.cbegin(), itr->second.cend(), itr->second.begin(), toupper);
    // Set flag based on TRUE/FALSE string
    if (itr->second == "TRUE") {
      flag = true;
    } else if (itr->second == "FALSE") {
      flag = false;
    } else {
      std::cerr << _KW << " is incorrect." << std::endl;
      return false;
    }
  }
  return true;
}

/**
 * @brief Helper function to set numeric value from input map
 * @tparam T Numeric type (int or double)
 * @param _info Map containing input parameters
 * @param _KW Keyword to look for  
 * @param _value Value to set
 * @return true if successful, false if keyword not found or invalid value
 */
template<typename T>
bool SetValue(std::map<std::string, std::string> _info, std::string _KW, T &_value) {
  // Find keyword in map
  std::map<std::string, std::string>::iterator itr;
  itr = _info.find(_KW);
  if (itr != _info.end()) {
    std::string svalue;
    svalue = itr->second;
    // Convert string to numeric value based on type
    if (typeid(T) == typeid(double))
      _value = std::stod(svalue);
    else if (typeid(T) == typeid(int))
      _value = std::stoi(svalue);
    else
      return false;
    return true;
  }
  return false;
}

/**
 * @brief Read input parameters from file
 * @param _InputFileName Name of input file
 * @param init Structure to store parameters
 * @return true if successful, false if error occurred
 */
bool GetInfo(
        const std::string _InputFileName,
        Init &init
) {
  // Map to store key-value pairs from input file
  std::map<std::string, std::string> info;
  
  // Open input file
  std::ifstream fin;
  fin.open(_InputFileName, std::ios::in);
  if (fin.fail()) {
    std::cerr << " Input file does not exist." << std::endl;
    return false;
  }

  // Read file line by line
  std::vector<std::string> reading_lines;
  std::string reading_line;
  while (!fin.eof()) {
    std::getline(fin, reading_line);
    reading_lines.push_back(reading_line);
  }

  // Parse each line into key-value pairs
  for (int i = 0; i < reading_lines.size(); i++) {
    std::istringstream stream(reading_lines[i]);
    std::string kwd, value;
    int iret = 0;
    while (std::getline(stream, reading_line, ' ')) {
      if (iret == 0) {
        // Convert keyword to uppercase
        std::transform(reading_line.cbegin(), reading_line.cend(), reading_line.begin(), toupper);
        kwd = reading_line;
      } else if (iret == 1) value = reading_line;
      iret++;
    }
    info.insert(std::pair<std::string, std::string>(kwd, value));
  }

  // Read required parameters
  if (!SetString(info, "X0_Q", init.FileNameX0_q)) {
    std::cout << "Fail to read X0_Q" << std::endl;
    return false;
  }
  if (!SetString(info, "X0_LOC", init.FileNameX0_loc)) {
    std::cout << "Fail to read X0_LOC" << std::endl;
    return false;
  }
  if (!SetString(info, "CALCTYPE", init.CalcType)) {
    std::cout << "Fail to read CALCTYPE" << std::endl;
    return false;
  }
  if (!SetValue(info, "SUFFIX", init.suffix)) {
    std::cout << "Fail to read SUFFIX" << std::endl;
    return false;
  }
  if (!SetString(info, "MATRIX_INFO", init.FileNameMatrixInfo)) {
    std::cout << "Fail to read MATRIX_INFO" << std::endl;
    return false;
  }
  if (!SetValue(info, "BETA", init.beta)) {
    std::cout << "Fail to read BETA" << std::endl;
    return false;
  }

  if (!SetString(info, "CHI0_LOC", init.FileNameChi0_loc)) {
    std::cout << "Fail to read CHI0_LOC" << std::endl;
    return false;
  }

  // Read calculation type specific parameters
  if (init.CalcType == "chi0"){
    // Chi0 mode - no additional parameters needed
  }
  else if (init.CalcType == "rpa") { //RPA mode
    if (!SetString(info, "GAMMA0", init.FileNameGamma0)) {
      std::cout << "Fail to read GAMMA0" << std::endl;
      return false;
    }
  }
  else if(init.CalcType == "bse") {
      if (!SetString(info, "X_LOC", init.FileNameX_loc)) {
        std::cout << "Fail to read X_LOC" << std::endl;
        return false;
      }
    if (!SetString(info, "CHI_LOC", init.FileNameChi_loc)) {
      std::cout << "Fail to read CHI_LOC" << std::endl;
      return false;
    }
  }

  if (init.CalcType == "rrpa") { //RRPA mode
    if (!SetString(info, "CHI_LOC", init.FileNameChi_loc)) {
      std::cout << "Fail to read CHI_LOC" << std::endl;
      return false;
    }
  } else if (init.CalcType == "scl") {
    if (!SetString(info, "PHI", init.FileNamePhi)) {
      std::cout << "Fail to read Phi" << std::endl;
      return false;
    }
    if (!SetString(info, "PHI_SUM", init.FileNamePhi_Sum)) {
      std::cout << "Fail to read Phi_SUM" << std::endl;
      return false;
    }
  }
  return true;
}

/**
 * @brief Save matrix data to file
 * @param filename Output filename
 * @param bm Block matrix to save
 * @param matrix_info Vector of matrix block sizes
 *
 * Saves matrix of type C where:
 * - block: [block], inner: [in]
 * Saves with superblock index i=(block1, in1), j=(block2, in2)
 */
void SaveData(std::string const &filename, bse::block_matrix<std::complex<double> > &bm, std::vector<int> &matrix_info){
  std::cout << "Save " << filename << std::endl;
  std::ofstream fout(filename);
  // Set maximum precision for floating point output
  fout << std::setprecision(std::numeric_limits<double>::max_digits10);
  
  // Track offsets for superblock indices
  int offset1 = 0;
  for (int b1 = 0; b1 < matrix_info.size(); b1++) {
    int offset2 = 0;
    for (int b2 = 0; b2 < matrix_info.size(); b2++) {
      if (bm.exists(b1, b2)) {
        auto bm_view = bm(b1, b2);
        // Write each element with superblock indices
        for (int in1=0; in1<matrix_info[b1]; in1++){
          for (int in2=0; in2<matrix_info[b2]; in2++){
            fout << in1+offset1 << ", " << in2+offset2 << ", "
                 << bm_view(in1, in2).real() << ", " << bm_view(in1, in2).imag()
                 << std::endl;
          }
        }
      }
      offset2 += matrix_info[b2];
    }
    offset1 += matrix_info[b1];
  }
  fout.close();
}

/**
 * @brief Print matrix size information
 * @param matrix_info Vector of matrix block sizes
 * @param str Description string
 */
void print_matrix_info(std::vector<int> &matrix_info, const std::string &str) {
  std::cout << "matrix_info  '" << str << "'" << std::endl;
  std::cout << "| n_blocks = " << matrix_info.size() << std::endl;
  std::cout << "|";
  for(int i=0; i<matrix_info.size(); i++){
    std::cout << " " << matrix_info[i];
  }
  std::cout << std::endl;
}

/**
 * @brief Main program entry point
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return 0 on success, non-zero on error
 */
int main(int argc, char *argv[]) {
  // Read Information from dat file
  int iret = 0;
  setInitial SI;
  Init init;

  // Parse command line arguments
  if (argc > 1) {
    if (argc != 2) iret = -1;
    else {
      if (GetInfo(argv[1], init) != true) {
        std::cout << "fail to read input file." << std::endl;
        exit(-1);
      }

    }
  }
  if (iret != 0) {
    std::cout << "too few arguments." << std::endl;
    exit(-1);
  }

  // Get numbers of inner dims without Matsubara
  //   innerDims[block] = norb*norb
  std::vector<int> innerDims = SI.GetMatrixInfo(init.FileNameMatrixInfo);
  print_matrix_info(innerDims, "innerDims");

  // Get MatsubaraSize from input files
  std::vector<std::string> sFileNames;  // Files to read
  sFileNames.push_back(init.FileNameX0_loc);
  sFileNames.push_back(init.FileNameX0_q);
  if (init.CalcType == "bse"){
    sFileNames.push_back(init.FileNameX_loc);
  }
  int nMatsubara = SI.GetMatsubaraSize(sFileNames);
  std::cout << "nMatsubara = " << nMatsubara << std::endl;

  // Make matrix_info arrays for different matrix types
  // matrix-A: block: [block], inner: [in, iw]
  std::vector<int> matrix_info_A(innerDims);
  for(auto &dim : matrix_info_A){
    dim *= nMatsubara;
  }
  
  // matrix-B: block: [iw, block], inner: [in]
  std::vector<int> matrix_info_B;
  for(int i=0; i<nMatsubara; i++){
    matrix_info_B.insert(matrix_info_B.end(), innerDims.begin(), innerDims.end());
  }
  
  // matrix-C: block: [block], inner: [in]
  std::vector<int> matrix_info_C(innerDims);

  // Print matrix info for A, B, C
  print_matrix_info(matrix_info_A, "matrix_info_A");
  print_matrix_info(matrix_info_B, "matrix_info_B");
  print_matrix_info(matrix_info_C, "matrix_info_C");

  std::cout << "=====================================================" << std::endl;
  std::cout << "===============Start BSE calculation=================" << std::endl;

  // Initialize matrices
  bse::block_matrix<std::complex<double> > X0_loc(matrix_info_B);
  bse::block_matrix<std::complex<double> > X0_q(matrix_info_B);
  bse::block_matrix<std::complex<double> > X_loc(matrix_info_A);

  bool verbose = true;  // verbose output on memory info
  // Initialize BSE solver
  bsesol::bse_solver<std::complex<double> > bse(init.beta, matrix_info_A, matrix_info_B, matrix_info_C, verbose);

  // Read input data
  std::cout << "Read X0_q.\n";
  SI.ReadData(init.FileNameX0_q, matrix_info_B, X0_q);
  bse.SetX0q(X0_q);

  std::cout << "Read X0_loc.\n";
  SI.ReadData(init.FileNameX0_loc, matrix_info_B, X0_loc);
  bse.SetX0Loc(X0_loc);

  std::cout << "Read Chi0_loc.\n";
  bse::block_matrix<std::complex<double> > chi0_loc_sum(matrix_info_C);
  SI.ReadData(init.FileNameChi0_loc, matrix_info_C, chi0_loc_sum);
  bse.SetChi0Loc(chi0_loc_sum);

  // Calculate Chi0_q
  bse::block_matrix<std::complex<double> > X0_q_sum;
  std::cout << "Start: Calc X0_q_sum.\n";
  bse.CalcChi0q();
  std::cout << "End  : Calc X0_q_sum.\n";

  auto chi0_q_sum = bse.GetChi0q();
  SaveData("chi0_q_"+std::to_string(init.suffix)+".tmp", chi0_q_sum, matrix_info_C);

  // Perform calculation based on type
  if (init.CalcType == "bse") { //BSE mode
    std::cout << "Read X_loc.\n";
    SI.ReadData(init.FileNameX_loc, matrix_info_A, X_loc);
    bse.SetXLoc(X_loc);

    std::cout << "Read chi_loc.\n";
    bse::block_matrix<std::complex<double> > chi_loc_sum(matrix_info_C);
    SI.ReadData(init.FileNameChi_loc, matrix_info_C, chi_loc_sum);
    bse.SetChiLoc(chi_loc_sum);

    std::cout << "Start: Calc chi_q_bse.\n";
    bse.CalcChiq_BSE();
    std::cout << "End  : Calc chi_q_bse.\n";

    // Save results
    auto chi_q_sum = bse.GetChiq_BSE();
    SaveData("chi_q_"+std::to_string(init.suffix)+".tmp", chi_q_sum, matrix_info_C);

    auto Iq = bse.GetIq();
    SaveData("Iq_"+std::to_string(init.suffix)+".tmp", Iq, matrix_info_C);

  } else if (init.CalcType == "rpa") { //RPA mode
    std::cout << "Read Gamma0.\n";
    std::vector<int> tmp_matrix = {0, 0};//N*N*N* ... * N = N^{Norb^2*Nspin^2}
    if (SI.GetMatrixInfoCore(init.FileNameGamma0, tmp_matrix) != true) {
      exit(-1);
    }
    bse::block_matrix<std::complex<double> > gamma0(matrix_info_C);
    SI.ReadData(init.FileNameGamma0, matrix_info_C, gamma0);
    bse.SetGamma0(gamma0);

    std::cout << "Calculate chi_q_rpa.\n";
    bse.CalcChiq_RPA();

    auto chi_q_rpa = bse.GetChiq_RPA();
    SaveData("chi_q_rpa_"+std::to_string(init.suffix)+".tmp", chi_q_rpa, matrix_info_C);

  } else if (init.CalcType == "rrpa") {
    std::cout << "Read Chi_loc.\n";
    std::vector<int> tmp_matrix = {0, 0};//N*N*N* ... * N = N^{Norb^2*Nspin^2}
    if (SI.GetMatrixInfoCore(init.FileNameChi_loc, tmp_matrix) != true) {
      exit(-1);
    }
    bse::block_matrix<std::complex<double> > chi_loc_sum(matrix_info_C);
    SI.ReadData(init.FileNameChi_loc, matrix_info_C, chi_loc_sum);
    bse.SetChiLoc(chi_loc_sum);

    bse::block_matrix<std::complex<double> > chi0_loc_sum(matrix_info_C);
    SI.ReadData(init.FileNameChi0_loc, matrix_info_C, chi0_loc_sum);
    bse.SetChi0Loc(chi0_loc_sum);

    std::cout << "Calculate chi_q_rrpa.\n";
    bse.CalcChiq_RRPA();

    auto chi_q_rrpa = bse.GetChiq_RRPA();
    SaveData("chi_q_rrpa_"+std::to_string(init.suffix)+".tmp", chi_q_rrpa, matrix_info_C);

  } else if (init.CalcType == "scl") {
    bse::block_matrix<std::complex<double> > Phi(matrix_info_B);
    std::cout << "Read Phi.\n";
    SI.ReadData(init.FileNamePhi, matrix_info_B, Phi);
    bse.SetPhi(Phi);

    std::cout << "Read Phi_Sum.\n";
    std::vector<int> tmp_matrix = {0, 0};//N*N*N* ... * N = N^{Norb^2*Nspin^2}
    if (SI.GetMatrixInfoCore(init.FileNamePhi_Sum, tmp_matrix) != true) {
      exit(-1);
    }
    bse::block_matrix<std::complex<double> > Phi_Sum(matrix_info_C);
    SI.ReadData(init.FileNamePhi_Sum, matrix_info_C, Phi_Sum);
    bse.SetPhi_Sum(Phi_Sum);

    std::cout << "Start: Calc chi_q_scl.\n";
    bse.CalcChiq_SCL();
    std::cout << "End  : Calc chi_q_scl.\n";

    auto chi_q_scl = bse.GetChiq_SCL();
    SaveData("chi_q_scl_"+std::to_string(init.suffix)+".tmp", chi_q_scl, matrix_info_C);

    auto Iq = bse.GetIq();
    SaveData("Iq_scl_"+std::to_string(init.suffix)+".tmp", Iq, matrix_info_C);
  }
  bse.Summary();
  std::cout << "=================End BSE calculation=================" << std::endl;
  std::cout << "=====================================================" << std::endl;
  return 0;
}
