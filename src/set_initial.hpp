#pragma once

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <complex>
#include <Eigen/Dense>
#include "./bse.hpp"

/**
 * @brief Class for initializing and reading input data
 *
 * This class handles reading matrix data from files and initializing
 * various parameters needed for BSE calculations.
 */
class setInitial {
public:
  /**
   * Default constructor
   * If an object is constructed via default constructor,
   * it must be initialized via copy before any use.
   */
  setInitial() {}

  /**
   * @brief Read diagonal matrix data from file
   * @param sFileName Input filename
   * @param matrix_info Vector containing matrix dimensions
   * @param matrix Block matrix to store the data
   * @return true if successful, false on error
   */
  bool ReadDataDiagonal(
          std::string sFileName,
          std::vector<int> &matrix_info,
          bse::block_matrix<std::complex<double> > &matrix
  ) {
    int itmp;
    std::ifstream fin(sFileName);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> bm_insert(matrix_info[0], matrix_info[0]);
    if (fin.fail()) {
      std::cerr << "Error: Read file:" << sFileName << std::endl;
      // return false;
      exit(-1);
    }
    //ignore oneline
    std::string str;
    int iMatrixDim;
    fin >> itmp >> iMatrixDim;
    getline(fin, str);
    // Initialize temporary matrix
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> bm(iMatrixDim, iMatrixDim);
    for(int i=0; i<iMatrixDim; i++) {
      for (int j = 0; j < iMatrixDim; j++) {
        bm(i, j) = 0;
      }
    }
    int idx = 0;
    int iblock1, iblock2;

    // Read data line by line
    while (getline(fin, str)) {
      std::string token;
      std::istringstream stream(str);
      if (idx % (iMatrixDim + 1) == 0) {
        // Read block indices
        int ivalue[2];
        int tmpidx = 0;
        while (getline(stream, token, ',')) {
          ivalue[tmpidx] = stoi(token);
          tmpidx++;
        }
        iblock1 = ivalue[0];
        iblock2 = ivalue[1];
        idx++;
      } else {
        // Read matrix elements
        double cvalue[2];
        int tmpidx = 0;
        while (getline(stream, token, ',')) {
          cvalue[tmpidx] = stod(token);
          tmpidx++;
        }
        bm(idx - 1, idx - 1) = std::complex<double>(cvalue[0], cvalue[1]);
        idx++;
      }

      // When block is complete, assign to matrix
      if (idx % (iMatrixDim + 1) == 0) {
        if (iMatrixDim > matrix_info[0]) {
          // Extract center block if input is larger
          int iLeftCorner = (iMatrixDim - matrix_info[0]) / 2;
          bm_insert = bm.block(iLeftCorner, iLeftCorner, matrix_info[0], matrix_info[0]);

        } else {
          if (iMatrixDim != matrix_info[0]) {
            return false;
          }
          bm_insert = bm;
        }
        matrix.assign(iblock1, iblock2, bm_insert);
        idx = 0;
      }
    }

    fin.close();
    return true;
  }

  /**
   * @brief Read full matrix data from file
   * @param sFileName Input filename
   * @param matrix_info Vector containing matrix dimensions
   * @param matrix Block matrix to store the data
   * @return true if successful, false on error
   */
  bool ReadDataMatrix(
          std::string sFileName,
          std::vector<int> &matrix_info,
          bse::block_matrix<std::complex<double> > &matrix
  ) {
    int itmp;
    std::ifstream fin(sFileName);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> bm_insert(matrix_info[0], matrix_info[0]);

    if (fin.fail()) {
      std::cerr << "Error: Read file:" << sFileName << std::endl;
      // return false;
      exit(-1);
    }
    //ignore oneline
    std::string str;
    int iMatrixDim;
    fin >> itmp >> iMatrixDim;
    getline(fin, str);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> bm(iMatrixDim, iMatrixDim);

    int idx = 0;
    int iblock1, iblock2;

    // Read data line by line
    while (getline(fin, str)) {
      std::string token;
      std::istringstream stream(str);
      if (idx % (iMatrixDim * iMatrixDim + 1) == 0) {
        // Read block indices
        int ivalue[2];
        int tmpidx = 0;
        while (getline(stream, token, ',')) {
          ivalue[tmpidx] = stoi(token);
          tmpidx++;
        }
        iblock1 = ivalue[0];
        iblock2 = ivalue[1];
        idx++;
      } else {
        // Read matrix elements
        double cvalue[2];
        int tmpidx = 0;
        while (getline(stream, token, ',')) {
          cvalue[tmpidx] = stod(token);
          tmpidx++;
        }
        bm((idx - 1) / iMatrixDim, (idx - 1) % iMatrixDim) = std::complex<double>(cvalue[0], cvalue[1]);
        idx++;
      }

      // When block is complete, assign to matrix
      if (idx % (iMatrixDim * iMatrixDim + 1) == 0) {
        if (iMatrixDim > matrix_info[0]) {
          // Extract center block if input is larger
          int iLeftCorner = (iMatrixDim - matrix_info[0]) / 2;
          bm_insert = bm.block(iLeftCorner, iLeftCorner, matrix_info[0], matrix_info[0]);
        } else {
          if (iMatrixDim != matrix_info[0]) {
            return false;
          }
          bm_insert = bm;
        }
        matrix.assign(iblock1, iblock2, bm_insert);
        idx = 0;
      }
    }
    fin.close();
    return true;
  }

  /**
   * @brief Read general matrix data from file
   * @param sFileName Input filename
   * @param matrix_info Vector containing matrix dimensions
   * @param matrix Block matrix to store the data
   * @return true if successful, false on error
   */
  bool ReadData(
          std::string sFileName,
          std::vector<int> &matrix_info,
          bse::block_matrix<std::complex<double> > &matrix
  ) {
    // int itmp;
    std::ifstream fin(sFileName);
    if (fin.fail()) {
      std::cerr << "Error: Read file:" << sFileName << std::endl;
      // return false;
      exit(-1);
    }
    std::string str;
    // int iMatrixDim;
    int nMatsubara;
    fin >> nMatsubara;
    // std::cout << nMatsubara << std::endl;

    // ignore one line
    getline(fin, str);

    // Read data block by block
    while (getline(fin, str)) {
      std::string token;
      std::istringstream stream(str);

      // Get block position
      int ivalue[2];
      {
        // check if the line expresses block position
        std::string key;
        stream >> key;
        if (key != "block"){
          std::cerr << "Error: invalid file structure " << sFileName << std::endl;
          // return false;
          exit(-1);
        }
        int tmpidx = 0;
        while (getline(stream, token, ',')) {
          ivalue[tmpidx] = stoi(token);
          tmpidx++;
        }
      }
      int iblock1 = ivalue[0];
      int iblock2 = ivalue[1];

      int iDim1 = matrix_info[iblock1];
      int iDim2 = matrix_info[iblock2];

      // Initialize block matrix
      Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> bm(iDim1, iDim2);

      // Read data elements
      for(int i=0; i<iDim1; i++){
        for(int j=0; j<iDim2; j++){
          getline(fin, str);
          std::string token;
          std::istringstream stream(str);
          double cvalue[2];
          int tmpidx = 0;
          while (getline(stream, token, ',')) {
            cvalue[tmpidx] = stod(token);
            tmpidx++;
          }
          bm(i, j) = std::complex<double>(cvalue[0], cvalue[1]);
          // idx++;
        }
      }

      // Assign data to block_matrix
      {
        matrix.assign(iblock1, iblock2, bm);
        // idx = 0;
      }
    }

    fin.close();
    return true;
  }


  std::vector<int>
  GetMatrixInfo(std::string sFileName_chi0_loc, std::string sFileName_chi_loc, std::string sFileName_ch0_q);

  std::vector<int>
  GetMatrixInfo(std::string sFileName_chi0_loc, std::string sFileName_ch0_q);

  bool GetMatrixInfoCore(std::string sFileName, std::vector<int> &matrix_info);

  int GetMatsubaraSize(std::vector<std::string> const& sFileNames);

  int GetMatsubaraSizeCore(std::string sFileName);

  std::vector<int> GetMatrixInfo(std::string sFileName);

private:
  std::string FileNameChi0_loc;  ///< Filename for local bare susceptibility
  std::string FileNameChi_loc;   ///< Filename for local full susceptibility
  std::string FileNameChi0_q;    ///< Filename for q-dependent bare susceptibility
};

