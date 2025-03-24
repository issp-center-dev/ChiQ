#pragma once

#include "./block_matrix.hpp"
#include <iostream>
#include <string>


namespace bse {

/**
 * @brief Class for tracking memory usage of block matrices
 * 
 * This class keeps track of total memory usage by block matrices
 * and can optionally print verbose memory information.
 */
//   template<typename Scalar>
    class memory_info {
    public:
        /**
         * @brief Constructor
         * @param verbose If true, prints memory info for each added matrix
         */
        memory_info(bool verbose=false)
            : mem_total(0), verbose(verbose) {
        }

        /**
         * @brief Add memory usage of a block matrix
         * @tparam Scalar Type of matrix elements
         * @param other Block matrix to track memory for
         * @param str Optional description string for verbose output
         */
        template<typename Scalar>
        void add(const block_matrix<Scalar> &other, const std::string &str=""){
            size_t bytes = other.get_bytes();
            mem_total += bytes;
            if (verbose){
                std::cout << "meminfo: " << bytes << " Bytes (" << str << ")" <<std::endl;
            }
        }

        /**
         * @brief Get total memory usage
         * @return Total memory usage in bytes
         */
        size_t get_total(){
            return mem_total;
        }

        /**
         * @brief Print summary of total memory usage
         */
        void summary(){
            std::cout << "meminfo: " << mem_total << " Bytes (total)" <<std::endl;
        }

    private:
        size_t mem_total;  ///< Total memory usage in bytes
        bool verbose;      ///< Whether to print verbose memory info
    };
}
