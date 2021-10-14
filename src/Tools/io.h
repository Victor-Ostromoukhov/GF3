//
// Created by lpaulin on 21/09/20.
//

#ifndef SOBOL_IO_H
#define SOBOL_IO_H

#include <string>
#include <iostream>
#include <fstream>

#include "../Samplers/SobolGenerator1D.h"

#define N_DIMS 2
struct t_pointND {double coord[N_DIMS] ; };

template <typename integer>
inline void save_dv_2d(const std::string& filename,
                       const SobolGenerator1D<integer>& dim1,
                       const SobolGenerator1D<integer>& dim2) {

    std::cout << "Writing 2d mk data into " << filename << std::endl;

    std::ofstream file(filename, std::ifstream::trunc);
    file << dim1.d << "\t" << dim1.s << "\t" << dim1.a << "\t";
    for (int k = 0; k < dim1.s; ++k) {
        file << dim1.m[k] << " ";
    }
    file << std::endl;
    file << dim2.d << "\t" << dim2.s << "\t" << dim2.a << "\t";
    for (int k = 0; k < dim2.s; ++k) {
        file << dim2.m[k] << " ";
    }
    file << std::endl;
    file.close();
}	// save_dv_2d

inline void save_mse_2d(const std::string& filename, const std::vector<t_pointND>& mse_data) {
    std::cout << "Writing 2d mse data into " << filename << std::endl;
    std::ofstream file(filename, std::ifstream::trunc);
    for (int i = 0; i < mse_data.size(); i++) {
        file << mse_data[i].coord[0] << " " << mse_data[i].coord[1] << std::endl;
    }
    file.close();
}	// save_mse_2d

inline void save_pts_2d(const std::string& filename, const std::vector<t_pointND>& pts) {
    std::cout << "Writing 2d pts into " << filename << std::endl;
    std::ofstream file(filename, std::ifstream::trunc);
    for (int i = 0; i < pts.size(); i++) {
        file << pts[i].coord[0] << " " << pts[i].coord[1] << std::endl;
    }
    file.close();
}	// save_pts_2d


#endif //SOBOL_IO_H
