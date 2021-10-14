//
// Created by lois on 07/03/2021.
//

#include <vector>
#include <cassert>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include <Eigen/Dense>
#include <chrono>
#include <functional>
#include <fstream>

#include <CLI11.hpp>
#include "Points/VecX.h"
#include "../Integration/Integration.h"
#include "Stats/Stats.h"

int rank( const std::vector< VecX<4> >& points, std::function<double(const VecX<4>&)> fun){

    int n = points.size();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix(n, n);

    for(int i = 0; i < n; i++) {
        VecX<4> point;
        point[0]= points[i][0];
        point[1]= points[i][1];
        for(int j = 0; j < n; j++) {
            point[2]= points[j][2];
            point[3]= points[j][3];
            // eval function
            matrix(i, j) = fun(point);
        }
    }

    // compute rank
    Eigen::FullPivLU< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > LU(matrix);
    LU.setThreshold(1e-5);

    return LU.rank();

}

std::function<double(const VecX<4>&)> getGauss(int i){

    std::function<double(const VecX<4>&)> f = [=](const VecX<4>& p)->double {
        return getMultivariateGaussian(4, p.data(), tab_Gauss4D[i].mu, tab_Gauss4D[i].mxCInv)
            + getMultivariateGaussian(4, p.data(), tab_Gauss4D[i+1].mu, tab_Gauss4D[i+1].mxCInv);
    };
    return f;
}

double heavyside(const VecX<4>& p, double* mu, double* normal){
    double dotProduct = 0;
    for (int idim = 0; idim < p.dim() ; idim++) {
        dotProduct += (p[idim] - mu[idim]) * normal[idim];
    }
    return  (dotProduct > 0) ? 1. : 0.;
}

std::function<double(const VecX<4>&)> getHeavy(int i){
    std::function<double(const VecX<4>&)> f = [=](const VecX<4>& p)->double {
        double* mu = tab_Heaviside4D[i].muDiscotinuity;
        double* normal = tab_Heaviside4D[i].normal;

        return heavyside(p, mu, normal);
    };
    return f;
}

int main(int argc, const char** argv){


    CLI::App app { "Generating Integrands Files filtered by rank"};

    int dim = 4;
    int nbGauss = 1;
    bool gaussian = false;
    app.add_option("--gauss", gaussian, "Filters Gaussians");
    bool heavyside = false;
    app.add_option("--heavy", heavyside, "Filters Heavysides");
    int targetRank = 10;
    app.add_option("-r", targetRank, "Target Rank (default 10)");
    int error = 300;
    app.add_option("-e", error, "Authorized distance to target Rank (default 0)");
    std::string output = "out.h";
    app.add_option("-o", output, "Output file name");

    CLI11_PARSE(app, argc, argv)

    //init output file
    std::ofstream out(output);
    out << std::setprecision(16);

    //Init points used to compute rank
    std::mt19937_64 gen(42);
    std::uniform_real_distribution<double> unif(0,1);
    std::vector<VecX<4>> pts(256);
    for (VecX<4>& p : pts){
        for (int i = 0; i < 4; ++i){
            p[i] = unif(gen);
        }
    }

    //init timer
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();

    //init stat computation
    Stats stats;
    stats.analytical(0.);

    //selected integrands indices
    std::vector<int> selected;


    for (int i = 0; i < nbGauss; ++i){
        if (i%5 == 0) {
            std::cout << std::setfill('0') << std::setw(5) << i << '\r';
            std::cout.flush();
        }
        //compute rank
        int v;
        if (gaussian){
            v = rank(pts, getGauss(i));
        }else if (heavyside){
            v = rank(pts, getHeavy(i));
        }
        //if rank in selection window then add to selected
        if (abs(targetRank - v) <= error){
            stats.addData(1, v);
            selected.push_back(i);
        }
    }
    std::cout << std::endl;

    //output time and stats
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    std::cerr << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
    std::cerr << stats << std::endl;


    if (gaussian){
        //output to file
        out << "#include \"../Integration/Integration.h\"\n"
               "\n"
               "t_Gauss4D tab_Gauss4D[" << selected.size() << "] = \n"
                                                              "{";

        for (int indSelected = 0; indSelected < selected.size(); ++indSelected){
            int v = selected[indSelected];
            if (indSelected != 0){
                out << ", " << std::endl;
            }
            out << "{" << tab_Gauss4D[v].integral << ", ";
            for (int i = 0; i < dim; ++i){
                out << tab_Gauss4D[v].mu[i] << ", ";
            }
            for (int i = 0; i < dim*dim; ++i){
                if (i != 0){
                    out << ", ";
                }
                out << tab_Gauss4D[v].mxCInv[i];
            }
            out << "}";
        }

        out << std::endl << "};" << std::endl;
    } else if (heavyside) {

        //output to file
        out << "#include \"../Integration/Integration.h\"\n"
               "\n"
               "t_Heaviside4D tab_Heaviside4D[" << selected.size() << "] = \n"
                                                              "{";

        for (int indSelected = 0; indSelected < selected.size(); ++indSelected){
            int v = selected[indSelected];
            if (indSelected != 0){
                out << ", " << std::endl;
            }
            out << "{" << tab_Heaviside4D[v].integral << ", ";
            for (int i = 0; i < dim; ++i){
                out << tab_Heaviside4D[v].muDiscotinuity[i] << ", ";
            }
            for (int i = 0; i < dim; ++i){
                if (i != 0){
                    out << ", ";
                }
                out << tab_Heaviside4D[v].normal[i];
            }
            out << "}";
        }

        out << std::endl << "};" << std::endl;
    }

    out.close();
}
