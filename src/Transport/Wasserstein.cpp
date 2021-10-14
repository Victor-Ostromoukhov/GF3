#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <limits>
#include <ctime>

#include <geogram/delaunay/delaunay_2d.h>
#include <geogram/basic/common.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/numerics/predicates.h>
#include <exploragram/optimal_transport/optimal_transport_2d.h>
#include <exploragram/optimal_transport/optimal_transport.h>
#ifndef _MSC_VER
#include <sys/time.h>
#endif
#include "../Points/VecX.h"

void usage(const char** argv){
    std::cerr << argv[0] << " -i <inputFileName> [-o <OutputFileName>]" << std::endl;
}

void init_zone_mesh(GEO::Mesh& m){
    GEO::vector<double> points(8);
    points[0]=0;
    points[1]=0;
    points[2]=0;
    points[3]=1;
    points[4]=1;
    points[5]=1;
    points[6]=1;
    points[7]=0;

    m.vertices.assign_points(points, 2, true);
    m.edges.create_edge(0, 1);
    m.edges.create_edge(1, 2);
    m.edges.create_edge(2, 3);
    m.edges.create_edge(3, 0);

    GEO::vector<GEO::index_t> facets(4);
    facets[0]=0;
    facets[1]=1;
    facets[2]=2;
    facets[3]=3;
    m.facets.create_polygon(facets);

}
void init_zone_mesh3(GEO::Mesh& m){
    GEO::vector<double> points(8);
    points[0]=0;
    points[1]=0;
    points[2]=0;
    points[3]=3;
    points[4]=3;
    points[5]=3;
    points[6]=3;
    points[7]=0;

    m.vertices.assign_points(points, 2, true);
    m.edges.create_edge(0, 1);
    m.edges.create_edge(1, 2);
    m.edges.create_edge(2, 3);
    m.edges.create_edge(3, 0);

    GEO::vector<GEO::index_t> facets(4);
    facets[0]=0;
    facets[1]=1;
    facets[2]=2;
    facets[3]=3;
    m.facets.create_polygon(facets);

}

double Wasserstein1(const GEO::vector<double>& points, const GEO::vector<double>& centroids){
    int DIM = 2;
    double d = 0;
    for (unsigned int i = 0; i < points.size(); i += DIM){
        double x = points[i] - centroids[i];
        double y = points[i+1] - centroids[i+1];
        d += std::sqrt(x*x + y*y);
    }
    return d;
}

bool read_points_from_file(std::istream& in, GEO::vector<double>& points){
    points.clear();
    std::string line;
    while(std::getline(in, line)){
        int c = line.find_first_not_of(" \t");
        if (line[c] != '#'){
            std::istringstream lineIn(line);
            double d;
            while (lineIn >> d){
                points.push_back(d);
            }
        } else {
            return true;
        }
    }
    return points.size() != 0;
}

//Compute wasserstein2 distance between p and the right triangle pAB with the right angle in B
double rightTriangleW2(const GEO::vec2& p, const GEO::vec2& A, const GEO::vec2& B){
/*
    double pB2 = p.distance2(B); //a2
    double pB = std::sqrt(pB2); //a
    double AB2 = A.distance2(B); //b2
    double AB = std::sqrt(AB2); //b

    double root = std::sqrt(pB2 + AB2);

    double res = AB * root * (5 * pB2 + 2 * AB2) + 3 * pB2*pB2 * (std::log(AB + root) - std::log(pB));

    res /= 24;

    return res;
*/
    double pB = p.distance(B); //h
    double AB = A.distance(B); //a

    double res = 3 * pB*pB*pB * AB + pB * AB*AB*AB;

    res /= 12;

    return res;

}

//Compute wasserstein2 distance between p and a polygon with vertices given in counter-clockwise order in points
double polygonWasserstein2(const GEO::vec2& p, const std::vector<GEO::vec2>& points){
    double res = 0.;
    int nbPoints = points.size();

    for (int i = 0; i < nbPoints; ++i){

        GEO::vec2 A = points[i];
        GEO::vec2 B = points[(i+1) % nbPoints];
        GEO::vec2 ab = GEO::normalize(B - A);

        double v = GEO::dot(ab, p - A);
        GEO::vec2 projP = A + v * ab;

        double resTriangle;
        double sign = GEO::PCK::orient_2d(GEO::vec2(p), GEO::vec2(A), GEO::vec2(projP));
        resTriangle = sign * rightTriangleW2(p, A, projP);
        sign = GEO::PCK::orient_2d(p, projP, B);
        resTriangle += sign * rightTriangleW2(p, B, projP);

        res += resTriangle;

    }

    return res;
}

double Wasserstein2PostLaguerre(const GEO::vector<double>& points, const GEO::Mesh& RVD){

    double res = 0;
    int nbskipped = 0;

    int* id = (int*)RVD.facets.attributes().find_attribute_store("chart")->data();
    for (unsigned int i = 0; i < RVD.facets.nb(); ++i){
        bool skip = false;
        if (id[i] < int(points.size()) / 2) {
            std::vector<GEO::vec2> vertices;
            for (unsigned int j = 0; j < RVD.facets.nb_vertices(i) && !skip; ++j) {
                vertices.emplace_back(RVD.vertices.point(RVD.facets.corner(i, j)).x,
                        RVD.vertices.point(RVD.facets.corner(i, j)).y);
                /*if (vertices.back().x == 0. || vertices.back().x == 1. || vertices.back().y == 0. || vertices.back().y == 1.){
                    skip = true;
                }*/
            }
            res += skip ? 0. : polygonWasserstein2(GEO::vec2(points[2 * id[i]], points[2 * id[i]+1]), vertices);
            nbskipped += skip ? 1 : 0;
        }
    }

    //std::cout << "Skipped: " << nbskipped << std::endl;

    return -res;
}

double Wasserstein2(const std::vector<VecX<2> >& inputPoints){

    GEO::Mesh m;
    init_zone_mesh(m);

    GEO::index_t nbPoints = GEO::index_t(inputPoints.size());
    GEO::vector<double> points(nbPoints*2);

    for (int i = 0; i < nbPoints; ++i){
        points[2*i] = inputPoints[i][0];
        points[2*i+1] = inputPoints[i][1];
    }

    GEO::Mesh RVD;
    GEO::vector<double> centroids(points.size());
    GEO::compute_Laguerre_centroids_2d(&m, nbPoints, points.data(), centroids.data(), &RVD, false);

    return Wasserstein2PostLaguerre(points, RVD);

}