//
// Created by lpaulin on 10/04/19.
//

#ifndef STATTOOLSSAMPLING_STATS_H
#define STATTOOLSSAMPLING_STATS_H

#include <vector>
#include <iostream>
#include <limits>
#include <iomanip>

class Stats {
private:

    int nbValues = 0;
    double nbPoints = 0;
    double mean = 0.;
    double mse = 0.;
    double var = 0.;
    double mini = std::numeric_limits<double>::max();
    double maxi = -std::numeric_limits<double>::max();
    double analyticalValue = 0.;
public:
    Stats(){}

    /** Call this function once after adding all Data. Tu compute MSE and var value **/
    void compute(){

        if(nbValues <= 1)
            var = 0;
        else
            var /= nbValues-1;
        mse /= nbValues;

    }

    inline void addData(int nbP, double value){
        nbValues += 1;

        double dnbp = double(nbP - nbPoints) / nbValues;
        nbPoints += dnbp;

        double dvalue = (value - mean) / nbValues;
        mean += dvalue;
        var += dvalue*(value - mean);
        mse += (value - analyticalValue) * (value - analyticalValue);

        if (value > maxi)
            maxi = value;
        if (value < mini)
            mini = value;

    }

    void analytical(double a){
        analyticalValue = a;
    }

    friend std::ostream& operator<<(std::ostream& out, const Stats& e);
};

std::ostream& operator<<(std::ostream& out, const Stats& stats){
    int prec = out.precision();
    out << "#Nbpts\t#Mean\t\t#Var\t#Min\t#Max\t#Analytical\t#MSE\t#NbPtsets" << std::endl;
    out << stats.nbPoints << "\t";
    out << std::setprecision(15) << stats.mean << "\t";
    out << stats.var << "\t";
    out << stats.mini << "\t";
    out << stats.maxi << "\t";
    out << stats.analyticalValue << "\t";
    out << stats.mse << "\t";
    out << stats.nbValues;
    out << std::setprecision(prec);
    return out;
}


#endif //STATTOOLSSAMPLING_STATS_H
