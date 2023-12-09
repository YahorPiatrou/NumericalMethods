#pragma once

#include <iostream>
#include <vector>
#include "GaussianMethod.h"

using namespace std;

vector<double> DifferentialOfFuncUOnT(vector<double> , vector<vector<double>> , vector<double> );
vector<double> ComposeVectorF(vector<double> , vector<double> , double, vector<vector<double>> , vector<double> );
vector<vector<double>> CalculateJacobianMatrix(vector<double> , vector<double> , double , double , double , vector<vector<double>> , vector<double> );
vector<double> SolveNewtonsMethod(vector<double> , vector<double> , double , vector<vector<double>> , vector<double> , double ,  const double );

