#include "ExplicitEulerMethod.h"

vector<double> SolveExplicitEulerMethod(vector<double> valueUAtZero, vector<double> parameterBoundariesT, double parameterOmega, double acceptableLocalError)
{
	vector<double> vectorOfRoots = valueUAtZero; // vector roots
	vector<double> vectorOfDifferentialEquations(valueUAtZero.size());
	double parameterT = parameterBoundariesT[0];
	double maxIntegrationStep = (parameterBoundariesT[1] - parameterBoundariesT[0]) * acceptableLocalError;
	double integrationStep = 0;
	while (parameterT < parameterBoundariesT[1])
	{
		vectorOfDifferentialEquations = CalculateDifferentialEquations(vectorOfRoots, parameterT, parameterOmega);
		integrationStep = CalculateIntegrationStep(vectorOfDifferentialEquations, maxIntegrationStep, acceptableLocalError);
		vectorOfRoots = CalculateSolutionOfEquations(vectorOfRoots, vectorOfDifferentialEquations, integrationStep);
		parameterT += integrationStep;
	}
	return vectorOfRoots;
}

double u0dt(vector<double> vectorOfValuesU, double parameterT)
{
	if (parameterT == 0)
		return 0;
	else
		return -1 * vectorOfValuesU[0] * vectorOfValuesU[1] + sin(parameterT) / parameterT;
}

double u1dt(vector<double> vectorOfValuesU, double parameterT, double parameterOmega)
{
	return -1 * pow(vectorOfValuesU[1], 2) + (((2.5 + parameterOmega / 40) * parameterT) / (1 + pow(parameterT, 2)));
}

vector<double> CalculateDifferentialEquations(vector<double> vectorOfValuesU, double parameterT, double parameterOmega)
{
	vector<double> vectorOfDifferentialEquations(vectorOfValuesU.size());
	vectorOfDifferentialEquations[0] = u0dt(vectorOfValuesU, parameterT);
	vectorOfDifferentialEquations[1] = u1dt(vectorOfValuesU, parameterT, parameterOmega);
	return vectorOfDifferentialEquations;
}

double CalculateIntegrationStep(vector<double> vectorOfDifferentialEquations, double maxIntegrationStep, double acceptableLocalError)
{
	vector<double> integrationStep(vectorOfDifferentialEquations.size());
	for (int i = 0; i < integrationStep.size(); i++)
		integrationStep[i] = acceptableLocalError / (abs(vectorOfDifferentialEquations[i]) + (acceptableLocalError / maxIntegrationStep));
	double minElement = integrationStep[0];
	for (int i = 1; i < integrationStep.size(); i++)
	{
		if (integrationStep[i] < minElement)
			minElement = integrationStep[i];
	}
	return minElement;
}

vector<double> CalculateSolutionOfEquations(vector<double> vectorOfValuesU, vector<double> vectorOfDifferentialEquations, double integrationStep)
{
	vector<double> vectorOfSolution(vectorOfDifferentialEquations.size());
	for (int i = 0; i < vectorOfDifferentialEquations.size(); i++)
		vectorOfSolution[i] = vectorOfValuesU[i] + vectorOfDifferentialEquations[i] * integrationStep;
	return vectorOfSolution;
}