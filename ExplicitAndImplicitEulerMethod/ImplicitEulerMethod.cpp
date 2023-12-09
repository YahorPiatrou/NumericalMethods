#include "ImplicitEulerMethod.h"

vector<vector<double>> FillMatrix(vector<double> lambdaValues) 
{
	vector<vector<double>> A(3, vector<double>(3));
	A[0] = { 2 * lambdaValues[0] + 4 * lambdaValues[1] , 2 * (lambdaValues[0] - lambdaValues[1]) , 2 * (lambdaValues[0] - lambdaValues[1]) };
	A[1] = { 2 * (lambdaValues[0] - lambdaValues[1]) , 2 * lambdaValues[0] + lambdaValues[1] + 3 * lambdaValues[2] , 2 * lambdaValues[0] + lambdaValues[1] - 3 * lambdaValues[2] };
	A[2] = { 2 * (lambdaValues[0] - lambdaValues[1]) , 2 * lambdaValues[0] + lambdaValues[1] - 3 * lambdaValues[2] , 2 * lambdaValues[0] + lambdaValues[1] + 3 * lambdaValues[2] };
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[i].size(); j++)
			A[i][j] = A[i][j] * 1 / 6;
	}
	return A;
}

vector<double> FillVector(vector<double> lambdaValues) 
{
	vector<double> b = { 4 * lambdaValues[0] + 2 * lambdaValues[1] , 4 * lambdaValues[0] - lambdaValues[1] - 9 * lambdaValues[2] , 4 * lambdaValues[0] - lambdaValues[1] + 9 * lambdaValues[2] };
	for (int i = 0; i < b.size(); i++)
		b[i] = b[i] *  - 1 / 6;
	return b;
}

vector<double> CalculateLocalError(vector<double> previousVectorYk, vector<double> vectorYk, vector<double> nextVectorYk, double integrationStep, double previousIntegrationStep) 
{
	vector<double> vectorOfLocalError(3);
	for (int i = 0; i < vectorOfLocalError.size(); i++)
	{
		vectorOfLocalError[i] = nextVectorYk[i] - vectorYk[i] - (vectorYk[i] - previousVectorYk[i]) * (integrationStep / previousIntegrationStep);
		vectorOfLocalError[i] = vectorOfLocalError[i] * (-integrationStep / (integrationStep + previousIntegrationStep));
	}
		return vectorOfLocalError;
}

bool CheckValidityForError(vector<double> vectorOfLocalError, double acceptableLocalError) 
{
	for (int i = 0; i < vectorOfLocalError.size(); i++) 
	{
		if (abs(vectorOfLocalError[i]) < acceptableLocalError) 
			return false;
	}
	return true;
}

double CalculateNextStepOfIntegration(vector<double> vectorOfLocalError, double acceptableLocalError, double integrationStep) 
{
	vector<double> vectorOfIntegrationSteps(vectorOfLocalError.size());
	double minValue;
	for (int i = 0; i < vectorOfLocalError.size(); i++)
		vectorOfIntegrationSteps[i] = sqrt(acceptableLocalError / abs(vectorOfLocalError[i])) * integrationStep;
	return min(vectorOfIntegrationSteps);
}

vector<double> SolveImplicitEulerMethod(vector<vector<double>> A,vector<double> b,vector<double> valueUAtZero, vector<double> parameterBoundariesT,double acceptableLocalError)
{
	double minIntegrationStep = abs((parameterBoundariesT[1] - parameterBoundariesT[0]) * acceptableLocalError);
	double maxIntegrationStep = abs((parameterBoundariesT[1] - parameterBoundariesT[0]) * pow(10, -2));

	double integrationStep = minIntegrationStep;
	double previousIntegrationStep = minIntegrationStep;
	double nextIntegrationStep = 0;
	double parameterT = parameterBoundariesT[0];
	double nextParameterT = 0;

	vector<double> vectorYk = valueUAtZero;
	vector<double> previousVectorYk = valueUAtZero;
	vector<double> nextVectorYk = valueUAtZero;
	vector<double> vectorOfLocalError(3); // those errors that we calculate using formula 3.16

	while (parameterT < parameterBoundariesT[1])
	{
		nextParameterT = parameterT + integrationStep;
		nextVectorYk = SolveNewtonsMethod(nextVectorYk, vectorYk, integrationStep, A, b, acceptableLocalError, 0.05);
		vectorOfLocalError = CalculateLocalError(previousVectorYk, vectorYk, nextVectorYk, integrationStep, previousIntegrationStep);

		if (CheckValidityForError(vectorOfLocalError, acceptableLocalError))
		{
			integrationStep = integrationStep / 2;
			nextParameterT = parameterT;
			nextVectorYk = vectorYk;
		}
		else 
		{
			nextIntegrationStep = CalculateNextStepOfIntegration(vectorOfLocalError, acceptableLocalError, integrationStep);
			if (nextIntegrationStep > maxIntegrationStep) 
				nextIntegrationStep = maxIntegrationStep;
			previousVectorYk = vectorYk;
			vectorYk = nextVectorYk;
			previousIntegrationStep = integrationStep;
			integrationStep = nextIntegrationStep;
			parameterT = nextParameterT;
		}
	}
	return vectorYk;
}