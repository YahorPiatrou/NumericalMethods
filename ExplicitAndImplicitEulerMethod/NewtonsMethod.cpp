#include "NewtonsMethod.h"

vector<double> DifferentialOfFuncUOnT(vector<double> vectorOfUValues, vector<vector<double>> A, vector<double> b) 
{
    vector<double> vectorOfDiff(A.size());
    for (int i = 0; i < vectorOfDiff.size(); i++)
    {
        vectorOfDiff[i] = 0;
        for (int j = 0; j < vectorOfDiff.size(); j++)
            vectorOfDiff[i] += A[i][j] * vectorOfUValues[j];
    }
    for (int i = 0; i < b.size(); i++)
        vectorOfDiff[i] = vectorOfDiff[i] - b[i];
    return vectorOfDiff;
}

// 3.14
vector<double> ComposeVectorF(vector<double> nextVectorYk, vector<double> vectorYk, double integrationStep, vector<vector<double>> A, vector<double> b) 
{
    vector<double> vectorOfDiff = DifferentialOfFuncUOnT(nextVectorYk, A, b);
    for (int i = 0; i < vectorOfDiff.size(); i++)
        vectorOfDiff[i] = nextVectorYk[i] - vectorYk[i] - vectorOfDiff [i] * integrationStep;
    return vectorOfDiff;
}

vector<vector<double>> CalculateJacobianMatrix(vector<double> nextVectorYk, vector<double> vectorOfYk, double integrationStep, double relIncrement, double localError, vector<vector<double>> A, vector<double> b) 
{
    vector<double> nextVectorYk_ = nextVectorYk;
    vector<double> firstVectorOfDiff(3);
    vector<double> secondVectorOfDiff(3);
    vector<vector<double>> jacobianMatrix(3, vector<double>(3));

    double differentialOfY = 0;

    for (int i = 0; i < jacobianMatrix.size(); i++) 
    {
        differentialOfY = relIncrement * nextVectorYk_[i];
        if (abs(differentialOfY) < localError)
            differentialOfY = localError;
        nextVectorYk_[i] += differentialOfY;
        firstVectorOfDiff = ComposeVectorF(nextVectorYk_, vectorOfYk, integrationStep, A, b);
        secondVectorOfDiff = ComposeVectorF(nextVectorYk, vectorOfYk, integrationStep, A, b);
        for (int j = 0; j < jacobianMatrix.size(); j++)
            jacobianMatrix[j][i] = (firstVectorOfDiff[j] - secondVectorOfDiff[j]) / differentialOfY;
        nextVectorYk_[i] -= differentialOfY;
    }
    return jacobianMatrix;
}

vector<double> SolveNewtonsMethod(vector<double> nextVectorYk, vector<double> vectorYk, double integrationStep, vector<vector<double>> A, vector<double> b, double localError, double relIncrement) 
{
    double firstDelta = max(abs(ComposeVectorF(nextVectorYk, vectorYk, integrationStep, A, b)));
    double secondDelta = 1;

    vector<double> residualVector(nextVectorYk.size());
    vector<double> vectorOfGaussSolution(nextVectorYk.size());
    vector<double> auxiliaryVector(nextVectorYk.size());
    vector<vector<double>> jacobianMatrix(nextVectorYk.size(), vector<double>(nextVectorYk.size()));

    while (firstDelta > localError || secondDelta > localError) 
    {
        residualVector = ComposeVectorF(nextVectorYk, vectorYk, integrationStep, A, b);
        for (int i = 0; i < residualVector.size(); i++)
            residualVector[i] = -residualVector[i];
        jacobianMatrix = CalculateJacobianMatrix(nextVectorYk, vectorYk, integrationStep, relIncrement, localError, A, b);
        vectorOfGaussSolution = SolveGaussMethod(jacobianMatrix, residualVector);
        auxiliaryVector = nextVectorYk;
        for (int i = 0; i < nextVectorYk.size(); i++)
            nextVectorYk[i] = nextVectorYk[i] + vectorOfGaussSolution[i];
        firstDelta = max(abs(residualVector));
        for (int i = 0; i < auxiliaryVector.size(); i++) 
            auxiliaryVector[i] = abs(auxiliaryVector[i]) < 1 ? abs(vectorOfGaussSolution[i]) : abs(vectorOfGaussSolution[i] / auxiliaryVector[i]);
        secondDelta = max(auxiliaryVector);
    }
    return nextVectorYk;
}
