#include "GaussianMethod.h"

vector<double> SolveGaussMethod(vector<vector<double>> A, vector<double> b)
{
	int vectorSize = A.size();
	for (int i = 0; i < vectorSize; i++)
	{
		int k = i;
		for (int j = i + 1; j < vectorSize; j++)
		{
			if (abs(A[j][i]) > abs(A[k][i]))
				k = j;
		}
		swap(A[i], A[k]);
		swap(b[i], b[k]);
		double div = A[i][i];
		for (int j = i; j < vectorSize; j++)
			A[i][j] /= div;
		b[i] /= div;
		for (int j = i + 1; j < vectorSize; j++)
		{
			double mult = A[j][i];
			for (int k = i; k < vectorSize; k++)
				A[j][k] -= mult * A[i][k];
			b[j] -= mult * b[i];
		}
	}
	vector<double> vectorOfRoots(vectorSize);
	for (int i = vectorSize - 1; i >= 0; i--)
	{
		vectorOfRoots[i] = b[i];
		for (int j = i + 1; j < vectorSize; j++)
			vectorOfRoots[i] -= A[i][j] * vectorOfRoots[j];
	}
	return vectorOfRoots;
}

vector<double> abs(vector<double> vectorOfElement)
{
	vector<double> result(vectorOfElement.size());
	for (int i = 0; i < vectorOfElement.size(); i++)
		result[i] = abs(vectorOfElement[i]);
	return result;
}
vector<double> sqrt(vector<double> vectorOfElement)
{
	vector<double> result(vectorOfElement.size());
	for (int i = 0; i < vectorOfElement.size(); i++)
		result[i] = sqrt(vectorOfElement[i]);
	return result;
}

double max(vector<double> vectorOfElement)
{
	double maxValue = vectorOfElement[0];
	for (int i = 1; i < vectorOfElement.size(); i++)
	{
		if (vectorOfElement[i] > maxValue)
			maxValue = vectorOfElement[i];
	}
	return maxValue;
}

double min(vector<double> vectorOfElement) 
{
	double minValue = vectorOfElement[0];
	for (int i = 1; i < vectorOfElement.size(); i++) 
	{
		if (vectorOfElement[i] < minValue) 
			minValue = vectorOfElement[i];
	}
	return minValue;
}