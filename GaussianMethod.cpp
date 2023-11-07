#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

bool operator / (vector<double> first_vector,vector<double> second_vector);
vector<double> operator * (vector<vector<double>> matrix, vector <double> vector);
bool CheckProportionality(vector<vector<double>> A, vector<double> b);
vector<double> GaussSolution(vector<vector<double>> A, vector<double> B);
double FindingResidualVector(vector<vector<double>> A, vector<double> b, vector<double> answer);
double ErrorCalculation(vector<double> x1, vector<double> x2);
void OutputMatrix(vector<vector<double>> A, vector<double> B);

bool operator / (vector<double> first_vector, vector<double> second_vector)
{
	int size = first_vector.size();
	double proportionalityNumber = first_vector[0] / second_vector[0];
	for (int i = 1; i < size; i++)
		if (proportionalityNumber == first_vector[i] / second_vector[i])
			continue;
		else
			return false;
	return true;
}

vector<double> operator * (vector<vector<double>> matrix, vector <double> vector_)
{
	vector<double> result(matrix.size());
	for (int i = 0; i < result.size(); i++) {
		result[i] = 0;
		for (int j = 0; j < result.size(); j++) {
			result[i] += matrix[i][j] * vector_[j];
		}
	}
	return result;
}

bool CheckProportionality(vector<vector<double>> A, vector<double> b)
{
	int n = A.size();
	vector<vector<double>> A_b(n, vector<double>(n + 1, 0));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n + 1; j++) {
			if (j == n)
				A_b[i][j] = b[i];
			else
				A_b[i][j] = A[i][j];
		}
	}
	OutputMatrix(A, b);
	bool proportional = false;
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			proportional = A_b[i] / A_b[j];
			if (proportional)
				return true;
		}
	}
	return false;
}

vector<double> GaussSolution(vector<vector<double>> A, vector<double> B)
{
	int n = A.size();
	for (int i = 0; i < n; i++)
	{
		int k = i;
		for (int j = i + 1; j < n; j++)
		{
			if (abs(A[j][i]) > abs(A[k][i]))
				k = j;
		}
		swap(A[i], A[k]);
		swap(B[i], B[k]);
		double div = A[i][i];
		for (int j = i; j < n; j++)
			A[i][j] /= div;
		B[i] /= div;
		for (int j = i + 1; j < n; j++)
		{
			double mult = A[j][i];
			for (int k = i; k < n; k++)
				A[j][k] -= mult * A[i][k];
			B[j] -= mult * B[i];
		}
	}
	vector<double> res(n);
	for (int i = n - 1; i >= 0; i--)
	{
		res[i] = B[i];
		for (int j = i + 1; j < n; j++)
			res[i] -= A[i][j] * res[j];
	}
	return res;
}

double FindingResidualVector(vector<vector<double>> A, vector<double> b, vector<double> answer)
{
	double maxNumber = 0;
	vector<double> errorVector(answer.size());
	for (int i = 0; i < answer.size(); i++)
	{
		for (int j = 0; j < answer.size(); j++)
			errorVector[i] += A[i][j] * answer[j];
		errorVector[i] -= b[i];
	}
	maxNumber = abs(errorVector[0]);
	for (int i = 1; i < errorVector.size(); i++)
	{
		if (maxNumber < abs(errorVector[i]))
			maxNumber = abs(errorVector[i]);
	}
	return  maxNumber;
}

double ErrorCalculation(vector<double> x1, vector<double> x2)
{
	int n = x1.size();
	double calculatedError = 0;
	double firstMaxNumber = 0, secondMaxNumber = 0;
	for (int i = 0; i < n; i++) {
		if (x2[i] - x1[i] > firstMaxNumber)
			firstMaxNumber = x2[i] - x1[i];
		if (x1[i] > secondMaxNumber)
			secondMaxNumber = x1[i];
	}
	calculatedError = firstMaxNumber / secondMaxNumber;
	return calculatedError;
}

void OutputMatrix(vector<vector<double>> A, vector<double> B) 
{
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A.size(); j++)
			cout << A[i][j] << setw(11);
		cout << setw(3) << "| " << B[i] << endl;
	}
}

int main()
{
	vector<vector<double>> enteredMatrix = { {0.14,0.24,-0.84},
											  {1.07,-0.83,0.56} ,
											{0.64,0.43,-0.38} };
	vector<double> extendedPart = { 1.11,0.48,-0.83};
	cout << "Entered matrix: \n";
	OutputMatrix(enteredMatrix, extendedPart);
	cout << "\nAnswer: \n";
	vector<double> firstResult = GaussSolution(enteredMatrix, extendedPart);
	for (int i = 0; i < firstResult.size(); i++)
		cout << "x" << i + 1 << "= " << firstResult[i] << " ";
	double normaRV = FindingResidualVector(enteredMatrix, extendedPart, firstResult);
	cout << "\n\nNorma of residual vector: \n" << normaRV << endl;
	vector<double> secondResult = GaussSolution(enteredMatrix, enteredMatrix * firstResult);
	cout << "\nSecond solution: \n";
	for (int i = 0; i < secondResult.size(); i++)
		cout << "x" << i + 1 << "= " << secondResult[i] << " ";
	double calculatedError = ErrorCalculation(firstResult, secondResult);
	cout << "\n\nRelative error estimate: " << calculatedError << endl;
	return 0;
}
