#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

void OutputMatrixInConsole(vector<vector<double>> A, vector<double> b);

enum UserChoice
{
	manualEntryMatrix = 1,
	defaultMatrix = 2,
};

bool CheckRowProportionality(vector<double> firstVector, vector<double> secondVector) 
{
    int size = firstVector.size();
    double ram = firstVector[0] / secondVector[0];
    for (int i = 1; i < size; i++) 
	{
        if (ram == firstVector[i] / secondVector[i])
            continue;
        else
            return false;
    }
    return true;
}

vector<double> MultiplyMatrixVector(vector<vector<double>> matrix, vector<double> vector_)
{
	vector<double> resultVector(matrix.size());
	for (int i = 0; i < resultVector.size(); i++) 
	{
		resultVector[i] = 0;
		for (int j = 0; j < resultVector.size(); j++) 
			resultVector[i] += matrix[i][j] * vector_[j];
	}
	return resultVector;
}

bool CheckProportionality(vector<vector<double>> A, vector<double> b)
{
	int vectorSize = A.size();
	vector<vector<double>> A_b(vectorSize, vector<double>(vectorSize + 1, 0));
	for (int i = 0; i < vectorSize; i++) 
	{
		for (int j = 0; j < vectorSize + 1; j++) 
		{
			if (j == vectorSize)
				A_b[i][j] = b[i];
			else
				A_b[i][j] = A[i][j];
		}
	}
	bool proportional = false;
	for (int i = 0; i < vectorSize - 1; i++)
	{
		for (int j = i + 1; j < vectorSize; j++)
		{
			proportional = CheckRowProportionality(A_b[i],A_b[j]);
			if (proportional)
				return true;
		}
	}
	return false;
}

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

double FindMaxInRV(vector<vector<double>> A, vector<double> b, vector<double> firstRoots)
{
	vector<double> vectorCalculationDiff(firstRoots.size());
	double maxInVCD = 0;
	for (int i = 0; i < firstRoots.size(); i++)
	{
		for (int j = 0; j < firstRoots.size(); j++)
			vectorCalculationDiff[i] += A[i][j] * firstRoots[j];
		vectorCalculationDiff[i] -= b[i];
	}
	maxInVCD = abs(vectorCalculationDiff[0]);
	for (int i = 1; i < vectorCalculationDiff.size(); i++)
	{
		if (maxInVCD < abs(vectorCalculationDiff[i]))
			maxInVCD = abs(vectorCalculationDiff[i]);
	}
	return  maxInVCD;
}

double CalculateError(vector<double> firstRoots, vector<double> x2)
{
	int vectorSize = firstRoots.size();
	double calcError = 0;
	double max1 = 0, max2 = 0;
	for (int i = 0; i < vectorSize; i++) 
	{
		if (x2[i] - firstRoots[i] > max1)
			max1 = x2[i] - firstRoots[i];
		if (firstRoots[i] > max2)
			max2 = firstRoots[i];
	}
	calcError = max1 / max2;
	return calcError;
}

void OutputMatrixInConsole(vector<vector<double>> A, vector<double> b)
{
	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A.size(); j++)
			cout << A[i][j] << setw(11);
		cout << setw(3) << "| " << b[i] << endl;
	}
}

int main()
{
	vector<vector<double>> A;
	vector<double> b;
	int userChoice;
	cout << "how do you want to initialize the matrix?\n"
		<< (int)UserChoice::manualEntryMatrix << " - enter the matrix manually\n"
		<< (int)UserChoice::defaultMatrix << " - use an already defined matrix\n"
		<< "any other value terminates the program!\n\nyour answer: ";
	cin >> userChoice;

	switch ((UserChoice)userChoice)
	{
	case UserChoice::manualEntryMatrix:
	{
		cout << "enter matrix order: ";
		int matrixSize;
		cin >> matrixSize;
		A.assign(matrixSize, vector<double>(matrixSize));
		b.resize(matrixSize);
		for (int i = 0; i < matrixSize; i++)
		{
			for (int j = 0; j < matrixSize; j++)
			{
				printf("enter element with index A[%d][%d]: ", i + 1, j + 1);
				cin >> A[i][j];
			}
		}
		for (int i = 0; i < matrixSize; i++)
		{
			printf("enter element b[%d]: ", i + 1);
			cin >> b[i];
		}
		break;
	}

	case UserChoice::defaultMatrix:
	{
		A = { {0.14,0.24,-0.84},
			{1.07,-0.83,0.56},
			{0.64,0.43,-0.38} };
		b = { 1.11,0.48,-0.83 };
		break;
	}

	default:
		cout << "the program has ended\n";
		return 0;
	}

	cout << "\nEntered matrix: \n";
	OutputMatrixInConsole(A, b);
	if (CheckProportionality(A, b))
	{
		cout << "rows are proportional, the roots of the system cannot be calculated";
		return 0;
	}
	cout << "\nAnswer: \n";
	vector<double> firstRoots = SolveGaussMethod(A, b);
	for (int i = 0; i < firstRoots.size(); i++)
		printf("x%d= %f  ", i + 1, firstRoots[i]);
	double maxNormaOfRV = FindMaxInRV(A, b, firstRoots);
	cout << "\n\nNorma of residual vector: \n" << maxNormaOfRV << endl;
	vector<double> result2 = SolveGaussMethod(A, MultiplyMatrixVector(A,firstRoots));
	cout << "\nSecond solution: \n";
	for (int i = 0; i < result2.size(); i++)
		printf("x%d= %f ", i + 1, result2[i]);
	cout << "\n\nCalculation error: \n";
	double calculationError = CalculateError(firstRoots, result2);
	cout << calculationError << endl;
	return 0;
}
