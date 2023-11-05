#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

bool operator / (const vector<double> first_vector, const vector<double> second_vector);
vector<double> operator * (const vector<vector<double>> matrix, const vector <double> vector);
bool checkProportionality(vector<vector<double>> A, vector<double> b);
vector<double> Gauss_solution(vector<vector<double>> A, vector<double> B);
double Finding_the_residual_vector(vector<vector<double>> A, vector<double> b, vector<double> answer);
double Error_calculation(vector<double> x1, vector<double> x2);
void Output_matrix(vector<vector<double>> A, vector<double> B);

bool operator / (const vector<double> first_vector, const vector<double> second_vector)
{
	int size = first_vector.size();
	double ram = first_vector[0] / second_vector[0];
	for (int i = 1; i < size; i++)
		if (ram == first_vector[i] / second_vector[i])
			continue;
		else
			return false;
	return true;
}

vector<double> operator * (const vector<vector<double>> matrix, const vector <double> vector_)
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

bool checkProportionality(vector<vector<double>> A, vector<double> b)
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
	Output_matrix(A, b);
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

vector<double> Gauss_solution(vector<vector<double>> A, vector<double> B)
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

double Finding_the_residual_vector(vector<vector<double>> A, vector<double> b, vector<double> answer)
{
	double max_number = 0;
	vector<double> inter_calculation(answer.size());
	for (int i = 0; i < answer.size(); i++)
	{
		for (int j = 0; j < answer.size(); j++)
			inter_calculation[i] += A[i][j] * answer[j];
		inter_calculation[i] -= b[i];
	}
	max_number = abs(inter_calculation[0]);
	for (int i = 1; i < inter_calculation.size(); i++)
	{
		if (max_number < abs(inter_calculation[i]))
			max_number = abs(inter_calculation[i]);
	}
	return  max_number;
}

double Error_calculation(vector<double> x1, vector<double> x2)
{
	int n = x1.size();
	double calc_error = 0;
	double max1 = 0, max2 = 0;
	for (int i = 0; i < n; i++) {
		if (x2[i] - x1[i] > max1)
			max1 = x2[i] - x1[i];
		if (x1[i] > max2)
			max2 = x1[i];
	}
	calc_error = max1 / max2;
	return calc_error;
}

void Output_matrix(vector<vector<double>> A, vector<double> B) 
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
	vector<vector<double>> entered_matrix = { {0.14,0.24,-0.84},
											  {1.07,-0.83,0.56} ,
											{0.64,0.43,-0.38} };
	vector<double> extended_part = { 1.11,0.48,-0.83};
	cout << "Entered matrix: \n";
	Output_matrix(entered_matrix, extended_part);
	cout << "Answer: \n";
	vector<double> result1 = Gauss_solution(entered_matrix, extended_part);
	for (int i = 0; i < result1.size(); i++)
		cout << "x" << i + 1 << "= " << result1[i] << " ";
	double NormaofRV = Finding_the_residual_vector(entered_matrix, extended_part, result1);
	cout << "\nNorma of residual vector: " << NormaofRV << endl;
	vector<double> result2 = Gauss_solution(entered_matrix, entered_matrix * result1);
	cout << "Second solution: \n";
	for (int i = 0; i < result2.size(); i++)
		cout << "x" << i + 1 << "= " << result2[i] << " ";
	double error = Error_calculation(result1, result2);
	cout << error;
	return 0;
}