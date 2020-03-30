#include <iostream>
#include <math.h>
#include <conio.h>

using namespace std;

double * divideIntervalEqually(int a, int b, int n)
{
	double * points = new double[n + 1];
	for (int i = 0; i <= n; i++)
		points[i] = a + i * ((b - a) / double(n));
	return points;
}

double * divideIntervalChebyshev(int a, int b, int n)
{
	double * points = new double[n + 1];
	for (int i = n; i >= 0; i--)
		points[n - i] = ((b - a) * cos(((2 * i + 1)*3.14159265358) / (2 * n + 2)) + (a + b)) / 2;
	return points;
}

double f(double x)
{
	return abs(x*x - 1);
}

double dividedDifference(double points[], int k)
{
	if (k == 1)
		return f(points[0]);
	double * pointsRight = new double[k - 1];
	double * pointsLeft = new double[k - 1];
	for (int i = 0; i < k; i++)
	{
		if (i != 0)
			pointsRight[i - 1] = points[i];		// All points without the first one
		if (i != k - 1)
			pointsLeft[i] = points[i];			// All points without the last one
	}
	return (dividedDifference(pointsRight, k - 1) - dividedDifference(pointsLeft, k - 1)) / (points[k - 1] - points[0]);
}

double * getNewtonCoefficients(double points[], int n)
{
	double * coefficients = new double[n + 1];
	double * pointsSubset;
	for (int i = 0; i <= n; i++)		// i represents the index of the last point taken into calculation of given coefficient, x0<, x0 x1<, x0 x1 x2<, etc.
	{
		pointsSubset = new double[i + 1];
		for (int j = 0; j <= i; j++)
			pointsSubset[j] = points[j];

		coefficients[i] = dividedDifference(pointsSubset, i + 1);
	}
	return coefficients;
}

double getNewtonPolynomialValueAt(double coefficients[], double points[], int n, double x)
{
	double polynomialValue = 0;
	double polynomialElement;
	for (int i = 0; i <= n; i++)		// i represents the index of the last point taken into calculation of given coefficient, x0<, x0 x1<, x0 x1 x2<, etc.
	{
		polynomialElement = coefficients[i];
		for (int j = 0; j < i; j++)
			polynomialElement *= (x - points[j]);			// f[x0, ..., xi] *= (x - xj), j < i

		polynomialValue += polynomialElement;
	}
	return polynomialValue;
}

double getLagrangePolynomialValueAt(double points[], int n, double x)
{
	double polynomialValue = 0;
	double auxilaryPolynomial;
	for (int i = 0; i <= n; i++)
	{
		auxilaryPolynomial = 1;
		for (int j = 0; j <= n; j++)
		{
			if (i == j)
				continue;
			auxilaryPolynomial *= (x - points[j]);
			auxilaryPolynomial /= (points[i] - points[j]);
		}
		auxilaryPolynomial *= f(points[i]);
		polynomialValue += auxilaryPolynomial;
	}
	return polynomialValue;
}

double estimateMaxErrorNewton(int a, int b, double activePoints[], int n, double coefficients[])
{
	double maxError = 0;
	for (double i = a; i <= b; i += 0.001)
		if (abs(getNewtonPolynomialValueAt(coefficients, activePoints, n, i) - f(i)) > maxError)
			maxError = abs(getNewtonPolynomialValueAt(coefficients, activePoints, n, i) - f(i));
	return maxError;
}

double estimateMaxErrorLagrange(int a, int b, double activePoints[], int n)
{
	double maxError = 0;
	for (double i = a; i <= b; i += 0.001)
		if (abs(getLagrangePolynomialValueAt(activePoints, n, i) - f(i)) > maxError)
			maxError = abs(getLagrangePolynomialValueAt(activePoints, n, i) - f(i));
	return maxError;
}

void initializeParameters(int & n, double & a, double & b)
{
	system("cls");
	cout << "Interpolation" << endl;
	cout << "Enter the interpolant degree: ";
	cin >> n;
	cout << "Enter the interpolating interval lower bound: ";
	cin >> a;
	cout << "Enter the interpolating interval upper bound: ";
	cin >> b;
	cout << endl;
}

int main()
{
	int n;
	double a, b;
	initializeParameters(n, a, b);

	double * eqDistPoints = divideIntervalEqually(a, b, n);
	double * chebyshevPoints = divideIntervalChebyshev(a, b, n);

	double * activePoints = new double[1];		// initialization only in conditional statements is illegal to the compiler
	bool run = true;
	bool repeat;

	while (run)
	{
		system("cls");
		repeat = false;
		cout << "[a, b] = [" << a << ", " << b << "], n = " << n << endl;
		cout << "1. Equally distributed nodes\n2. Chebyshev nodes\n3. Change parameters\n> ";
		switch (_getch())
		{
		case 49:
			activePoints = eqDistPoints;
			cout << "Equally distributed nodes <";
			break;
		case 50:
			activePoints = chebyshevPoints;
			cout << "Chebyshev nodes <";
			break;
		case 51:
			system("cls");
			initializeParameters(n, a, b);
			eqDistPoints = divideIntervalEqually(a, b, n);
			chebyshevPoints = divideIntervalChebyshev(a, b, n);
			repeat = true;
			break;
		default: repeat = true;
			break;
		case 27: run = false;
			break;
		}
		if (!run || repeat)
			continue;

		cout << "\n[";
		for (int i = 0; i < n; i++)
			cout << activePoints[i] << ", ";
		cout << activePoints[n] << "]\n";

		double * newtonCoefficients = getNewtonCoefficients(activePoints, n);

		double arg;
		do
		{
			cout << "\nEnter the argument for wanted interpolant value: ";
			cin >> arg;
		} while (arg < a || arg > b);
		cout << "Lagrange Interpolating Polynomial value at " << arg << " is " << getLagrangePolynomialValueAt(activePoints, n, arg) << endl;
		cout << "Newton Interpolating Polynomial value at " << arg << " is " << getNewtonPolynomialValueAt(newtonCoefficients, activePoints, n, arg) << endl;


		cout << "Max error (Lagrange): " << estimateMaxErrorLagrange(a, b, activePoints, n) << endl;
		cout << "Max error (Newton): " << estimateMaxErrorNewton(a, b, activePoints, n, newtonCoefficients) << endl;
		system("pause");
	}

	return 0;
}