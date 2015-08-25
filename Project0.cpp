// Project0.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <time.h>
using namespace std;


// intialize functions
double f_1(double a, double b, double x);
double f_2(double c, double d, double x);
double f_approx(double e, double f, double x);
double f_target(double x);
int num_coeff(int num_func);
int number_coeff;
double PI = 3.1415926535897;
double a;
double b;
double c;
double d;


// creates solution class
class soln
{
public:
	vector <double> coeff;
	double x_val;
	double fitness;
	double test(double, double);
	double calc_fit(double);
	double ans_f_approx;
	double ans_f_target;
	double diff;
	double get_fitness();
};

//calculates the number of coeffs based on the number of primitaive functions
int num_coeff()
{
	int num_func;
	cout << "Enter number of primitaive functions" << endl;
	cin >> num_func;
	int num_coeff = num_func * 2;
	return num_coeff;
}

//function 1
double f_1(double a, double b, double x)
{
	double ans_f_1 = a * (x * x) + b;
	return ans_f_1;
}

//function 2
double f_2(double c, double d, double x)
{
	double ans_f_2 = c * sin(d * x);
	return ans_f_2;
}

//function approximation
double f_approx(double ans_f_1, double ans_f_2)
{
	double ans_f_approx = ans_f_1 + ans_f_2;
	return ans_f_approx;
}

//target function calculation
double f_target(double x)
{
	double ans_f_target = 1 * sin(1 * x);
	return ans_f_target;
}

//calculates the difference between the approximate function and target function
double soln::test(double ans_f_approx, double ans_f_target)
{
	diff = ans_f_approx - ans_f_target;
return diff;
}

//calculates the fitness function
double soln::calc_fit(double diff)
{
	fitness = -diff;
return fitness;
}

double soln::get_fitness()
{
	for (int i = 0; i < number_coeff; i++)
	{
		double co = ((double)rand() / RAND_MAX) * 0.001;		//creates random coeff between 0 and 0.001
		coeff.push_back(co);
		cout << coeff.at(i) << "\t";
	}
	
	cout << "\n";
	//cout << coeff.empty() << endl;
	//cout << coeff.size() << endl;
	a = coeff.at(0);
	b = coeff.at(1);
	c = coeff.at(2);
	d = coeff.at(3);
	//cout << "check1" << endl;
	for (int i = 0; i < 1; i++)
	{
		x_val = ((double)rand() / RAND_MAX) * (2 * PI);			//creates random x value between 0 and 2pi
		cout << x_val << "\t" << endl;
	}
	ans_f_approx = f_approx(f_1(a, b, x_val), f_2(c, d, x_val));
	ans_f_target = f_target(x_val);
	cout << "Begin function calculations";
	cout << f_1(a, b, x_val) << endl;							//writes function 1 value
	cout << f_2(c, d, x_val) << endl;							//writes function 2 value
	cout << ans_f_approx << endl;								//writes approximate function value
	cout << ans_f_target << endl;								//writes target function value
	cout << test(ans_f_approx, ans_f_target) << endl;			//writes difference betweeen the approximate function and target function values
	cout << calc_fit(diff) << endl;								//writes the fitness
	cout << "\n";
	
	return fitness;
}


int main()
{
	soln S;
	srand(time(NULL));
	number_coeff = num_coeff();
	cout << "The number of coefficients is ";
	cout << number_coeff << endl;								//displays number of coeff based on number of prim functions
	cout << S.get_fitness() << endl;							//displays fitness

   return 0;
}