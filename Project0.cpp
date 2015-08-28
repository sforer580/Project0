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
int q = 1;
class soln;
double replicate(soln& S, soln& S2);


// creates solution class
class soln
{
public:
	vector <double> coeff;
	double fitness;
	double test(double, double);
	double calc_fit(double);
	double ans_f_approx;
	double ans_f_target;
	double diff;
	double get_fitness();
	double coefficients();
	double get_x_val(double);
	double x_val;
	double mutation();
	double co;
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
	diff = fabs(ans_f_approx - ans_f_target);
return diff;
}

//calculates the fitness function
double soln::calc_fit(double diff)
{
	fitness = -diff;
return fitness;
}

double soln::coefficients()
{
	for (int i = 0; i < number_coeff; i++)
	{
		co = ((double)rand() / RAND_MAX) * 0.001;		//creates random coeff between 0 and 0.001
		coeff.push_back(co);
	}
	for (int i = 0; i < coeff.size(); i++)
	{
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
	return 1;
}

double soln::get_x_val(double x_value)
{
	for (int i = 0; i < 1; i++)
	{
		x_val = x_value;			//creates random x value between 0 and 2pi
		//cout << x_val << "\t" << endl;
	}
	return x_val;
}

double soln::get_fitness()
{
	ans_f_approx = f_approx(f_1(a, b, x_val), f_2(c, d, x_val));
	ans_f_target = f_target(x_val);
	cout << "Begin function calculations" << endl;
	cout << f_1(a, b, x_val) << endl;							//writes function 1 value
	cout << f_2(c, d, x_val) << endl;							//writes function 2 value
	cout << ans_f_approx << endl;								//writes approximate function value
	cout << ans_f_target << endl;								//writes target function value
	cout << test(ans_f_approx, ans_f_target) << endl;			//writes difference betweeen the approximate function and target function values
	cout << calc_fit(diff) << endl;								//writes the fitness
	//cout << "\n";
	
	return fitness;
}


double replicate(soln& S, soln& S2)
{
	cout << "replicate in" << endl;
	if (S.fitness < S2.fitness)
	{
		S = S2;
		cout << "case 1" << endl;
	}
	else
	{
		S2 = S;
		cout << "case 2" << endl;
	}
	cout << S.fitness << endl;
	cout << S2.fitness << endl;
	cout << "replicate out" << endl;
	return 1;
}

double soln::mutation()
{
	cout << "mutation in" << endl;
	cout << coeff.size() << endl;
	for (int i = 0; i < number_coeff; i++)
	{
		//coeff.erase(coeff.begin()+3);
		coeff.at(i) = coeff.at(i) + ((((double)rand() / RAND_MAX) * 0.001) - (((double)rand() / RAND_MAX) * 0.001));		//creates random coeff between 0 and 0.001
	}
	cout << coeff.size() << endl;
	cout << "mutation out" << endl;
	return 1;
}




int main()
{
	soln S;
	soln S2;
	srand(time(NULL));
	number_coeff = num_coeff();
	cout << "The number of coefficients is ";
	cout << number_coeff << endl;								//displays number of coeff based on number of prim functions
	
		double xvalue = ((double)rand() / RAND_MAX) * (2 * PI);
		cout << S.coefficients() << endl;
		cout << S.get_x_val(xvalue) << endl;
		S.get_fitness();							//displays fitness
		cout << "\n" << endl;
		cout << "next set" << endl;
		//cout << S.coeff.size() << endl;
		cout << S2.coefficients() << endl;
		cout << S2.get_x_val(xvalue) << endl;
		S2.get_fitness();						//displays fitness
		//cout << S2.coeff.size() << endl;
		cout << "\n" << endl;
		cout << S.fitness << endl;
		cout << S2.fitness << endl;
		//Determines which solution stays and replicates
		cout << replicate(S, S2) << endl;

		for (int i = 0; i < 1; i++)
		{
			cout << S2.mutation() << endl;
			cout << "\n" << endl;
			cout << "next set" << endl;
			cout << S2.coeff.size() << endl;
			cout << S2.coefficients() << endl;
			cout << S2.get_x_val(xvalue) << endl;
			S2.get_fitness();						//displays fitness
													//cout << S2.coeff.size() << endl;
			cout << "\n" << endl;
			cout << S.fitness << endl;
			cout << S2.fitness << endl;
			//Determines which solution stays and replicates
			cout << replicate(S, S2) << endl;
		}

   return 0;
}