// Project0.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
using namespace std;


// intialize functions
double f_1(double a, double b, double x);
double f_2(double a, double b, double x);
double f_approx(double a, double b, double x);
double f_target(double x);
int num_coeff(int num_func);
double PI = 3.1415926535897;


// creates solution class
class soln
{
public:
	vector <double> coeff;
	double fitt;
	double test(double, double);
	double calc_fitt(double);

};

int num_coeff()
{
	int num_func;
	cout << "Enter number of primitaive functions" << endl;
	cin >> num_func;
	int num_coeff = num_func * 2;
	return num_coeff;
}


double f_1(double a, double b, double x)
{
	double ans_f_1 = a * (x * x) + b;
	return ans_f_1;
}
double f_2(double a, double b, double x)
{
	double ans_f_2 = a * sin(b * x);
	return ans_f_2;
}
double f_approx(double ans_f_1, double ans_f_2, double x)
{
	double ans_f_approx = ans_f_1 + ans_f_2;
	return ans_f_approx;
}


double soln::test(double approx, double target)
{
	double diff;
	diff = approx - target;
return diff;
}

double soln::calc_fitt(double diff)
{
	double fittness;
	fittness = diff;
return fittness;
}


int main()
{
	int number_coeff = num_coeff();
	cout << "The number of coefficients is ";
	cout << number_coeff << endl;				//displays number of coeff based on number of prim functions
	soln S;
	for (int i = 0; i < number_coeff; i++)
	{
		double a = ((double)rand() / RAND_MAX) * 0.001;		//creates random coeff between 0 and 0.001
		S.coeff.push_back(a);
		cout << S.coeff.at(i) << "\t";
	}
	for (int i = 0; i < 1; i++)
	{
		double a = ((double)rand() / RAND_MAX) * (2*PI);		//creates random coeff between 0 and 2pi
		S.?//.push_back(a);
		cout << S.//.at(i) << "\t";
	}
	cout << "\n";
	cout << S.coeff.empty() << endl;
	cout << S.coeff.size() << endl;
	cout << S.calc_fitt << endl;
   return 0;
}

//while (f_approx =! f_target) 
//{
//	
//	
//	
//	
//	
//	
//	f_approx = f_1 + f_2;
//
//
//
//
//
//}
//
//
//}
//
////double build_vector_a(double)
////
////double build_vector_a()
////{
////vector<double> vector_a; 
////	vector_a.push_back(double a)
////	
////}
////vector<double> vector_b();
////{
////cout << "enter values for vector b" << endl;
////cin >> vector_b endl;
////}
////
////
////double f_1(double a, double b, double x)
////{
////	double answer_f_1 = a * (x ^ 2) + b;
////}
////double f_2(double a, double b, double x)
////{
////	double answer_f_2 = a * sin(b * x);
////}
////double f_approx(double a, double b, double x)
////{
////
////}
////double f_target(double a, double b, double x)
////{
////	double answer_f_target = 1 * sin(1 * x);
////}
////