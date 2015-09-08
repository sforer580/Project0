// Project0.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdio.h>
using namespace std;


//Intialize functions and global variables
int num_coeff(int num_func);
int number_coeff;
double PI = 3.1415926535897;
int n = 100;
int samples = 5;
double f_1(double a, double b, double x);
double f_2(double c, double d, double x);
double f_approx(double e, double f, double x);
double f_target(double x);
class soln;
double replicate(soln& S, soln& S2);


//Creates solution class
class soln
{
public:
	vector <double> coeff;
	double coefficients();
	double co;
	double check_coeff();
	double assign_coeff();
	double a;
	double b;
	double c;
	double d;
	double create_x_value(int);
	double set_x_value(double q);
	void interp_x_value(int index);
	//double get_x_val(double);
	double x_value;
	double x_val;
	double ans_f_approx;
	double ans_f_target;
	double diff;
	double test(double, double);
	double calc_fit(double);
	double get_fitness();
	double fitness;
	void reset_fitness();
	double mutation();
};


//Calculates the number of coeffs based on the number of primitive functions
int num_coeff()
{
	int num_func;
	cout << "Enter number of primitive functions" << endl;
	cin >> num_func;
	int num_coeff = num_func * 2;
	return num_coeff;
}


//Creates random coefficients between 0 and 0.001
double soln::coefficients()
{
	for (int i = 0; i < number_coeff; i++)
	{
		co = ((double)rand() / RAND_MAX) * 0.001+1;
		coeff.push_back(co);
	}
	return 1;
}


//Displays coefficients
double soln::check_coeff()
{
	for (int i = 0; i < coeff.size(); i++)
	{
		cout << coeff.at(i) << "\t";
	}
	cout << "\n" << endl;
	return 1;
}


//Assigns coefficients to their respective index
double soln::assign_coeff()
{
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


//Creates x value between 0 and 2pi
double soln::create_x_value(int i)
{
	double X = -1.0;
	static int first = 0;
	vector<double> Xstart;
	if (first == 0)
	{
		for (int z = 0; z < samples; z++)
		{
			Xstart.push_back(((double)rand() / RAND_MAX)*(2*PI));
		}
		first = 1;
	}
	const static vector<double> x_value = Xstart;
	return x_value.at(i);
}



//Sets the x value 
//double soln::set_x_value(double q)
//{
//	x_value = q;
//	return 1;
//}


//interp x values
void soln::interp_x_value(int index)
{
	//take in an index
	//interp x value at that index
	x_val = create_x_value(index);
	cout << x_val << endl;
}

//Assings the x value to a variable so that both instances of soln class will be evaluated at the same point
//double soln::get_x_val(double x_value)
//{
//	for (int i = 0; i < 1; i++)
//	{
//		x_val = x_value;
//	//cout << x_val << "\t" << endl;
//	}
//	return x_val;
//}

//Function 1
double f_1(double a, double b, double x)
{
	double ans_f_1 = a * (x * x) + b;
	return ans_f_1;
}


//Function 2
double f_2(double c, double d, double x)
{
	double ans_f_2 = c * sin(d * x);
	return ans_f_2;
}


//Function approximation
double f_approx(double ans_f_1, double ans_f_2)
{
	double ans_f_approx = ans_f_1 + ans_f_2;
	return ans_f_approx;
}


//Target function calculation
double f_target(double x)
{
	double ans_f_target = 1 * sin(1 * x);
	return ans_f_target;
}


//Calculates the difference between the approximate function and target function
double soln::test(double ans_f_approx, double ans_f_target)
{
	diff = fabs(ans_f_approx - ans_f_target);
return diff;
}


//Calculates the fitness function
double soln::calc_fit(double diff)
{
	fitness = fitness - diff;
return fitness;
}


//Calculates each function and the fitness
double soln::get_fitness()
{
	ans_f_approx = f_approx(f_1(a, b, x_val), f_2(c, d, x_val));
	ans_f_target = f_target(x_val);
	cout << "Begin function calculations" << endl;
	cout << "The calculated value for function 1 is " << endl;
	cout << f_1(a, b, x_val) << endl;							//writes function 1 value
	cout << "The calculated value for function 2 is " << endl;
	cout << f_2(c, d, x_val) << endl;							//writes function 2 value
	cout << "The calculated approximation is " << endl;
	cout << ans_f_approx << endl;								//writes approximate function value
	cout << "The calculated target fucntion is " << endl;
	cout << ans_f_target << endl;								//writes target function value
	test(ans_f_approx, ans_f_target);							//writes difference betweeen the approximate function and target function values
	cout << "The fitness is " << endl;
	cout << calc_fit(diff) << endl;								//writes the fitness
	//cout << "\n";
	return 1;
}


//Resets the fitness values for both instacne of the soln class
void soln::reset_fitness()
{
	fitness = 0;
}


//Decides which instance of the soln class survives and then replicates the surviver
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


//Mutates the coefficients of the replication of the surviver
double soln::mutation()
{
	cout << "mutation in" << endl;
	//cout << coeff.size() << endl;
	for (int i = 0; i < number_coeff; i++)
	{
		double range = 0.1;
		//coeff.erase(coeff.begin()+3);
		coeff.at(i) = coeff.at(i) + ((((double)rand() / RAND_MAX) * range) - (((double)rand() / RAND_MAX) * range));		//creates random coeff between 0 and 0.001
	}
	for (int i = 0; i < coeff.size(); i++)
	{
		cout << coeff.at(i) << "\t";
	}
	cout << "\n" << endl;
	//cout << coeff.size() << endl;
	cout << "mutation out" << endl;
	return 1;
}


int main()
{
	soln S;															//Creates the first instance of soln class
	soln S2;														//Creates the second instance of soln class
	srand(time(NULL));
	number_coeff = num_coeff();										//Assings a varibale to the output of the num_coeff function
	cout << "The number of coefficients is ";
	cout << number_coeff << endl;									//displays number of coeff based on number of prim functions
	S.coefficients();												//Creates random coefficients the first instance of soln class
	S2.coefficients();												//Creates random coefficients the second instance of soln class
	FILE * pfile;
	errno_t err;
	err = fopen_s(&pfile, "test.txt", "w");
	for (int i = 0; i < n; i++)
	{
		cout << "First instance coefficients" << endl;
		S.check_coeff();											//Displays coefficients the first instance of soln class
		cout << "Second instance coefficients" << endl;
		S2.check_coeff();											//Displys coefficients the second instance of soln class
		S.assign_coeff();											//Assigns coefficients to their respective index in the first instance of soln class
		S2.assign_coeff();											//Assigns coefficients to their respective index in the second instance of soln class
		S.reset_fitness();
		S2.reset_fitness();
		for (int i = 0; i < samples; i++)
		{
			S.create_x_value(i);
			cout << "First instance of soln class" << endl;
			cout << "The x value is " << endl;
			S.interp_x_value(i);					//Gets the x value and assigns it to the first instance of soln class
			S.get_fitness();										//Calculates and displays fitness for the first instance of soln class
			cout << "\n" << endl;
			cout << "Second instance of soln class" << endl;
			cout << "The x value is " << endl;
			S2.interp_x_value(i);					//Gets the x value and assigns it to the second instance of soln class
			S2.get_fitness();										//Calculates and displays fitness for the second instance of soln class
			cout << "\n" << endl;
		}
		fprintf(pfile, "%.5f\t%.5f\n", S.fitness, S2.fitness);
		replicate(S, S2);											//Determines which solution stays and replicates
		cout << "\n" << endl;
		S2.mutation();												//Mutates the losing instance of the soln class
		cout << "\n" << endl;
		cout << "next set" << endl;
		
	}
	cout << "Solution coefficients" << endl;
	S.check_coeff();												//Displays the solution coefficients
	fclose(pfile);
	return 0;
}

//two fucntions for putting the target fuction to a text file