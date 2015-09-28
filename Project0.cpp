// Project0.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <algorithm>
using namespace std;


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Intialize functions and global variables
int num_coeff(int num_func);
int number_coeff;
double PI = 3.1415926535897;
int max_generations = 300;
int max_samples = 100;
double f_1(double a, double b, double x);
double f_2(double c, double d, double x);
double f_approx(double e, double f, double x);
double f_target(double x);
class soln;
double replicate(soln& S, soln& S2);
void generation(vector <soln>* pPop);
int binary_select(vector<soln>* pPop);
void display_x_values(FILE*pFile2);
class scoreboard;


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Creates solution class
class soln
{
public:
	vector <double> coeff;
	void coefficients();
	double co;
	double check_coeff();
	double assign_coeff();
	double a;
	double b;
	double c;
	double d;
	double create_x_value(int);
	void interp_x_value(int index);
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
	int age;
	vector<double> get_final_coefficients();
	void plot_to_file3(FILE* pFile3);
};


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//creates the scoreboard class
class scoreboard
{
public:
	void reorder_solutions(vector<soln>* pPop);
	//void obtain_solns(vector<soln>* pPop);
	//void sort(less_than_fitness());
};


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Calculates the number of coefficients based on the number of primitive functions
int num_coeff()
{
	int num_func;
	cout << "Enter number of primitive functions" << endl;
	cin >> num_func;
	int num_coeff = num_func * 2;
	return num_coeff;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Creates random coefficients between 0 and 0.001
void soln::coefficients()
{
	age = 0;
	for (int i = 0; i < number_coeff; i++)
	{
		co = ((double)rand() / RAND_MAX) * 0.001+1;
		coeff.push_back(co);
	}
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Assigns coefficients to their respective index
double soln::assign_coeff()
{
	//cout << "\n";
	a = coeff.at(0);
	b = coeff.at(1);
	c = coeff.at(2);
	d = coeff.at(3);
	return 1;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Creates static x value between 0 and 2pi
double soln::create_x_value(int i)
{
	static int first = 0;
	vector<double> Xstart;
	if (first == 0)
	{
		for (int z = 0; z < max_samples; z++)
		{
			Xstart.push_back(((double)rand() / RAND_MAX)*(2*PI));
		}
		first = 1;
	}
	const static vector<double> x_value = Xstart;
	return x_value.at(i);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//interp x values
void soln::interp_x_value(int index)
{
	x_val = create_x_value(index);
	//cout << x_val << endl;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Function 1
double f_1(double a, double b, double x)
{
	double ans_f_1 = a * (x) + b;
	return ans_f_1;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Function 2
double f_2(double c, double d, double x)
{
	double ans_f_2 = c * sin(d * x);
	return ans_f_2;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Function approximation
double f_approx(double ans_f_1, double ans_f_2)
{
	double ans_f_approx = ans_f_1 + ans_f_2;
	return ans_f_approx;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Target function calculation
double f_target(double x)
{
	double ans_f_target = 1 * sin(1 * x);
	return ans_f_target;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Calculates the difference between the approximate function and target function
double soln::test(double ans_f_approx, double ans_f_target)
{
	diff = fabs(ans_f_approx - ans_f_target);
return diff;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Calculates the fitness based on the total fitness for each x value for each solution
double soln::calc_fit(double diff)
{
	fitness = fitness - diff;
return fitness;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Calculates each function and the fitness
double soln::get_fitness()
{
	ans_f_approx = f_approx(f_1(a, b, x_val), f_2(c, d, x_val));
	ans_f_target = f_target(x_val);
	f_1(a, b, x_val);							//writes function 1 value
	f_2(c, d, x_val);							//writes function 2 value
	ans_f_approx;								//writes approximate function value
	ans_f_target;								//writes target function value
	test(ans_f_approx, ans_f_target);							//writes difference betweeen the approximate function and target function values
	calc_fit(diff);								//writes the fitness
	return 1;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Resets the fitness values for each instacne of the soln class
void soln::reset_fitness()
{
	fitness = 0;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Performs calculations to determine the outcome of each generation
void generation(vector <soln>* pPop)
{
	for (int mem = 0; mem < pPop->size(); mem++)
	{
		pPop->at(mem).age++;
		pPop->at(mem).assign_coeff();										//assigns the coefficients to the fucntions
		pPop->at(mem).reset_fitness();										//resets the fitness
		for (int samples = 0; samples < max_samples; samples++)
		{
			pPop->at(samples).interp_x_value(samples);						//interprets the x value
			pPop->at(samples).get_fitness();								//get the fitness for each soln at each x value
		}
	}
	// Downselect
	int to_kill = pPop->size() / 2;
	for (int i = 0; i < to_kill; i++)
	{
		int kill;
		kill = binary_select(pPop);
		pPop->erase(pPop->begin() + kill);
	}
	//replicates the surviving solutions and then mutates them
	int to_replicate = to_kill;
	for (int i = 0; i < to_replicate; i++)
	{
		soln A;
		int spot = rand() % pPop->size();
		A = pPop->at(spot);
			A.mutation();
		pPop->push_back(A);
	}
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Radomly comapares two solutions and kills one until half the solutions have been killed off
int binary_select(vector<soln>* pPop)
{
	int loser;
	int index_1 = rand() % pPop->size();
	int index_2 = rand() % pPop->size();
	double fit_1 = pPop->at(index_1).fitness;
	double fit_2 = pPop->at(index_2).fitness;
	if (fit_1 > fit_2) 
	{
		//then fit 1 survives
		 loser = index_2;
	}
	else
	{
		loser = index_1;
	}
	return loser;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Mutates the coefficients of the replication of the survivers
double soln::mutation()
{
	age = 0;
	for (int i = 0; i < number_coeff; i++)
	{
		if (rand() % 2)
		{
			continue;
		}
		double range = 0.1;
		coeff.at(i) = coeff.at(i) + ((((double)rand() / RAND_MAX) * range) - (((double)rand() / RAND_MAX) * range));		//creates random coeff between 0 and range
	}
	return 1;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//
struct less_than_fitness
{
	inline bool operator() (const soln& A1, const soln& A2)
	{
		return (A1.fitness > A2.fitness);
	}
};

//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//reorders vector of solutions
void scoreboard::reorder_solutions(vector<soln>* pPop)
{
	sort(pPop->begin(), pPop->end(), less_than_fitness());
}





//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//writes the fitness for each generation of the best, median, and worst solutions to a text file
void display_soln_fitness(bool display, vector<soln>* pPop)
{
	FILE*pFile1;
	errno_t err;
	err = fopen_s(&pFile1, "solns_fitness.txt", "w");
	if (display)
	{
		double best_fitness = pPop->at(0).fitness;
		double median_fitness = pPop->at(pPop->size() / 2).fitness;
		double worst_fitness = pPop->at(pPop->size()-1).fitness;
		int best_age = pPop->at(0).age;
		int median_age = pPop->at(pPop->size() / 2).age;
		int worst_age = pPop->at(pPop->size() - 1).age;
		cout << best_age << endl;
		cout << median_age << endl;
		cout << worst_age << endl;
		fprintf(pFile1, "best solution age and fitness \n");
		fprintf(pFile1, "%d\t", best_age);
		fprintf(pFile1, "%.4f\n", best_fitness);
		fprintf(pFile1, "median solution age and fitness \n");
		fprintf(pFile1, "%d\t", median_age);
		fprintf(pFile1, "%.4f\n", median_fitness);
		fprintf(pFile1, "wort solution age and fitness \n");
		fprintf(pFile1, "%d\t", median_age);
		fprintf(pFile1, "%.4f", worst_fitness);
	}
	fclose(pFile1);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//writes the coefficients for the best median and worst solutions to a text file
void display_soln_coefficients(bool display, vector<soln>* pPop)
{
	FILE * pFile2;
	errno_t err;
	err = fopen_s(&pFile2, "soln_coefficients.txt", "w");
	if (display)
	{
		vector<double> best_coeff = pPop->at(0).get_final_coefficients();
		vector<double> median_coeff = pPop->at(pPop->size() / 2).get_final_coefficients();
		vector<double> worst_coeff = pPop->at(pPop->size()-1).get_final_coefficients();
		fprintf(pFile2, "best coefficients \n");
		for (int i = 0; i < best_coeff.size(); i++)
		{
			double best_final_coeff = best_coeff.at(i);
			fprintf(pFile2, "%.4f\t", best_final_coeff);
		}
		fprintf(pFile2, "\n");
		fprintf(pFile2, "median coefficients \n");
		for (int i = 0; i < median_coeff.size(); i++)
		{
			double median_final_coeff = median_coeff.at(i);
			fprintf(pFile2, "%.4f\t", median_final_coeff);
		}
		fprintf(pFile2, "\n");
		fprintf(pFile2, "worst coefficients \n");
		for (int i = 0; i < worst_coeff.size(); i++)
		{
			double worst_final_coeff = worst_coeff.at(i);
			fprintf(pFile2, "%.4f\t", worst_final_coeff);
		}
	}
	fclose(pFile2);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//gets the coefficients for each solution
vector<double> soln::get_final_coefficients()
{
	return coeff;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//writes the y values for the best median and worst solutions to a text file
void display_soln_y_values(bool display, vector<soln>* pPop)
{
	FILE * pFile3;
	errno_t err;
	err = fopen_s(&pFile3, "y-values.txt", "w");
	if (display)
	{
		fprintf(pFile3, "best solution y values \n");
		pPop->at(0).plot_to_file3(pFile3);
		fprintf(pFile3, "\n");
		fprintf(pFile3, "median solution y values \n");
		pPop->at(pPop->size() / 2).plot_to_file3(pFile3);
		fprintf(pFile3, "\n");
		fprintf(pFile3, "worst solution y values \n");
		pPop->at(pPop->size()-1).plot_to_file3(pFile3);
	}
	fclose(pFile3);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//calculates y values for the best median and worst solutions to a text file
void soln::plot_to_file3(FILE* pFile3)
{
	double xx = 0;
	while (xx < 2 * PI)
	{
		double yy = 0;
		yy += f_1(a, b, xx);
		yy += f_2(c, d, xx);
		fprintf(pFile3, "%.2f\t", yy);
		xx += .1;
	}
	fprintf(pFile3, "\n");
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
//writes the static x values to a text file
void display_x_values(bool display, vector<soln>* pPop)
{
	FILE*pFile4;
	errno_t err;
	err = fopen_s(&pFile4, "x-values.txt", "w");
	if (display)
	{
		fprintf(pFile4, "x values \n");
		for (int e = 0; e < max_samples; e++)
		{
			double xx = pPop->at(e).create_x_value(e);
			fprintf(pFile4, "%.2f\t", xx);
		}
	}
	fclose(pFile4);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------
int main()
{
	srand(time(NULL));
	number_coeff = num_coeff();										//Assings a varibale to the output of the num_coeff function
	cout << "The number of coefficients is ";
	cout << number_coeff << endl;									//displays number of coeff based on number of prim functions
	vector<soln> Pop(100);											//population size
	vector<soln>* pPop = &Pop;										//creates a vector of solutions
	scoreboard SB;
	for (int i = 0; i < pPop->size(); i++)
	{
		pPop->at(i).coefficients();									//creates coefficients
	}
	for (int gen = 0; gen < max_generations; gen++)					//calculates the outcome of the generation until the max generations has been reached
	{		
		generation(pPop);
		if (gen % 10 == 0)
		{
			cout << gen << endl;
		}
	}
	SB.reorder_solutions(pPop);
	display_soln_fitness(true, pPop);
	display_soln_coefficients(true, pPop);
	display_soln_y_values(true, pPop);
	display_x_values(true, pPop);
	return 0;
}



//need to plot more data (plot the best median and worst reults)
// write a text file that shows the fitness for each generation