#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDKDDKVer.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <complex>
#include "Eigen/Eigenvalues"


#pragma warning(disable : 4996)


void PrintMatrix(const Eigen::MatrixXd &Mat)
{
	std::cout << Mat << '\n';
}

void Eigenproblem(const Eigen::MatrixXd& Mat)
{
	Eigen::EigenSolver<Eigen::MatrixXd> ces;
	ces.compute(Mat);
}

long double MaxEigenValue(const Eigen::MatrixXd& Mat)
{
	Eigen::EigenSolver<Eigen::MatrixXd> ces;
	ces.compute(Mat);

	//double Result = ces.eigenvalues().maxCoeff();
	Eigen::VectorXcd Vec = ces.eigenvalues();

	int j;
	double Result;
	for (int i = 0; i < Vec.cols(); i++)
	{
		if (std::imag(Vec[i]) == 0)
		{
			Result = std::real(Vec[i]);
			break;
		}
		if (i == Vec.cols() - 1) return -1;
	}
	
	for (int i = 0; i < Vec.cols(); i++)
	{
		if (std::imag(Vec[i]) == 0 && 
			(std::real(Vec[i])> Result))
			Result = std::real(Vec[i]);
	}

	return Result;
}

void FillHamiltonian( Eigen::MatrixXd& Mat, double h , double J, double  Jd, double Jt)
{
	double H, S1, S2, S3, S4;
	int i, j;

	Mat.fill(0);

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			Mat(i, j) = 0;
			if (j >= 2) S1 = -1;
			else S1 = 1;
			if (j % 2 == 0) S2 = 1;
			else S2 = -1;

			if (i >= 2) S3 = -1;
			else S3 = 1;
			if (i % 2 == 0) S4 = 1;
			else S4 = -1;


			H = 0;
			H += (J*(S1 * S2) + J*(S2 * S3) + Jd*(S1 * S3) + Jt*(S3 * S4));
			H += -0.5*(S1*0.5 + S2 + S3 + S4*0.5);	
			Mat(i, j) = H;
		}
	}
}

void FillTransferMatrix(const Eigen::MatrixXd& H, Eigen::MatrixXd& Mat, double Beta = 0)
{
	int i, j;

	Mat.fill(0);

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < 4; j++)
		{
			Mat(i, j) = exp(-Beta*H(i,j));
		}
	}
}
void FillTransferMatrix(Eigen::MatrixXd& Mat, double h, double  J, double  Jd, double  Jt, double Beta = 0)
{
	int i, j;

	Mat.fill(0);

	// Ising 2x2
	/*double M11 = exp(Beta * J + Beta * h);
	double M12 = exp(-Beta * J);
	double M21 = exp(-Beta * J);
	double M22 = exp(Beta * J - Beta * h);
	
	Mat(0, 0) = M11;
	Mat(1, 0) = M12;
	Mat(0, 1) = M21;
	Mat(1, 1) = M22;*/

	//triangle

	double M11 = exp(Beta * (2 * J + Jd + Jt + 3 * h))
		+ exp(Beta * (-2 * J + Jd + Jt + h)) 
		+ exp(Beta * (-Jd - Jt + h)) + exp(Beta * (-Jd - Jt - h));

	double M12 = exp(Beta * (2 * J + Jd - Jt + 3 * h))
		+ exp(Beta * (-2 * J + Jd - Jt + h))
		+ exp(Beta * (-Jd + Jt + h)) + exp(Beta * (-Jd + Jt - h));
		
	double M21 = +exp(Beta * (-Jd + Jt + h)) + exp(Beta * (-Jd + Jt - h)) +
		+ exp(Beta * (-2 * J + Jd - Jt - h))
		+ exp(Beta * (2 * J + Jd - Jt - 3 * h));

	double M22 = +exp(Beta * (-Jd - Jt + h)) + exp(Beta * (-Jd - Jt - h)) +
		+ exp(Beta * (-2 * J + Jd + Jt - h))
		+ exp(Beta * (2 * J + Jd + Jt - 3 * h));
		Mat(0, 0) = M11;
		Mat(1, 0) = M12;
		Mat(0, 1) = M21;
		Mat(1, 1) = M22;

		//saw tooth
/*	double M11 = exp(Beta * (2 * J + Jd + 2 * h))
		+ exp(Beta * (- Jd ));
	double M12 = +exp(Beta * (-Jd + h)) + exp(Beta * (-2 * J + Jd - h));

	double M21 = +exp(Beta * (-Jd - h)) + exp(Beta * (-2 * J + Jd + h));

	double M22 = exp(Beta * (2 * J + Jd - 2 * h))
		+ exp(Beta * (-Jd));
	Mat(0, 0) = M11;
	Mat(1, 0) = M12;
	Mat(0, 1) = M21;
	Mat(1, 1) = M22;*/
}

void TextOut(const std::string& Text, int Number)
{
	// #include <stdio.h>
	FILE* stream;

	stream =
		fopen("DataEnh.txt", "a");

	char* p;
	p = (char*)Text.c_str();
	if (stream != NULL)
	{
		fprintf(stream, p, "DataEnh.txt");
		fclose(stream);
	}
	return;
}

double FindMagetization(double FreeNrg, double FreeNrgOld, double h_Step)
{
	// F = - kT ln( lambda_max) free energy
	// M = dF/dh

	double Res = -(FreeNrg - FreeNrgOld)/ h_Step;
	return Res;
}

double CalculateFreeNrg(long double Value, double Beta)
{
	// F = - kT ln( lambda_max) free energy
	long double Res = -(1/Beta) * log(Value);

	return Res;
}

double FindFreeNrg(double h, double Beta, const double h_Step , double J, double  Jd, double Jt)
{
	const int n = 2;
	double FreeNrg;

	Eigen::MatrixXd TransferM(n, n);


	FillTransferMatrix(TransferM, h, J, Jd, Jt, Beta);


	long double Val = MaxEigenValue(TransferM);
	FreeNrg = CalculateFreeNrg(Val, Beta);

	return FreeNrg;
}

double Step(double h , double Beta, double FreeNrgOld, const double h_Step , double J, double  Jd, double Jt)
{
	// returns free energy
	double FreeNrg;

	FreeNrg = FindFreeNrg(h, Beta, h_Step , J, Jd, Jt);

	double M = FindMagetization( FreeNrg, FreeNrgOld, h_Step);
	//std::cout << FreeNrgOld  << ' '<< FreeNrg << "\n";


	std::string Text;
	Text = std::to_string(h) + "  ";
	Text += std::to_string(M) + "\n";
	TextOut(Text, 0);

	return FreeNrg;
}


int main()
{
	const int Steps = 100;
	double h = 0.1, Beta = 10;
	double h_Step = 0.1;
	double FreeNrg, FreeNrgOld;

	double J = -1.0, Jd = -0.0, Jt = -1.0;

	FreeNrgOld = FindFreeNrg(h, Beta, h_Step, J, Jd, Jt);
	
	int i;
	for (i = 0; i < Steps; i++)
	{
		FreeNrg = Step(h, Beta, FreeNrgOld, h_Step, J, Jd, Jt);
		FreeNrgOld = FreeNrg;
		h = h + h_Step;
	}

	std::cin;
	return 0;
}

