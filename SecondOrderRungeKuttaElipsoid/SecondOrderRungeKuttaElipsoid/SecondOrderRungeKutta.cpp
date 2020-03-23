#include "pch.h"
#include "SecondOrderRungeKutta.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <fstream>

using namespace std;

SecondOrderRungeKutta::SecondOrderRungeKutta()
{
}


SecondOrderRungeKutta::~SecondOrderRungeKutta()
{
}

void SecondOrderRungeKutta::Solve(float stepSize, int totalSteps, float mass)
{
	//Define arrays
	float wx[800], wy[800], wz[800];

	//Define initial array values
	wx[0] = 1.0f;
	wy[0] = 1.0f;
	wz[0] = 1.0f;

	//Define axis, not sure where these values come from for cone?
	float a = 3.0f, b = 2.0f, c = 1.0f;

	//Define moments of inertia (will be different for cone)
	float I1 = mass * (pow(b, 2) + pow(c, 2)) / 5.0f;
	float I2 = mass * (pow(a, 2) + pow(c, 2)) / 5.0f;
	float I3 = mass * (pow(a, 2) + pow(b, 2)) / 5.0f;

	//Define appropriate gamma values
	float g1 = (I3 - I2) / I1;
	float g2 = (I1 - I3) / I2;
	float g3 = (I2 - I1) / I3;

	//Define and open text file to be written to
	ofstream SORKOutputData;

	SORKOutputData.open("SORKOutputData.txt");

	//Solve Eulers equations using fourth-order Runge-Kutta algorithm (second-order at the moment)
	for (int i = 0; i < totalSteps - 1; i++)
	{
		//Calculate k1 values
		float kx1 = -stepSize * g1 * wy[i] * wz[i];
		float ky1 = -stepSize * g2 * wx[i] * wz[i];
		float kz1 = -stepSize * g3 * wx[i] * wy[i];

		//Calculate k2 values
		float kx2 = -stepSize * g1 * (wy[i] + 0.5 * ky1) * (wz[i] + 0.5 * kz1);
		float ky2 = -stepSize * g2 * (wx[i] + 0.5 * kx1) * (wz[i] + 0.5 * kz1);
		float kz2 = -stepSize * g3 * (wx[i] + 0.5 * kx1) * (wy[i] + 0.5 * ky1);

		//Calculate next array values (will be different for fourth order)
		wx[i + 1] = wx[i] + kx2;
		wy[i + 1] = wy[i] + ky2;
		wz[i + 1] = wz[i] + kz2;

		//Round values to two decimal places
		wx[i] = roundf(wx[i] * 100.0f) / 100.0f;
		wy[i] = roundf(wy[i] * 100.0f) / 100.0f;
		wz[i] = roundf(wz[i] * 100.0f) / 100.0f;

		//Write data to structured text file
		SORKOutputData << wx[i] << "       " << wy[i] << "       " << wz[i] << endl;
	}

	SORKOutputData.close();

	cout << "Arrays populated successfully!" << endl << endl;
}