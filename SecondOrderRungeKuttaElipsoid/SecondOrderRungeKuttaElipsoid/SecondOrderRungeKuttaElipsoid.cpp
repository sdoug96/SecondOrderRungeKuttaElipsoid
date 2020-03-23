#include "pch.h"
#include <iostream>
#include "SecondOrderRungeKutta.h"

int main()
{
	//Define and initialise second-order Runge-Kutta point generator
	SecondOrderRungeKutta* SORK;
	SORK = new SecondOrderRungeKutta;

	SORK->Solve(0.08f, 800, 3.14f);
}