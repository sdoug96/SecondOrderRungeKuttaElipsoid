#pragma once

class SecondOrderRungeKutta
{

public:

	SecondOrderRungeKutta();
	~SecondOrderRungeKutta();

	void Solve(float stepSize, int totalSteps, float mass);
};

