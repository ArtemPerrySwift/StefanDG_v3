#pragma once
#include <stdio.h>
#include <cstdint>

struct MaterialPhase
{
	static enum State : uint8_t { SOLID, LIQUID };

	const char* name;
	State state;
	double thermalConductivity;
	double heatCapacity;
	double density;


	static void read(FILE* file, MaterialPhase materialsPhases[], uint8_t nMaterialsPhases);
	static void freeNamesMemory(MaterialPhase materialsPhases[], uint8_t nMaterialsPhases);
};

