#pragma once
#include "MaterialPhase.h"
#include "Volume.h"
#include "Boundary.h"
#include "ConformInterface.h"
#include "NonconformInterface.h"

#include <utility>

struct Model
{
public:
	const MaterialPhase* materialPhases;
	const Volume* volumes;
	const Boundary* volumesBoundaries;
	const Boundary::Condition* conditions;
	const ConformInterface* conformInterfaces;
	const NonconformInterface* nonconformInterfaces;

	unsigned nMaterialPhases;
	unsigned nVolumes;
	unsigned nVolumesBoundaries;
	unsigned nConformInterfaces;
	unsigned nNonconformInterfaces;
	unsigned nSharedBoundaries;
	unsigned nConditions;

	void initilizeByCurrentGMSHModel();
	void clear();

private:
	void determineBoundaries(Boundary*& boundariesMap, unsigned int& nBoundaries);
	void determineVolumesMaterials(std::pair<int, const MaterialPhase*>* const materialPhaseByVolumeTag, const unsigned int nVolumes);
};
