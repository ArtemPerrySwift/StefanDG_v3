#include "MaterialPhase.h"
#include <stdexcept>

void MaterialPhase::read(FILE* file, MaterialPhase materialPhaseIt[], uint8_t nMaterialsPhases)
{
	for (uint8_t i = 0; i < nMaterialsPhases; ++i)
	{
		if (feof(file))
		{
			throw std::runtime_error("Failed to read material string");
		}
		else
		{
			uint8_t nameLength = 0;
			int currentChar = fgetc(file);

			while (currentChar != ' ')
			{
				if (currentChar == EOF)
				{
					throw std::runtime_error("Unexpected end of material file reading material name");
				}
				else
				{
					nameLength++;
					currentChar = fgetc(file);
				}
			}

			fseek(file, -nameLength, SEEK_CUR);

			char* materialName = new char[nameLength + 1];
			char state;

			if (fscanf(file, "%s%c%lf%lf%lf", materialName, &state, &materialPhaseIt->thermalConductivity, &materialPhaseIt->heatCapacity, &materialPhaseIt->density) != 4)
			{
				throw std::runtime_error("Faled to read material properties");
			}

			switch (state)
			{
			case 's':
			{
				materialPhaseIt->state = SOLID;
				break;
			}
			case 'l':
			{
				materialPhaseIt->state = LIQUID;
				break;
			}
			default:
			{
				throw std::runtime_error("Wrong agregation state code");
			}
			}


			materialPhaseIt->name = materialName;
			++materialPhaseIt;
		}
	}
}

void MaterialPhase::freeNamesMemory(MaterialPhase materialPhaseIt[], uint8_t nMaterialsPhases)
{
	for (uint8_t i = 0; i < nMaterialsPhases; ++i)
	{
		delete[] materialPhaseIt->name;
		++materialPhaseIt;
	}
}
