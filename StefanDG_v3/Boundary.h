#pragma once
#include <cstdint>

struct Boundary
{
	struct Condition
	{
		static enum class Type { DIRICHLET, NEWMAN, STEFAN, CONFORM_INTERFACE, NONCONFORM_INTERFACE, HOMOGENEOUS_NEWMAN, UNDEFINED };

		Type type;
		unsigned int index;
	};
	//static enum class ConditionType { DIRICHLET, NEWMAN, STEFAN, CONFORM_INTERFACE, NONCONFORM_INTERAFACE, HOMOGENEOUS_NEWMAN, UNDEFINED };
	
	int tag;
	bool isPlane;
	Condition condition;
	//ConditionType condition;
};