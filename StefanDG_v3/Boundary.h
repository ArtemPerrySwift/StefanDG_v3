#pragma once
#include <cstdint>
#include "Coordinates.h"

struct Boundary
{
	struct ICondition
	{
		static enum class Type { DIRICHLET_VALUE, DIRICHLET_FUNCTION, NEWMAN_FUNCTION, NEWMAN_VALUE, STEFAN, CONFORM_INTERFACE, NONCONFORM_INTERFACE, HOMOGENEOUS_NEWMAN, UNDEFINED };
		const Type type;

	protected:
		ICondition(Type type) : type{ type }
		{

		}
	};

	template<ICondition::Type conditionType> struct Condition;

	template<>
	struct Condition<ICondition::Type::DIRICHLET_VALUE> : ICondition
	{
		Condition() : ICondition(Type::DIRICHLET_VALUE), value{ 0.0 }
		{

		}
		double value;
	};

	template<>
	struct Condition<ICondition::Type::DIRICHLET_FUNCTION> : ICondition
	{
		Condition() : ICondition(Type::DIRICHLET_FUNCTION)
		{

		}

		static void getValues(const Coordinates points[], uint8_t nPoints, double values[])
		{
			for (uint8_t i = 0; i < nPoints; ++i)
			{
				*values = 0.0;

				++points;
				++values;
			}
		}
	};

	template<>
	struct Condition<ICondition::Type::NEWMAN_VALUE> : ICondition
	{
		Condition() : ICondition(Type::NEWMAN_VALUE), value{ 0.0 }
		{

		}
		double value;
	};

	template<>
	struct Condition<ICondition::Type::NEWMAN_FUNCTION> : ICondition
	{
		Condition() : ICondition(Type::NEWMAN_FUNCTION)
		{

		}

		static void getValues(const Coordinates points[], uint8_t nPoints, double values[])
		{
			for (uint8_t i = 0; i < nPoints; ++i)
			{
				*values = 0.0;

				++points;
				++values;
			}
		}
	};

	template<>
	struct Condition<ICondition::Type::STEFAN> : ICondition
	{
		Condition() : ICondition(Type::STEFAN), value{ 0.0 }
		{

		}
		double value;
	};

	template<>
	struct Condition<ICondition::Type::CONFORM_INTERFACE> : ICondition
	{
		Condition() : ICondition(Type::CONFORM_INTERFACE), index{ 0 }
		{

		}
		unsigned int index;
	};

	template<>
	struct Condition<ICondition::Type::NONCONFORM_INTERFACE> : ICondition
	{
		Condition() : ICondition(Type::NONCONFORM_INTERFACE), index{ 0 }
		{

		}
		unsigned int index;
	};

	template<>
	struct Condition<ICondition::Type::HOMOGENEOUS_NEWMAN> : ICondition
	{
		Condition() : ICondition(Type::HOMOGENEOUS_NEWMAN)
		{

		}
	};

	template<>
	struct Condition<ICondition::Type::STEFAN> : ICondition
	{
		Condition() : ICondition(ICondition::Type::STEFAN), value{0.0}
		{

		}
		double value;
	};

	struct ConditionsSet
	{
		const Condition<ICondition::Type::DIRICHLET_VALUE>* valuesDirichlet;
		unsigned int nValuesDirichlet;
		const Condition<ICondition::Type::NEWMAN_VALUE>* valuesNewman;
		unsigned int nValuesNewman;
		const Condition<ICondition::Type::NONCONFORM_INTERFACE>* nonconformInterfaces;
		unsigned int nNonconformInterfaces;

		Condition<ICondition::Type::STEFAN> stefan;

		static const Condition<ICondition::Type::HOMOGENEOUS_NEWMAN> homogeneousNewmanCondition;
		static const Condition<ICondition::Type::NEWMAN_FUNCTION> newmanFunctionCondition;
		static const Condition<ICondition::Type::DIRICHLET_FUNCTION> dirichletFunctionCondition;
	};

	int tag;
	bool isPlane;
	const ICondition* condition;
};
