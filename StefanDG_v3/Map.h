#pragma once
#include <stdexcept>

namespace Map
{
	template<typename TagType, class ValueType>
	struct Element
	{
		TagType tag;
		ValueType value;
		Element* next;

		Element(const TagType tag, const ValueType value) : tag{ tag }, next{ nullptr }, value(value)
		{
		}

		Element() : tag{ 0 }, next{nullptr}, value()
		{
		}

		~Element()
		{
			delete next;
		}
	};

	template<typename TagType, class ValueType>
	void addElement(const TagType tag, const ValueType element, Element<TagType, ValueType> map[], const TagType n)
	{
		const TagType index = tag % n;
		Element<TagType, ValueType>*  mapElement = map + index;

		if (mapElement->tag == 0)
		{
			mapElement->tag = tag;
			mapElement->value = element;

			return;
		}
		
		while (mapElement->next != nullptr)
		{
			mapElement = mapElement->next;
		}

		mapElement->next = new Element<TagType, ValueType>(tag, element);
	}

	template<typename TagType, class ValueType>
	const Element<TagType, ValueType>* tryElementAdding(const TagType tag, const ValueType element, Element<TagType, ValueType> map[], const TagType n)
	{
		const TagType index = tag % n;
		Element<TagType, ValueType>* mapElement = map + index;

		if (mapElement->tag == 0)
		{
			mapElement->tag = tag;
			mapElement->value = element;
			return mapElement;
		}
		
		while (mapElement->next != nullptr)
		{
			if (mapElement->tag == tag)
			{
				return mapElement;
			}
			mapElement = mapElement->next;
		}

		mapElement->next = new Element<TagType, ValueType>(tag, element);
		return mapElement->next;
	}

	template<typename TagType, class ValueType>
	ValueType& getElement(const TagType tag, Element<TagType, ValueType> map[], const TagType n)
	{
		const TagType index = tag % n;
		Element<TagType, ValueType>* mapElement = map + index;

		if (mapElement->tag != tag)
		{
			while (mapElement->next != nullptr)
			{
				if (mapElement->tag == tag)
				{
					return mapElement->value;
				}

				mapElement = mapElement->next;
			}

			if (mapElement->tag != tag)
			{
				throw std::runtime_error("There is no element in map with such tag");
			}
		}
		return mapElement->value;
	}

	template<typename TagType, class ValueType>
	ValueType& getExistedValue(const TagType tag, const Element<TagType, ValueType>* element)
	{
		while (element->tag != tag)
		{
			element = element->next;
		}
		return element->value;
	}
}


