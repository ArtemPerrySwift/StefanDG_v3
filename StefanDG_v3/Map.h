#pragma once
#include <stdexcept>

namespace Map
{
	template<typename TagType, class ElementType>
	struct MapElement
	{
		TagType tag;
		ElementType elem;
		MapElement* next;

		MapElement(const TagType tag, const ElementType elem) : tag{ tag }, next{ nullptr }, elem(elem)
		{
		}

		MapElement() : tag{ 0 }, next{nullptr}, elem()
		{
		}

		~MapElement()
		{
			delete next;
		}
	};

	template<typename TagType, class ElementType>
	void addElement(const TagType tag, const ElementType element, MapElement<TagType, ElementType> map[], const TagType n)
	{
		const TagType index = tag % n;
		MapElement<TagType, ElementType>*  mapElement = map + index;

		if (mapElement->tag == 0)
		{
			mapElement->tag = tag;
			mapElement->elem = element;

			return;
		}
		
		while (mapElement->next != nullptr)
		{
			mapElement = mapElement->next;
		}

		mapElement->next = new MapElement<TagType, ElementType>(tag, element);
	}

	template<typename TagType, class ElementType>
	const MapElement<TagType, ElementType>* tryElementAdding(const TagType tag, const ElementType element, MapElement<TagType, ElementType> map[], const TagType n)
	{
		const TagType index = tag % n;
		MapElement<TagType, ElementType>* mapElement = map + index;

		if (mapElement->tag == 0)
		{
			mapElement->tag = tag;
			mapElement->elem = element;
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

		mapElement->next = new MapElement<TagType, ElementType>(tag, element);
		return mapElement->next;
	}

	template<typename TagType, class ElementType>
	ElementType& getElement(const TagType tag, MapElement<TagType, ElementType> map[], const TagType n)
	{
		const TagType index = tag % n;
		MapElement<TagType, ElementType>* mapElement = map + index;

		if (mapElement->tag != tag)
		{
			while (mapElement->next != nullptr)
			{
				if (mapElement->tag == tag)
				{
					return mapElement->elem;
				}

				mapElement = mapElement->next;
			}

			if (mapElement->tag != tag)
			{
				throw std::runtime_error("There is no element in map with such tag");
			}
		}
		return mapElement->elem;
	}


}


