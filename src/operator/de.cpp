#include "operator/de.h"

#include "random/random.h"

namespace emoc
{

	void DE(Individual *parent1, Individual *parent2, Individual *parent3, Individual *offspring,
			DecisionSpace dec_space, CrossoverParameter &cross_para)
	{
		double value = 0.0;
		int dec_num = parent1->dec_.size();
		double F = cross_para.index1;

		for (int i = 0; i < dec_num; ++i)
		{
			if (std::holds_alternative<double>(parent1->dec_[i]))
			{
				double yl = dec_space.GetLowerBound(i);
				double yu = dec_space.GetUpperBound(i);

				if (randomperc() < cross_para.pro)
				{
					value = std::get<double>(parent1->dec_.at(i)) + F * (std::get<double>(parent2->dec_.at(i)) - std::get<double>(parent3->dec_.at(i)));
					value = (value > yu) ? yu : (value < yl) ? yl
															 : value;
				}
				else
				{
					value = std::get<double>(parent1->dec_.at(i));
				}
				offspring->dec_[i] = value;
			}
			else
			{
				// if the variable is not real, copy the value from parent1. Is it reasonable?
				offspring->dec_[i] = parent1->dec_[i];
			}
		}
	}

}