#pragma once

#include <vector>
#include <core/variable.h>

namespace emoc
{

	// Class for each solution in optimization process, and
	// it holds some basic datas of each solution.
	class Individual
	{
	public:
		Individual(int dec_num, int obj_num);
		~Individual();

	public:
		int rank_;
		double fitness_;

		// std::vector<double> dec_;
		std::vector<VariantValue> dec_;
		std::vector<double> obj_;
		std::vector<double> con_;
	};

}