#pragma once
#include <iostream>
#include <vector>

namespace emoc {

	double** UniformPoint(int num, int* weight_num, int obj_num);
	std::vector<std::vector<double>> UniformPointWrap(int num, int* weight_num, int obj_num);


}


