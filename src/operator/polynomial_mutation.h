#pragma once
#include <vector>

// #include "core/global.h"
#include "core/individual.h"
#include "core/emoc_utility_structures.h"

namespace emoc {

	void PolynomialMutationPopulation(Individual **pop, int pop_num, std::vector<double>& lower_bound, std::vector<double>& upper_bound, MutationParameter& mutation_para);
	void PolynomialMutationIndividual(Individual *ind, std::vector<double>& lower_bound, std::vector<double>& upper_bound, MutationParameter& mutation_para);

}