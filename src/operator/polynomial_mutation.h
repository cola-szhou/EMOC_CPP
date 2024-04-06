#pragma once
#include <vector>

// #include "core/global.h"
#include "core/individual.h"
#include "core/emoc_utility_structures.h"

namespace emoc
{

	void PolynomialMutationPopulation(Individual **pop, int pop_num, DecisionSpace dec_space, MutationParameter &mutation_para);
	void PolynomialMutationIndividual(Individual *ind, DecisionSpace dec_space, MutationParameter &mutation_para);

}