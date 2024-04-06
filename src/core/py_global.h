#pragma once
#include <string>
#include <vector>

#include "core/individual.h"
#include "problem/problem.h"
#include "algorithm/algorithm.h"
#include <pybind11/pybind11.h>
#include "core/variable.h"

namespace py = pybind11;

namespace emoc
{
	typedef struct
	{
		std::vector<Individual *> pop;
		double runtime;
		int iteration;
	} Record;
	// Global class holds all necessary parameter settings and datas for algorithms to run and
	// provides some useful foundmental functions.
	class Py_Global
	{
	public:
		// operator parameter structures
		typedef struct
		{
			double crossover_pro;
			double eta_c;
		} SBXPara;

		typedef struct
		{
			double crossover_pro;
			double F;
			double K;
		} DEPara;

		typedef struct
		{
			double muatation_pro;
			double eta_m;
		} PolyMutationPara;

	public:
		Py_Global();
		~Py_Global();

		void Init();
		void Start();
		void Restart();
		void ResetPopulationSize(int pop_num);

		void InitializePopulation(Individual **pop, int pop_num, Problem *problem_); // initialize given population, i.e. set the decision variables' value
		void InitializeIndividual(Individual *ind, Problem *problem_);				 // initialize given individual, i.e. set the decision variables' value

		// for python dlls
		void SetCustomInitialPop(std::vector<std::vector<double>> &initial_pop);
		void SetParam(int dec_num, int obj_num, DecisionSpace dec_space, int population_num, int output_interval, int max_evaluation);
		void RecordPop(int real_popnum, double runtime);

	public:
		int dec_num_;
		int obj_num_;
		int population_num_;
		int max_evaluation_;
		int output_interval_ = 1;
		bool record_X_ = false;

		std::vector<Individual *> parent_population_;
		std::vector<Individual *> offspring_population_;
		std::vector<Individual *> mixed_population_;
		std::vector<Record *> record_;

		// std::vector<double> dec_lower_bound_; // set by problem's lower bound
		// std::vector<double> dec_upper_bound_; // set by problem's lower bound
		DecisionSpace dec_space_; // set by problem's decision variable space

		SBXPara sbx_parameter_;
		DEPara de_parameter_;
		PolyMutationPara pm_parameter_;

		int iteration_num_;
		int current_evaluation_;

		int run_id_;
		int thread_id_;

		// for python dlls
		bool is_customized_init_pop_;
		bool is_customized_problem_;

	private:
		void SetDecBound();
		void InitializeProblem();
		void InitializeAlgorithm();
		void AllocateMemory();
		void DestroyMemory();
	};
}