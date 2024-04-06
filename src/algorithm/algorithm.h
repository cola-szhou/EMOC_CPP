#pragma once
#include <ctime>

#include "problem/problem.h"
#include "core/individual.h"
#include "core/emoc_utility_structures.h"
#include "core/py_global.h"

namespace emoc
{

	class Py_Global;
	// Basic class of all algorithms, it provides some useful functions to derived classes.
	// All derived classes need to override Run() which implement the real algorithm details
	class Algorithm
	{
	public:
		Algorithm(Py_Global *global, Problem *problem);
		virtual ~Algorithm();

		// void PrintPop();
		void PrintResult();
		virtual void Solve() = 0;

		inline int GetRealPopNum() { return real_popnum_; }
		inline double GetRuntime() { return runtime_; }

	protected:
		bool IsTermination();
		void EvaluatePop(Individual **pop, int pop_num);
		void EvaluateInd(Individual *ind);

		// utility functions for population operation
		int MergePopulation(Individual **pop_src1, int pop_num1, Individual **pop_src2, int pop_num2, Individual **pop_dest);
		void CopyIndividual(Individual *ind_src, Individual *ind_dest);
		void SwapIndividual(Individual *ind1, Individual *ind2);

		void TrackPopulation();
		// void PlotPopulation(Individual** pop, int gen);
		// bool CheckStopAndPause();

	public:
		Py_Global *g_GlobalSettings; // pointer to current run's global settings
		int thread_id_;				 // current thread id
		int real_popnum_;			 // real population number
		double runtime_;			 // total time for the algorithm to optimize the problem
		CrossoverParameter cross_para_;
		MutationParameter mutation_para_;
		Problem *problem_;

		clock_t start_, end_;
	};

}