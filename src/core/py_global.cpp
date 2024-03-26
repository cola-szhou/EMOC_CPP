#include "core/py_global.h"

#include <iostream>
#include <cstdlib>

#include "random/random.h"
#include <typeinfo>
#include <pybind11/pybind11.h>
#include "problem/problem.h"

namespace py = pybind11;

namespace emoc
{

	Py_Global::Py_Global()
	{
	}

	Py_Global::~Py_Global()
	{
		DestroyMemory();
	}

	void Py_Global::InitializePopulation(Individual **pop, int pop_num, Problem *problem_)
	{
		// When the initial population is set before algorithm starts, return directly.
		if (iteration_num_ == 0 && is_customized_init_pop_)
			return;

		for (int i = 0; i < pop_num; ++i)
		{
			InitializeIndividual(pop[i], problem_);
		}
	}

	void Py_Global::InitializeIndividual(Individual *ind, Problem *problem_)
	{
		if (problem_->encoding_ == Problem::REAL)
		{
			for (int i = 0; i < dec_num_; ++i)
			{
				ind->dec_[i] = rndreal(dec_lower_bound_[i], dec_upper_bound_[i]);
			}
		}
		else if (problem_->encoding_ == Problem::BINARY)
		{
			for (int i = 0; i < dec_num_; ++i)
				ind->dec_[i] = rnd(0, 1);
		}
		else if (problem_->encoding_ == Problem::INTEGER)
		{
			for (int i = 0; i < dec_num_; ++i)
				ind->dec_[i] = rnd(dec_lower_bound_[i], dec_upper_bound_[i]);
		}
		else if (problem_->encoding_ == Problem::PERMUTATION)
		{
			std::vector<int> perm(dec_num_);
			random_permutation(perm.data(), perm.size());
			for (int i = 0; i < dec_num_; ++i)
				ind->dec_[i] = perm[i];
		}
	}

	void Py_Global::Restart()
	{
		iteration_num_ = 0;
		current_evaluation_ = 0;
	}

	void Py_Global::ResetPopulationSize(int pop_num)
	{
		if (population_num_ < pop_num)
		{
			for (int i = population_num_; i < pop_num; ++i)
			{
				parent_population_.push_back(new Individual(dec_num_, obj_num_));
				offspring_population_.push_back(new Individual(dec_num_, obj_num_));
				mixed_population_.push_back(new Individual(dec_num_, obj_num_));
				mixed_population_.push_back(new Individual(dec_num_, obj_num_));
			}
		}
		else if (population_num_ > pop_num)
		{
			for (int i = pop_num; i < population_num_; ++i)
			{
				delete parent_population_[i];
				delete offspring_population_[i];
				delete mixed_population_[i];
				delete mixed_population_[i + pop_num];
				parent_population_[i] = nullptr;
				offspring_population_[i] = nullptr;
				mixed_population_[i] = nullptr;
				mixed_population_[i + pop_num] = nullptr;
			}
		}
		population_num_ = pop_num;
	}

	void Py_Global::SetParam(int dec_num, int obj_num, std::vector<double> lower_bound, std::vector<double> upper_bound, int population_num, int output_interval, int max_evaluation)
	{
		dec_num_ = dec_num;
		obj_num_ = obj_num;
		population_num_ = population_num;
		iteration_num_ = 0;
		current_evaluation_ = 0;
		is_customized_init_pop_ = false;
		output_interval_ = output_interval;
		max_evaluation_ = max_evaluation;
		dec_lower_bound_ = std::vector<double>(dec_num_);
		dec_upper_bound_ = std::vector<double>(dec_num_);
		for (int i = 0; i < dec_num_; i++)
		{
			dec_lower_bound_[i] = lower_bound[i];
			dec_upper_bound_[i] = upper_bound[i];
		}

		// reserve population space
		parent_population_.reserve(population_num_);
		offspring_population_.reserve(population_num_);
		mixed_population_.reserve(population_num_ * 2);

		// set operator parameters
		sbx_parameter_.crossover_pro = 1.0;
		sbx_parameter_.eta_c = 20.0;
		de_parameter_.crossover_pro = 1.0;
		de_parameter_.F = 0.5;
		de_parameter_.K = 0.5;
		pm_parameter_.muatation_pro = 1.0 / (double)dec_num_;
		pm_parameter_.eta_m = 20.0;

		// allocate memory for all population
		AllocateMemory();
		randomize();
	}

	void Py_Global::SetCustomInitialPop(std::vector<std::vector<double>> &initial_pop)
	{
		int initial_pop_num = initial_pop.size();
		int initial_dec_dim = initial_pop[0].size();
		if (initial_pop_num > population_num_)
			throw std::runtime_error("initial population number is larger than the setted parameter!\n");
		if (initial_dec_dim != dec_num_)
			throw std::runtime_error("initial population decision dimensions is not equal to the setted parameter!\n");

		for (int i = 0; i < initial_pop_num; i++)
			for (int j = 0; j < initial_dec_dim; j++)
				parent_population_[i]->dec_[j] = initial_pop[i][j];

		is_customized_init_pop_ = true;
	}

	void Py_Global::RecordPop(int real_popnum, double runtime)
	{
		Record *rec_ = new Record();
		for (int i = 0; i < real_popnum; i++)
		{
			Individual *ind = new Individual(dec_num_, obj_num_);
			for (int j = 0; j < dec_num_; j++)
				ind->dec_[j] = parent_population_[i]->dec_[j];
			for (int j = 0; j < obj_num_; j++)
				ind->obj_[j] = parent_population_[i]->obj_[j];
			ind->rank_ = parent_population_[i]->rank_;
			ind->fitness_ = parent_population_[i]->fitness_;
			rec_->pop.push_back(ind);
		}
		rec_->runtime = runtime;
		rec_->iteration = iteration_num_;
		record_.push_back(rec_);
	}

	void Py_Global::AllocateMemory()
	{
		for (int i = 0; i < population_num_; ++i)
		{
			parent_population_.push_back(new Individual(dec_num_, obj_num_));
			offspring_population_.push_back(new Individual(dec_num_, obj_num_));
			mixed_population_.push_back(new Individual(dec_num_, obj_num_));
			mixed_population_.push_back(new Individual(dec_num_, obj_num_));
		}
	}

	void Py_Global::DestroyMemory()
	{
		int size = parent_population_.size();
		for (int i = 0; i < size; ++i)
		{
			delete parent_population_[i];
			delete offspring_population_[i];
			delete mixed_population_[i];
			delete mixed_population_[i + size];
			parent_population_[i] = nullptr;
			offspring_population_[i] = nullptr;
			mixed_population_[i] = nullptr;
			mixed_population_[i + size] = nullptr;
		}
		size = record_.size();
		for (int i = 0; i < size; ++i)
		{
			delete record_[i];
			record_[i] = nullptr;
		}
	}

}