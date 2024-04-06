#include <iostream>
#include <memory>
#include <sstream>

#include "problem/problem.h"
#include "problem/problem_head_collect.h"
#include "core/individual.h"
#include "core/emoc_utility_structures.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
#include "algorithm/algorithm_head_collect.h"
#include "core/py_global.h"
#include "operator/bit_mutation.h"
#include "operator/de.h"
#include "operator/order_crossover.h"
#include "operator/polynomial_mutation.h"
#include "operator/swap_mutation.h"
#include "operator/tournament_selection.h"
#include "operator/uniform_crossover.h"
#include "operator/sbx.h"
#include "core/nd_sort.h"
#include "core/uniform_point.h"
#include "core/utility.h"
#include "metric/metric_head_collect.h"
#include "core/variable.h"
#include "random/random.h"

// using emoc::EMOCManager;
// using emoc::Individual;
// using emoc::NSGA2, emoc::CMOEAD;
// using emoc::Problem;
// using emoc::Py_Global;
// using emoc::ZDT1, emoc::ZDT2;

namespace py = pybind11;

namespace emoc
{
	// for polymorphism of Problem Class in python
	class PyProblem : public emoc::Problem
	{
	public:
		using emoc::Problem::Problem;

		/* Trampoline (need one for each virtual function) */
		void CalObj(emoc::Individual *ind) override
		{
			PYBIND11_OVERRIDE_PURE(
				void,		   /* Return type */
				emoc::Problem, /* Parent class */
				CalObj,		   /* Name of function in C++ (must match Python name) */
				ind			   /* Argument(s) */
			);
		}

		void CalCon(emoc::Individual *ind) override
		{
			PYBIND11_OVERLOAD_PURE(
				void,
				emoc::Problem,
				CalCon,
				ind);
		}
	};

	pybind11::object VariantValueToPyObject(const emoc::VariantValue &value)
	{
		return std::visit([](auto &&arg) -> pybind11::object
						  { return pybind11::cast(arg); },
						  value);
	}

	VariantValue PyObjectToVariantValue(const py::handle &obj)
	{
		if (py::isinstance<py::bool_>(obj))
		{
			return obj.cast<bool>();
		}
		else if (py::isinstance<py::int_>(obj))
		{
			return obj.cast<int>();
		}
		else if (py::isinstance<py::float_>(obj))
		{
			return obj.cast<double>();
		}
		else if (py::isinstance<py::str>(obj))
		{
			return obj.cast<std::string>();
		}
		else if (py::isinstance<py::list>(obj))
		{
			return obj.cast<std::vector<int>>();
		}
		throw std::runtime_error("Type not supported.");
	}

	class ArrayWrapper
	{
	public:
		ArrayWrapper(std::vector<double> &vec) : m_vector(vec)
		{
		}
		py::object get(py::slice slice)
		{
			size_t start, stop, step, slicelength;
			if (!slice.compute(m_vector.size(), &start, &stop, &step, &slicelength))
				throw py::error_already_set();
			std::vector<double> result;
			result.reserve(slicelength);
			for (size_t i = 0; i < slicelength; ++i)
			{
				result.push_back(m_vector[start]);
				start += step;
			}
			return py::cast(result);
		}

		py::object get(size_t index)
		{
			check_index(index);
			return py::cast(m_vector[index]);
		}

		py::object get()
		{
			std::cout << "get all" << std::endl;
			return py::cast(m_vector);
		}

		void set(py::slice slice, const std::vector<double> &value)
		{
			size_t start, stop, step, slicelength;
			if (!slice.compute(m_vector.size(), &start, &stop, &step, &slicelength))
				throw py::error_already_set();
			if (slicelength != value.size())
				throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
			for (size_t i = 0; i < slicelength; ++i)
			{
				m_vector[start + i * step] = value[i];
			}
		}

		void set(size_t index, double value)
		{
			check_index(index);
			m_vector[index] = value;
		}

		void set(const std::vector<double> &value)
		{
			m_vector = value;
		}

		int size()
		{
			return m_vector.size();
		}

		std::string to_string() const
		{
			std::stringstream ss;
			ss << "[";
			for (size_t i = 0; i < m_vector.size(); ++i)
			{
				ss << m_vector[i];
				if (i != m_vector.size() - 1)
				{
					ss << ", ";
				}
			}
			ss << "]";
			return ss.str();
		}

		void append(double value)
		{
			m_vector.push_back(value);
		}

	private:
		std::vector<double> &m_vector;

		void check_index(size_t index)
		{
			if (index >= m_vector.size())
				throw std::out_of_range("Index out of bounds");
		}
	};

	class VariantArrayWrapper
	{
	public:
		VariantArrayWrapper(std::vector<emoc::VariantValue> &vec) : m_vector(vec) {}
		py::object get(py::slice slice)
		{
			size_t start, stop, step, slicelength;
			if (!slice.compute(m_vector.size(), &start, &stop, &step, &slicelength))
				throw py::error_already_set();
			std::vector<emoc::VariantValue> result;
			result.reserve(slicelength);
			for (size_t i = 0; i < slicelength; ++i)
			{
				result.push_back(m_vector[start]);
				start += step;
			}
			return py::cast(result);
		}

		py::object get(size_t index)
		{
			check_index(index);
			return py::cast(m_vector[index]);
		}

		py::object get()
		{
			std::cout << "get all" << std::endl;
			return py::cast(m_vector);
		}

		void set(py::slice slice, const std::vector<emoc::VariantValue> &value)
		{
			size_t start, stop, step, slicelength;
			if (!slice.compute(m_vector.size(), &start, &stop, &step, &slicelength))
				throw py::error_already_set();
			if (slicelength != value.size())
				throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
			for (size_t i = 0; i < slicelength; ++i)
			{
				m_vector[start + i * step] = value[i];
			}
		}

		void set(size_t index, emoc::VariantValue value)
		{
			check_index(index);
			m_vector[index] = value;
		}

		void set(const std::vector<emoc::VariantValue> &value)
		{
			if (value.size() != m_vector.size())
				throw std::runtime_error("Left and right hand size of slice assignment have different sizes!");
			m_vector = value;
		}

		int size()
		{
			return m_vector.size();
		}

		std::string to_string() const
		{
			std::stringstream ss;
			ss << "[";
			for (size_t i = 0; i < m_vector.size(); ++i)
			{
				if (i > 0)
					ss << ", ";
				ss << variantValueToString(m_vector[i]);
			}
			ss << "]";
			return ss.str();
		}

		void append(emoc::VariantValue value)
		{
			m_vector.push_back(value);
		}

	private:
		std::vector<emoc::VariantValue> &m_vector;

		void check_index(size_t index)
		{
			if (index >= m_vector.size())
				throw std::out_of_range("Index out of bounds");
		}

		std::string variantValueToString(const VariantValue &value) const
		{
			return std::visit([](auto &&arg) -> std::string
							  {
            using T = std::decay_t<decltype(arg)>;
            if constexpr (std::is_same_v<T, bool>) {
                return arg ? "true" : "false";
            } else if constexpr (std::is_same_v<T, int> || std::is_same_v<T, double>) {
                return std::to_string(arg);
            } else if constexpr (std::is_same_v<T, std::string>) {
                return arg; // Directly return the string
            } else if constexpr (std::is_same_v<T, std::vector<int>>) {
                std::stringstream ss;
                ss << "[";
                for (size_t i = 0; i < arg.size(); i++) {
                    ss << arg[i];
                    if (i != arg.size() - 1) ss << ", ";
                }
                ss << "]";
                return ss.str();
            } else {
                return "unknown type";
            } },
							  value);
		}
	};
}
void bindOperators(py::module &m)
{
	py::class_<emoc::MutationParameter>(m, "MutationParameter")
		.def(py::init<>())
		.def_readwrite("pro", &emoc::MutationParameter::pro)
		.def_readwrite("index1", &emoc::MutationParameter::index1)
		.def_readwrite("index2", &emoc::MutationParameter::index2);

	py::class_<emoc::CrossoverParameter>(m, "CrossoverParameter")
		.def(py::init<>())
		.def_readwrite("pro", &emoc::CrossoverParameter::pro)
		.def_readwrite("index1", &emoc::CrossoverParameter::index1)
		.def_readwrite("index2", &emoc::CrossoverParameter::index2);

	// define bit flip mutation operator
	m.def("BitFlipMutation", &emoc::BitFlipMutation, "Perform bit flip mutation for an individual");

	// define differential evolution operator
	m.def("DE", &emoc::DE, "Perform differential evolution operation");

	// define order mutation operator
	m.def("OrderCrossover", &emoc::OrderCrossover, "Perform Order Crossover on two individuals");

	// define polynomial mutation operator
	m.def("PolynomialMutationIndividual", &emoc::PolynomialMutationIndividual, "Perform polynomial mutation for an individual");

	// define simulated binary crossover operator
	m.def("SBX", &emoc::SBX, "Perform simulated binary crossover on two individuals");

	// define swap mutation operator
	m.def("SwapMutation", &emoc::SwapMutation, "Perform swap mutation for an individual");

	// define tournament selection operators
	m.def("TournamentByRank", &emoc::TournamentByRank, "Perform tournament selection by rank");

	m.def("TournamentByFitness", &emoc::TournamentByFitness, "Perform tournament selection by fitness");

	m.def("TournamentByCustom", &emoc::TournamentByCustom, "Perform custom tournament selection");

	// define uniform crossover operator
	m.def("UniformCrossover", &emoc::UniformCrossover,
		  "Perform uniform crossover on two individuals");

	// define non-dominated sort
	m.def(
		"NonDominatedSort", [](py::list individuals, int pop_num, int obj_num, bool is_consider_cons)
		{
        // translate py::list to emoc::Individual** and call the NonDominatedSort function
        std::vector<emoc::Individual*> individual_ptrs;
        for (size_t i = 0; i < individuals.size(); i++) {
            emoc::Individual* ind = individuals[i].cast<emoc::Individual*>();
            individual_ptrs.push_back(ind);
        }
        emoc::NonDominatedSort(individual_ptrs.data(), pop_num, obj_num, is_consider_cons); },
		"Non-dominated sort function");

	m.def("UniformPoint", &emoc::UniformPointWrap, "Generate uniform points");

	m.def("CalInverseChebycheff", &emoc::CalInverseChebycheffWrap, "Calculate inverse Chebycheff");
}

void bindAlgorithms(py::module &m)
{
	// define NSGA2 algorithm
	py::class_<emoc::NSGA2>(m, "NSGA2")
		.def(py::init<emoc::Py_Global *, emoc::Problem *, double, double, double, double>())
		.def("Solve", &emoc::NSGA2::Solve)
		.def("PrintResult", &emoc::NSGA2::PrintResult);

	// define CMOEAD algorithm
	py::class_<emoc::CMOEAD>(m, "CMOEAD")
		.def(py::init<emoc::Py_Global *, emoc::Problem *>())
		.def("Solve", &emoc::CMOEAD::Solve)
		.def("PrintResult", &emoc::CMOEAD::PrintResult);
}

void bindProblems(py::module &m)
{
	// define Problem base class, so we can write custom problem derived from this class in python
	py::class_<emoc::Problem, emoc::PyProblem /* <--- trampoline*/>(m, "Problem")
		.def(py::init<int, int>())
		.def("CalObj", &emoc::Problem::CalObj)
		.def("CalCon", &emoc::Problem::CalCon)
		.def_readwrite("dec_num_", &emoc::Problem::dec_num_)
		.def_readwrite("obj_num_", &emoc::Problem::obj_num_)
		.def_readwrite("dec_space_", &emoc::Problem::dec_space_)
		.def_property(
			"lower_bound_",
			[](emoc::Problem &p)
			{
				return emoc::ArrayWrapper(p.lower_bound_);
			},
			[](emoc::Problem &p, const std::vector<double> &value)
			{
				p.lower_bound_ = value;
			})
		.def_property(
			"upper_bound_",
			[](emoc::Problem &p)
			{
				return emoc::ArrayWrapper(p.upper_bound_);
			},
			[](emoc::Problem &p, const std::vector<double> &value)
			{
				p.upper_bound_ = value;
			})
		.def_readwrite("encoding_", &emoc::Problem::encoding_)
		.def_readwrite("problem_name_", &emoc::Problem::problem_name_);

	py::enum_<EncodingType>(m, "EncodingType")
		.value("REAL", REAL)
		.value("BINARY", BINARY)
		.value("INTEGER", INTEGER)
		.value("CATEGORICAL", CATEGORICAL)
		.value("PERMUTATION", PERMUTATION)
		.value("MIXED", MIXED);

	// define ZDT problems
	py::class_<emoc::ZDT1, emoc::Problem>(m, "ZDT1")
		.def(py::init<int, int>())
		.def("cal_obj", &emoc::ZDT1::CalObj);

	py::class_<emoc::ZDT2, emoc::Problem>(m, "ZDT2")
		.def(py::init<int, int>())
		.def("cal_obj", &emoc::ZDT2::CalObj);
}

void bindMetrics(py::module &m)
{
	// define HV metric

	// define IGD metric
	m.def(
		"CalculateIGD", [](py::list individuals, std::vector<std::vector<double>> pf_data_py)
		{
		std::vector<emoc::Individual*> individual_ptrs;
        for (size_t i = 0; i < individuals.size(); i++) {
            emoc::Individual* ind = individuals[i].cast<emoc::Individual*>();
            individual_ptrs.push_back(ind);
        }
		int pf_size = pf_data_py.size();
		int obj_num = pf_data_py[0].size();
		double **pf_data = new double *[pf_size];
		for (int i = 0; i < pf_size; ++i)
		{
			pf_data[i] = new double[obj_num];
			for(int j = 0; j < obj_num; ++j){
				pf_data[i][j] = pf_data_py[i][j];
			}
		}
		double igd = emoc::CalculateIGD(individual_ptrs.data(), individuals.size(), obj_num, pf_data, pf_size);

		for (int i = 0; i < pf_size; ++i)
		{
			delete[] pf_data[i];
		}
		delete[] pf_data;
		return igd; },
		"Calculate IGD metric");

	// defind GD metric
	m.def(
		"CalculateGD", [](py::list individuals, std::vector<std::vector<double>> pf_data_py)
		{
		std::vector<emoc::Individual*> individual_ptrs;
        for (size_t i = 0; i < individuals.size(); i++) {
            emoc::Individual* ind = individuals[i].cast<emoc::Individual*>();
            individual_ptrs.push_back(ind);
        }
		int pf_size = pf_data_py.size();
		int obj_num = pf_data_py[0].size();
		double **pf_data = new double *[pf_size];
		for (int i = 0; i < pf_size; ++i)
		{
			pf_data[i] = new double[obj_num];
			for(int j = 0; j < obj_num; ++j){
				pf_data[i][j] = pf_data_py[i][j];
			}
		}
		double gd = emoc::CalculateGD(individual_ptrs.data(), individuals.size(), obj_num, pf_data, pf_size);
		for (int i = 0; i < pf_size; ++i)
		{
			delete[] pf_data[i];
		}
		delete[] pf_data;
		return gd; },
		"Calculate GD metric");

	// define IGD+ metric
	m.def(
		"CalculateIGDPlus", [](py::list individuals, std::vector<std::vector<double>> pf_data_py)
		{
		std::vector<emoc::Individual*> individual_ptrs;
        for (size_t i = 0; i < individuals.size(); i++) {
            emoc::Individual* ind = individuals[i].cast<emoc::Individual*>();
            individual_ptrs.push_back(ind);
        }
		int pf_size = pf_data_py.size();
		int obj_num = pf_data_py[0].size();
		double **pf_data = new double *[pf_size];
		for (int i = 0; i < pf_size; ++i)
		{
			pf_data[i] = new double[obj_num];
			for(int j = 0; j < obj_num; ++j){
				pf_data[i][j] = pf_data_py[i][j];
			}
		}
		double igdplus = emoc::CalculateIGDPlus(individual_ptrs.data(), individuals.size(), obj_num, pf_data, pf_size);
		for (int i = 0; i < pf_size; ++i)
		{
			delete[] pf_data[i];
		}
		delete[] pf_data;
		return igdplus; },
		"Calculate IGD+ metric");

	// define GD+ metric
	m.def(
		"CalculateGDPlus", [](py::list individuals, std::vector<std::vector<double>> pf_data_py)
		{
		std::vector<emoc::Individual*> individual_ptrs;
        for (size_t i = 0; i < individuals.size(); i++) {
            emoc::Individual* ind = individuals[i].cast<emoc::Individual*>();
            individual_ptrs.push_back(ind);
        }
		int pf_size = pf_data_py.size();
		int obj_num = pf_data_py[0].size();
		double **pf_data = new double *[pf_size];
		for (int i = 0; i < pf_size; ++i)
		{
			pf_data[i] = new double[obj_num];
			for(int j = 0; j < obj_num; ++j){
				pf_data[i][j] = pf_data_py[i][j];
			}
		}
		double gdplus = emoc::CalculateGDPlus(individual_ptrs.data(), individuals.size(), obj_num, pf_data, pf_size);
		for (int i = 0; i < pf_size; ++i)
		{
			delete[] pf_data[i];
		}
		delete[] pf_data;
		return gdplus; },
		"Calculate GDPlus metric");

	// define Spacing metric
	m.def(
		"CalculateSpacing", [](py::list individuals)
		{
		std::vector<emoc::Individual*> individual_ptrs;
		for (size_t i = 0; i < individuals.size(); i++) {
			emoc::Individual* ind = individuals[i].cast<emoc::Individual*>();
			individual_ptrs.push_back(ind);
		}
		double spacing = emoc::CalculateSpacing(individual_ptrs.data(), individuals.size(), individual_ptrs[0]->dec_.size());
		return spacing; },
		"Calculate Spacing metric");

	// define HV metric
}

PYBIND11_MODULE(EMOC, m)
{
	// define Individual class
	py::class_<emoc::Individual>(m, "Individual")
		.def(py::init<int, int>())
		.def_property(
			"dec_", [](emoc::Individual &ind)
			{ return emoc::VariantArrayWrapper(ind.dec_); },
			[](emoc::Individual &ind, const std::vector<emoc::VariantValue> &value)
			{ ind.dec_ = value; })
		.def_property(
			"obj_", [](emoc::Individual &ind)
			{ return emoc::ArrayWrapper(ind.obj_); },
			[](emoc::Individual &ind, const std::vector<double> &value)
			{ ind.obj_ = value; })
		.def_property(
			"con_", [](emoc::Individual &ind)
			{ return emoc::ArrayWrapper(ind.con_); },
			[](emoc::Individual &ind, const std::vector<double> &value)
			{ ind.con_ = value; })
		.def_readwrite("rank_", &emoc::Individual::rank_)
		.def_readwrite("fitness_", &emoc::Individual::fitness_);

	py::class_<emoc::Record>(m, "Record")
		.def(py::init<>())
		.def_readwrite("pop_", &emoc::Record::pop)
		.def_readwrite("runtime_", &emoc::Record::runtime)
		.def_readwrite("iteration_", &emoc::Record::iteration);

	// define Global class to save all parameters
	py::class_<emoc::Py_Global>(m, "Py_Global")
		.def(py::init<>())
		.def("SetParam", &emoc::Py_Global::SetParam)
		.def("Restart", &emoc::Py_Global::Restart)
		.def("ResetPopulationSize", &emoc::Py_Global::ResetPopulationSize)
		.def("RecordPop", &emoc::Py_Global::RecordPop)
		.def_readwrite("parent_population_", &emoc::Py_Global::parent_population_)
		.def_readwrite("offspring_population_", &emoc::Py_Global::offspring_population_)
		.def_readwrite("mixed_population_", &emoc::Py_Global::mixed_population_)
		.def_readwrite("iteration_num_", &emoc::Py_Global::iteration_num_)
		.def_readwrite("current_evaluation_", &emoc::Py_Global::current_evaluation_)
		.def_readwrite("is_customized_init_pop_", &emoc::Py_Global::is_customized_init_pop_)
		.def_readwrite("population_num_", &emoc::Py_Global::population_num_)
		.def_readwrite("dec_num_", &emoc::Py_Global::dec_num_)
		.def_readwrite("obj_num_", &emoc::Py_Global::obj_num_)
		.def_readwrite("max_evaluation_", &emoc::Py_Global::max_evaluation_)
		.def_readwrite("output_interval_", &emoc::Py_Global::output_interval_)
		.def_readwrite("record_", &emoc::Py_Global::record_);
	// .def_property(
	// 	"dec_lower_bound_", [](emoc::Py_Global &g)
	// 	{ return emoc::ArrayWrapper(g.dec_lower_bound_); },
	// 	[](emoc::Py_Global &g, const std::vector<double> &value)
	// 	{
	// 		g.dec_lower_bound_ = value;
	// 	})
	// .def_property(
	// 	"dec_upper_bound_", [](emoc::Py_Global &g)
	// 	{ return emoc::ArrayWrapper(g.dec_upper_bound_); },
	// 	[](emoc::Py_Global &g, const std::vector<double> &value)
	// 	{
	// 		g.dec_upper_bound_ = value;
	// 	});

	py::class_<emoc::ArrayWrapper>(m, "ArrayWrapper")
		.def("__getitem__", [](emoc::ArrayWrapper &a, py::handle h)
			 {
			if (h.is_none()) {
				return a.get();
			}
			else if (py::isinstance<py::slice>(h)) {
				return a.get(h.cast<py::slice>());
			}
			else {
				return a.get(h.cast<size_t>());
			} })
		.def("__setitem__", [](emoc::ArrayWrapper &a, py::slice slice, const std::vector<double> &value)
			 { a.set(slice, value); })
		.def("__setitem__", [](emoc::ArrayWrapper &a, size_t index, double value)
			 { a.set(index, value); })
		.def("__setitem__", [](emoc::ArrayWrapper &a, const std::vector<double> &value)
			 { a.set(value); })
		.def("__len__", &emoc::ArrayWrapper::size)
		.def("__repr__", &emoc::ArrayWrapper::to_string)
		.def("__str__", &emoc::ArrayWrapper::to_string)
		.def("append", &emoc::ArrayWrapper::append);

	py::class_<emoc::VariantArrayWrapper>(m, "VariantArrayWrapper")
		.def("__getitem__", [](emoc::VariantArrayWrapper &a, py::handle h)
			 {
			if (h.is_none()) {
				return a.get();
			}
			else if (py::isinstance<py::slice>(h)) {
				return a.get(h.cast<py::slice>());
			}
			else {
				return a.get(h.cast<size_t>());
			} })
		.def("__setitem__", [](emoc::VariantArrayWrapper &a, py::slice slice, const std::vector<emoc::VariantValue> &value)
			 { a.set(slice, value); })
		.def("__setitem__", [](emoc::VariantArrayWrapper &a, size_t index, emoc::VariantValue value)
			 { a.set(index, value); })
		.def("__setitem__", [](emoc::VariantArrayWrapper &a, const std::vector<emoc::VariantValue> &value)
			 { a.set(value); })
		.def("__len__", &emoc::VariantArrayWrapper::size)
		.def("__repr__", &emoc::VariantArrayWrapper::to_string)
		.def("__str__", &emoc::VariantArrayWrapper::to_string)
		.def("append", &emoc::VariantArrayWrapper::append);

	py::class_<emoc::Variable, std::shared_ptr<emoc::Variable>>(m, "Variable")
		.def(py::init<std::string>())
		.def_readwrite("encoding_", &emoc::Variable::encoding_)
		.def_readwrite("name_", &emoc::Variable::name_);

	py::class_<emoc::BinaryVariable, emoc::Variable, std::shared_ptr<emoc::BinaryVariable>>(m, "BinaryVariable")
		.def(py::init<>())
		.def_readwrite("encoding_", &emoc::BinaryVariable::encoding_)
		.def_readwrite("name_", &emoc::BinaryVariable::name_)
		.def("Sample", &emoc::BinaryVariable::Sample);

	py::class_<emoc::RealVariable, emoc::Variable, std::shared_ptr<emoc::RealVariable>>(m, "RealVariable")
		.def(py::init<double, double, std::string>())
		.def_readwrite("encoding_", &emoc::RealVariable::encoding_)
		.def_readwrite("name_", &emoc::RealVariable::name_)
		.def("Sample", &emoc::RealVariable::Sample);

	py::class_<emoc::IntegerVariable, emoc::Variable, std::shared_ptr<emoc::IntegerVariable>>(m, "IntegerVariable")
		.def(py::init<int, int, std::string>())
		.def_readwrite("encoding_", &emoc::IntegerVariable::encoding_)
		.def_readwrite("name_", &emoc::IntegerVariable::name_)
		.def("Sample", &emoc::IntegerVariable::Sample);

	py::class_<emoc::CategoricalVariable, emoc::Variable, std::shared_ptr<emoc::CategoricalVariable>>(m, "CategoricalVariable")
		.def(py::init<std::vector<std::string>, std::string>())
		.def_readwrite("encoding_", &emoc::CategoricalVariable::encoding_)
		.def_readwrite("name_", &emoc::CategoricalVariable::name_)
		.def("Sample", &emoc::CategoricalVariable::Sample);

	py::class_<emoc::PermutationVariable, emoc::Variable, std::shared_ptr<emoc::PermutationVariable>>(m, "PermutationVariable")
		.def(py::init<int, std::string>())
		.def_readwrite("encoding_", &emoc::PermutationVariable::encoding_)
		.def_readwrite("name_", &emoc::PermutationVariable::name_)
		.def("Sample", &emoc::PermutationVariable::Sample);

	m.def("randomize", &randomize, "Randomize the seed");

	py::class_<emoc::DecisionSpace>(m, "DecisionSpace")
		.def(py::init<>())
		.def(py::init<std::vector<emoc::VariantType>>())
		.def("AddVariable", &emoc::DecisionSpace::AddVariable)
		.def("RemoveVariable", &emoc::DecisionSpace::RemoveVariable)
		.def("ModifyVariable", &emoc::DecisionSpace::ModifyVariable)
		.def("GetLowerBound", &emoc::DecisionSpace::GetLowerBound)
		.def("GetUpperBound", &emoc::DecisionSpace::GetUpperBound)
		.def("Sample", &emoc::DecisionSpace::Sample)
		.def("GetCategoricalSpace", &emoc::DecisionSpace::GetCategoricalSpace)
		.def_readwrite("space_", &emoc::DecisionSpace::space_);

	// *******Operators*******
	bindOperators(m);

	// *******Problem*******
	bindProblems(m);

	// *******Algorithm*******
	bindAlgorithms(m);

	// *******Metric*******
	bindMetrics(m);
}