#include "core/variable.h"
#include "random/random.h"
#include <iostream>

namespace emoc
{
    Variable::Variable(std::string name) : name_(name) {}
    Variable::~Variable() {}

    BinaryVariable::BinaryVariable(std::string name) : name_(name) {}
    BinaryVariable::~BinaryVariable() {}
    bool BinaryVariable::Sample()
    {
        return rnd(0, 2) != 0;
    }
    bool BinaryVariable::CheckValue(bool value) { return true; }

    RealVariable::RealVariable(double lower_bound, double upper_bound, std::string name) : lower_bound_(lower_bound), upper_bound_(upper_bound), name_(name) {}

    RealVariable::~RealVariable() {}

    double RealVariable::Sample()
    {
        return rndreal(lower_bound_, upper_bound_);
    }

    bool RealVariable::CheckValue(double value)
    {
        if (value >= lower_bound_ && value <= upper_bound_)
            return true;
        return false;
    }

    IntegerVariable::IntegerVariable(int lower_bound, int upper_bound, std::string name) : lower_bound_(lower_bound), upper_bound_(upper_bound), name_(name) {}

    IntegerVariable::~IntegerVariable() {}

    int IntegerVariable::Sample()
    {
        return rnd(lower_bound_, upper_bound_);
    }

    bool IntegerVariable::CheckValue(int value)
    {
        if (value >= lower_bound_ && value <= upper_bound_)
            return true;
        return false;
    }

    CategoricalVariable::CategoricalVariable(std::vector<std::string> categories, std::string name) : categories_(categories), name_(name) {}

    CategoricalVariable::~CategoricalVariable() {}

    std::string CategoricalVariable::Sample()
    {
        return categories_[rnd(0, categories_.size())];
    }

    bool CategoricalVariable::CheckValue(std::string value)
    {
        for (const auto &category : categories_)
        {
            if (category == value)
                return true;
        }
        return false;
    }

    PermutationVariable::PermutationVariable(int size, std::string name) : size_(size), name_(name) {}

    PermutationVariable::~PermutationVariable() {}

    std::vector<int> PermutationVariable::Sample()
    {
        std::vector<int> perm(size_);
        random_permutation(perm.data(), size_);
        return perm;
    }

    bool PermutationVariable::CheckValue(std::vector<int> value)
    {
        if (value.size() != size_)
            return false;

        std::vector<int> flag(size_, 0);
        for (const auto &val : value)
        {
            if (val < 0 || val >= size_)
                return false;
            if (flag[val] == 1)
                return false;
            flag[val] = 1;
        }
        return true;
    }

    DecisionSpace::DecisionSpace() {}
    DecisionSpace::DecisionSpace(std::vector<VariantType> space) : space_(space) {}
    DecisionSpace::~DecisionSpace() {}
    void DecisionSpace::AddVariable(VariantType variable) { space_.push_back(variable); }
    void DecisionSpace::RemoveVariable(int index)
    {
        if (index >= 0 && index < space_.size())
            space_.erase(space_.begin() + index);
        else
            throw std::out_of_range("Index out of DecisionSpace range");
    }
    void DecisionSpace::ModifyVariable(int index, VariantType variable)
    {
        if (index >= 0 && index < space_.size())
            space_[index] = variable;
        else
            throw std::out_of_range("Index out of DecisionSpace range");
    }
    double DecisionSpace::GetLowerBound(int index)
    {
        if (index < 0 || index >= space_.size())
            throw std::out_of_range("Index out of DecisionSpace range");
        return std::visit([](auto &arg) -> double
                          { using T = std::decay_t<decltype(arg)>;
                          if constexpr(std::is_same_v<T, RealVariable> || std::is_same_v<T, IntegerVariable>)
                            return (double)(arg.lower_bound_);
                          else
                            throw std::invalid_argument("Variable type does not have lower bound"); },
                          space_[index]);
    }
    double DecisionSpace::GetUpperBound(int index)
    {
        if (index < 0 || index >= space_.size())
            throw std::out_of_range("Index out of DecisionSpace range");
        return std::visit([](auto &arg) -> double
                          { using T = std::decay_t<decltype(arg)>;
                          if constexpr(std::is_same_v<T, RealVariable> || std::is_same_v<T, IntegerVariable>)
                            return (double)(arg.upper_bound_);
                          else
                            throw std::invalid_argument("Variable type does not have upper bound"); },
                          space_[index]);
    }
    std::vector<std::string> DecisionSpace::GetCategoricalSpace(int index)
    {
        if (index < 0 || index >= space_.size())
            throw std::out_of_range("Index out of DecisionSpace range");
        return std::visit([](auto &arg) -> std::vector<std::string>
                          { using T = std::decay_t<decltype(arg)>;
                          if constexpr(std::is_same_v<T, CategoricalVariable>)
                            return arg.categories_;
                          else
                            throw std::invalid_argument("Variable type is not categorical"); },
                          space_[index]);
    }

    VariantValue DecisionSpace::Sample(int index)
    {
        if (index < 0 || index >= space_.size())
            throw std::out_of_range("Index out of DecisionSpace range");
        return std::visit([](auto &arg) -> VariantValue
                          { return arg.Sample(); },
                          space_[index]);
    }

}