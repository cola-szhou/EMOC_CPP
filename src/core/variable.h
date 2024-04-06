#pragma once
#include <map>
#include <variant>
#include <string>
#include <vector>

enum EncodingType
{
    REAL = 1,
    BINARY = 2,
    INTEGER = 3,
    CATEGORICAL = 4,
    PERMUTATION = 5,
    MIXED = 6,
};

namespace emoc
{

    class Variable
    {
    public:
        Variable(std::string name = "");
        virtual ~Variable();
        EncodingType encoding_;
        std::string name_;
    };

    class BinaryVariable : public Variable
    {
    public:
        BinaryVariable(std::string name = "");
        ~BinaryVariable();
        bool Sample();
        bool CheckValue(bool value);

        EncodingType encoding_ = BINARY;
        std::string name_;
    };

    class RealVariable : public Variable
    {
    public:
        RealVariable(double lower_bound, double upper_bound, std::string name = "");
        ~RealVariable();
        double Sample();
        bool CheckValue(double value);

        double lower_bound_;
        double upper_bound_;
        EncodingType encoding_ = REAL;
        std::string name_;
    };

    struct IntegerVariable : public Variable
    {
    public:
        IntegerVariable(int lower_bound, int upper_bound, std::string name = "");
        ~IntegerVariable();
        int Sample();
        bool CheckValue(int value);

        int lower_bound_;
        int upper_bound_;
        EncodingType encoding_ = INTEGER;
        std::string name_;
    };

    struct CategoricalVariable : public Variable
    {
    public:
        CategoricalVariable(std::vector<std::string> categories, std::string name = "");
        ~CategoricalVariable();
        std::string Sample();
        bool CheckValue(std::string value);

        std::vector<std::string> categories_;
        EncodingType encoding_ = CATEGORICAL;
        std::string name_;
    };

    struct PermutationVariable : public Variable
    {
    public:
        PermutationVariable(int size, std::string name = "");
        ~PermutationVariable();
        std::vector<int> Sample();
        bool CheckValue(std::vector<int> value);

        int size_;
        EncodingType encoding_ = PERMUTATION;
        std::string name_;
    };

    using VariantValue = std::variant<bool, int, double, std::string, std::vector<int>>;
    using VariantType = std::variant<BinaryVariable, RealVariable, IntegerVariable, CategoricalVariable, PermutationVariable>;

    class DecisionSpace
    {
    public:
        DecisionSpace();
        DecisionSpace(std::vector<VariantType> space);
        ~DecisionSpace();
        void AddVariable(VariantType variable);
        void RemoveVariable(int index);
        void ModifyVariable(int index, VariantType variable);
        double GetLowerBound(int index);
        double GetUpperBound(int index);
        VariantValue Sample(int index);
        std::vector<std::string> GetCategoricalSpace(int index);

        std::vector<VariantType> space_;
    };
}