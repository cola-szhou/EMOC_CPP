#include "problem/zdt/zdt.h"

#include <cmath>

#include "core/macro.h"
#include "core/variable.h"
// #include "core/global.h"

namespace emoc
{

	ZDT1::ZDT1(int dec_num, int obj_num) : Problem(dec_num, obj_num)
	{
		for (int i = 0; i < dec_num; ++i)
		{
			dec_space_.AddVariable(RealVariable(0.0, 1.0));
		}
		problem_name_ = "ZDT1";
	}

	ZDT1::~ZDT1()
	{
	}

	void ZDT1::CalObj(Individual *ind)
	{
		double f1 = 0, f2 = 0;
		double g = 0, h = 0;

		f1 = std::get<double>(ind->dec_.at(0));
		for (int i = 1; i < dec_num_; i++)
		{
			g += std::get<double>(ind->dec_.at(i));
		}

		g = 1 + g * 9.0 / (dec_num_ - 1);
		h = 1 - sqrt((double)(f1 / g));
		f2 = g * h;

		ind->obj_[0] = f1;
		ind->obj_[1] = f2;
	}

	ZDT2::ZDT2(int dec_num, int obj_num) : Problem(dec_num, obj_num)
	{
		for (int i = 0; i < dec_num; ++i)
		{
			dec_space_.AddVariable(RealVariable(0.0, 1.0));
		}

		problem_name_ = "ZDT2";
	}

	ZDT2::~ZDT2()
	{
	}

	void ZDT2::CalObj(Individual *ind)
	{
		double f1 = 0, f2 = 0;
		double g = 0, h = 0;

		f1 = std::get<double>(ind->dec_.at(0));
		for (int i = 1; i < dec_num_; i++)
		{
			g += std::get<double>(ind->dec_.at(i));
		}

		g = 1 + g * 9.0 / (dec_num_ - 1);
		h = 1 - (f1 / g) * (f1 / g);
		f2 = g * h;

		ind->obj_[0] = f1;
		ind->obj_[1] = f2;
	}

	ZDT3::ZDT3(int dec_num, int obj_num) : Problem(dec_num, obj_num)
	{
		for (int i = 0; i < dec_num; ++i)
		{
			dec_space_.AddVariable(RealVariable(0.0, 1.0));
		}
		problem_name_ = "ZDT3";
	}

	ZDT3::~ZDT3()
	{
	}

	void ZDT3::CalObj(Individual *ind)
	{
		double f1 = 0, f2 = 0;
		double g = 0, h = 0;

		f1 = std::get<double>(ind->dec_.at(0));
		for (int i = 1; i < dec_num_; i++)
		{
			g += std::get<double>(ind->dec_.at(i));
		}

		g = 1 + g * 9.0 / (dec_num_ - 1);
		h = 1.0 - sqrt(f1 / g) - (f1 / g) * sin(10.0 * PI * f1);
		f2 = g * h;

		ind->obj_[0] = f1;
		ind->obj_[1] = f2;
	}

	ZDT4::ZDT4(int dec_num, int obj_num) : Problem(dec_num, obj_num)
	{
		dec_space_.AddVariable(RealVariable(0.0, 1.0));
		for (int i = 1; i < dec_num; i++)
		{
			dec_space_.AddVariable(RealVariable(-5.0, 5.0));
		}
		problem_name_ = "ZDT4";
	}

	ZDT4::~ZDT4()
	{
	}

	void ZDT4::CalObj(Individual *ind)
	{
		double f1 = 0, f2 = 0;
		double g = 0, h = 0;

		f1 = std::get<double>(ind->dec_.at(0));
		for (int i = 1; i < dec_num_; i++)
		{
			double tmp = std::get<double>(ind->dec_.at(i));
			g += tmp * tmp - 10.0 * cos(4.0 * PI * tmp);
		}

		g += 10.0 * (dec_num_ - 1) + 1.0;
		h = 1.0 - sqrt(f1 / g);
		f2 = g * h;

		ind->obj_[0] = f1;
		ind->obj_[1] = f2;
	}

	ZDT6::ZDT6(int dec_num, int obj_num) : Problem(dec_num, obj_num)
	{
		for (int i = 0; i < dec_num; ++i)
		{
			dec_space_.AddVariable(RealVariable(0.0, 1.0));
		}
	}

	ZDT6::~ZDT6()
	{
	}

	void ZDT6::CalObj(Individual *ind)
	{
		double f1 = 0, f2 = 0;
		double g = 0, h = 0;

		f1 = 1.0 - exp(-4.0 * std::get<double>(ind->dec_.at(0))) * pow(sin(6.0 * PI * std::get<double>(ind->dec_.at(0))), 6.0);
		for (int i = 1; i < dec_num_; i++)
		{
			g += std::get<double>(ind->dec_.at(i));
		}

		g = 9.0 * pow(g / (dec_num_ - 1), 0.25) + 1.0;
		h = 1.0 - (f1 / g) * (f1 / g);
		f2 = g * h;

		ind->obj_[0] = f1;
		ind->obj_[1] = f2;
	}

}