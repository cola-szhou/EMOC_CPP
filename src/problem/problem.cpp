#include "problem/problem.h"

namespace emoc
{

	Problem::Problem(int dec_num, int obj_num) : dec_num_(dec_num), obj_num_(obj_num), problem_name_("Problem"), dec_space_(DecisionSpace())
	{
		// the default encoding type is REAL
		encoding_ = REAL;
	}

	Problem::~Problem()
	{
	}

	void Problem::CalCon(Individual *ind)
	{
		// do nothing
	}

}
