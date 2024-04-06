#include "operator/order_crossover.h"

#include <vector>
#include <cmath>

#include "random/random.h"

// OrderCrossover is for permutation encoding
namespace emoc
{

	void OrderCrossover(Individual *parent1, Individual *parent2, Individual *offspring1, Individual *offspring2)
	{
		auto &p1_vec = std::get<std::vector<int>>(parent1->dec_[0]);
		auto &p2_vec = std::get<std::vector<int>>(parent2->dec_[0]);

		int dec_num = p1_vec.size();
		std::vector<int> ht1(dec_num, 0), ht2(dec_num, 0);
		int k = rnd(0, dec_num - 1);

		offspring1->dec_[0] = std::vector<int>(dec_num);
		offspring2->dec_[0] = std::vector<int>(dec_num);
		auto &o1_vec = std::get<std::vector<int>>(offspring1->dec_[0]);
		auto &o2_vec = std::get<std::vector<int>>(offspring2->dec_[0]);

		// Copy the first k elements from parent1 and parent2 to offspring1 and offspring2
		for (int i = 0; i <= k; ++i)
		{
			o1_vec[i] = p1_vec[i];
			o2_vec[i] = p2_vec[i];
			ht1[p1_vec[i]] = 1;
			ht2[p2_vec[i]] = 1;
		}

		// Complete the offspring with elements from the other parent that haven't been added yet
		int index1 = k + 1, index2 = k + 1;
		for (int i = 0; i < dec_num; ++i)
		{
			if (ht1[p2_vec[i]] == 0)
			{
				o1_vec[index1++] = p2_vec[i];
				ht1[p2_vec[i]] = 1;
			}
			if (ht2[p1_vec[i]] == 0)
			{
				o2_vec[index2++] = p1_vec[i];
				ht2[p1_vec[i]] = 1;
			}
		}
	}
}
