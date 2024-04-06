#include "operator/bit_mutation.h"

#include <cmath>
#include "random/random.h"

// BitFlipMutation is for binary encoding

namespace emoc
{

	void BitFlipMutation(Individual *ind, MutationParameter &mutation_para)
	{

		for (auto &value : ind->dec_)
		{
			std::visit([&](auto &&dec)
					   {
				using T = std::decay_t<decltype(dec)>;
				if constexpr(std::is_same_v<T, bool>){
					if (randomperc() < mutation_para.pro)
						dec = !dec;
				} },
					   value);
		}
	}

}