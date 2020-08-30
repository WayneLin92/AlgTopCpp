#include "database.h"
#include <chrono>

int main_test(int argc, char** argv)
{
	auto start = std::chrono::system_clock::now();

	array gen_degs;
	int n_max = 8;
	for (int d = 1; d <= n_max; d++) {
		for (int i = 0; i <= n_max - d; i++) {
			int j = i + d;
			gen_degs.push_back((1 << j) - (1 << i));
		}
	}
	array3d rels;
	for (int d = 2; d <= n_max; d++) {
		for (int i = 0; i <= n_max - d; i++) {
			int j = i + d;
			array2d rel;
			for (int k = i + 1; k < j; k++) {
				int a = (1 << k) - (1 << i);
				int b = (1 << j) - (1 << k);
				auto p1 = std::find(gen_degs.begin(), gen_degs.end(), a);
				auto p2 = std::find(gen_degs.begin(), gen_degs.end(), b);
				int index1 = int(p1 - gen_degs.begin());
				int index2 = int(p2 - gen_degs.begin());
				rel = add(rel, index1 < index2 ? array2d{ { index1, 1, index2, 1 } } : array2d{ { index2, 1, index1, 1 } });
			}
			rels.push_back(std::move(rel));
		}
	}

	array3d gb;
	add_rels(gb, rels, gen_degs, -1);

	/*for (size_t i = 0; i < gb.size(); ++i)
		std::cout << i << ':' << gb[i].size() << '\n';*/

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Elapsed time: " << elapsed.count() << "s\n";
	return 0;
}