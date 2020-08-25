#include "database.h"
//#include <chrono>

int main_test(int argc, char** argv)
{
	array gen_degs = { 1, 2, 4, 8, 11, 12, 16, 18, 18, 21 };
	//array3d rels = { {{0, 3}}, {{1, 2}, {0, 2}} };
	array3d gb = { {{0, 1, 1, 1}},
{{1, 3}, {0, 2, 2, 1}},
{{1, 1, 2, 1}},
{{0, 3, 2, 1}},
{{0, 1, 2, 2}},
{{2, 3}, {1, 2, 3, 1}},
{{2, 1, 3, 1}},
{{0, 1, 4, 1}},
{{2, 1, 4, 1}},
{{1, 2, 4, 1}},
{{1, 1, 3, 2}},
{{3, 1, 4, 1}},
{{3, 1, 5, 1}, {1, 1, 8, 1}, {0, 2, 7, 1}},
{{2, 2, 5, 1}, {0, 2, 8, 1}},
{{1, 1, 7, 1}} };
	//add_rels(gb, rels, gen_degs, -1);
	array3d polys = { {{0, 1}} };
	array4d result = ann_seq(gb, polys, gen_degs, 21);
	std::cout << result << '\n';

	//auto start = std::chrono::system_clock::now();

	//auto end = std::chrono::system_clock::now();
	//std::chrono::duration<double> elapsed = end - start;
	//std::cout << "Elapsed time: " << elapsed.count() << "s\n";
	return 0;
}