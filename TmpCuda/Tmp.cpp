#include "groebner.h"
#include "benchmark.h"
#include "myio.h"
#include <iostream>

int benchmark1()
{
	array gen_degs = { 1, 1, 1, 1 }; /* x, y, z, w */
	Poly rel0 = { {{0, 4}} }; /* x^4 */
	Poly rel1 = { {{0, 2}, {1, 4}} }; /* x^2y^4 */
	Poly rel2 = { {{0, 1}, {2, 2}}, {{0, 2}, {1, 1}} }; /* x^2y + xz^2 */
	Poly rel3 = { {{3, 4}}, {{2, 4}} }; /* z^2 + w^4 */
	Poly1d gb;
	grbn::add_rels(gb, { rel0, rel1, rel2, rel3 }, gen_degs, -1);

	Poly2d ann = grbn::ann_seq_v2(gb, { Poly{ {{0, 1}} } }, gen_degs, -1);
	
	std::cout << gb << '\n';
	std::cout << ann << '\n';
	return 0;
}

int main()
{
	return benchmark1();
}