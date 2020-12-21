#include "groebner.h"
#include "benchmark.h"
#include <iostream>

int benchmark1()
{
	array gen_degs;
	int n_max = 8;
	for (int d = 1; d <= n_max; d++) {
		for (int i = 0; i <= n_max - d; i++) {
			int j = i + d;
			gen_degs.push_back((1 << j) - (1 << i));
		}
	}
	Poly1d rels;
	for (int d = 2; d <= n_max; d++) {
		for (int i = 0; i <= n_max - d; i++) {
			int j = i + d;
			Poly rel;
			for (int k = i + 1; k < j; k++) {
				int a = (1 << k) - (1 << i);
				int b = (1 << j) - (1 << k);
				auto p1 = std::find(gen_degs.begin(), gen_degs.end(), a);
				auto p2 = std::find(gen_degs.begin(), gen_degs.end(), b);
				int index1 = int(p1 - gen_degs.begin());
				int index2 = int(p2 - gen_degs.begin());
				rel = add(rel, index1 < index2 ? Poly{ {{index1, 1}, {index2, 1}} } : Poly{ {{index2, 1}, {index1, 1}} });
			}
			rels.push_back(std::move(rel));
		}
	}

	Poly1d gb1;
	Poly1d gb2;
	std::sort(rels.begin(), rels.end(), [&gen_degs](const Poly& p1, const Poly& p2) {
		return get_deg(p1, gen_degs) < get_deg(p2, gen_degs); });

	Timer timer;
	grbn::add_rels(gb2, rels, gen_degs, -1);
	timer.print("add_rels: ");
	grbn::AddRels(gb1, rels, gen_degs, -1);
	timer.print("AddRels: ");
	gb1.clear();
	gb2.clear();
	grbn::add_rels(gb2, rels, gen_degs, -1);
	timer.print("add_rels: ");
	grbn::AddRels(gb1, rels, gen_degs, -1);
	timer.print("AddRels: ");
	gb1.clear();
	gb2.clear();
	grbn::add_rels(gb2, rels, gen_degs, -1);
	timer.print("add_rels: ");
	grbn::AddRels(gb1, rels, gen_degs, -1);
	timer.print("AddRels: ");


	std::cout << gb1.size() << "==" << gb2.size() << '\n';
	return 0;
}

int main()
{
	return benchmark1();
}