#include "mymath.h"

bool compare(const Poly& poly1, const Poly& poly2) { return deg(poly1) < deg(poly2); }

bool gcd_nonzero(const Mon& m1, const Mon& m2)
{
	size_t len = std::min(m1.size(), m2.size());
	for (size_t i = 0; i < len; i++)
		if (m1[i] || m2[i])
			return true;
	return false;
}

Poly& simplify(const Relations& rels, Poly& poly)
{
	Poly result;
	auto pbegin = poly.begin();
	auto pend = poly.end();
	while (pbegin != pend) {
		bool irreducible = true;
		for (const Poly& rel : rels) {
			if (divides(rel[0], *pbegin)) {
				Mon q = *pbegin / rel[0];
				Poly rel1 = rel * q; // TODO: improve
				Poly tmp;
				std::set_symmetric_difference(pbegin, pend, rel1.begin(), rel1.end(),
					std::back_inserter(tmp), std::greater<Mon>());
				poly = std::move(tmp);
				pbegin = poly.begin();
				pend = poly.end();
				irreducible = false;
				break;
			}
		}
		if (irreducible) {
			result.push_back(std::move(*pbegin));
			pbegin++;
		}
	}
	poly = std::move(result);
	return poly;
}

void add_rel(Relations& rels, const Poly& rel)
{
	if (rel.empty())
		return;
	Relations heap = { rel };
	while (!heap.empty()) {
		std::pop_heap(heap.begin(), heap.end(), compare);
		Poly r = std::move(heap.back());
		heap.pop_back();
		simplify(rels, r);
		if (!r.empty()) {
			for (Poly& r1 : rels) {
				if (gcd_nonzero(r[0], r1[0])) {
					if (divides(r[0], r1[0]))
						r1.clear();
					else {
						Mon mlcm = lcm(r[0], r1[0]);
						Mon q = mlcm / r[0];
						Mon q1 = mlcm / r1[0];
						Poly new_rel = r * q + r1 * q1;
						if (!new_rel.empty()) {
							heap.push_back(std::move(new_rel));
							std::push_heap(heap.begin(), heap.end(), compare);
						}
					}
				}
			}
			rels.erase(std::remove(rels.begin(), rels.end(), Poly()), rels.end());
			rels.push_back(r);
		}
	}
}

void test()
{
	Poly rel1 = { {10}, {7, 3} };
	Poly rel2 = { {7, 4}, {2, 9} };
	Relations rels;
	add_rel(rels, rel1);
	add_rel(rels, rel2);
	for (auto r : rels)
		std::cout << r << std::endl;
}