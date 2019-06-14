#include "mymath.h"

bool gcd_nonzero(const Mon& m1, const Mon& m2)
{
	size_t len = std::min(m1.size(), m2.size());
	for (size_t i = 0; i < len; ++i)
		if (m1[i] && m2[i])
			return true;
	return false;
}

Poly& simplify(Poly& poly, const Relations& rels)
{
	Poly result;
	auto pbegin = poly.begin(); auto pend = poly.end();
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
				pbegin = poly.begin(); pend = poly.end();
				irreducible = false;
				break;
			}
		}
		if (irreducible) {
			result.push_back(std::move(*pbegin));
			++pbegin;
		}
	}
	poly = std::move(result);
	return poly;
}

struct rel_heap_t {
	Poly poly;
	int deg;
};

inline bool compare(const rel_heap_t& s1, const rel_heap_t& s2) { return s1.deg > s2.deg; }

void add_rel(const Poly& rel, Relations& rels, std::vector<int>* pgen_degs)
{
	if (rel.empty())
		return;
	std::vector<rel_heap_t> heap = { rel_heap_t{rel, pgen_degs == NULL ? 
		deg(rel) : deg(rel, *pgen_degs)} };
	while (!heap.empty()) {
		std::pop_heap(heap.begin(), heap.end(), compare);
		rel_heap_t heap_ele = std::move(heap.back());
		heap.pop_back();

		Poly& r = heap_ele.poly;
		simplify(r, rels);
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

						heap.push_back(rel_heap_t{ new_rel, pgen_degs == NULL ? 
							deg(mlcm) : deg(mlcm, *pgen_degs) });
						std::push_heap(heap.begin(), heap.end(), compare);
					}
				}
			}
			rels.erase(std::remove(rels.begin(), rels.end(), Poly()), rels.end());
			rels.push_back(std::move(r));
		}
	}
}

void test()
{
	Mon m1 = { 1, 0, 1 };
	Mon m2 = { 0, 1, 0 };
	std::cout << gcd_nonzero(m1, m2) << std::endl;
}