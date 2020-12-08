#include "groebner.h"
#include <algorithm>
#include <iterator>

/******** Groebner Basis ********/

Poly pow(const Poly& poly, int n, const Poly1d& gb)
{
	Poly result = { {} };
	if (n == 0)
		return result;
	Poly power = poly;
	while (n) {
		if (n & 1)
			result = reduce(mul(result, power), gb);
		n >>= 1;
		if (n)
			power = reduce(mul(power, power), gb);
	}
	return result;
}

Poly reduce(Poly poly, const Poly1d& gb)
{
	Poly result;
	auto pbegin = poly.begin(); auto pend = poly.end();
	while (pbegin != pend) {
		auto pGb = gb.begin();
		for (; pGb != gb.end(); ++pGb)
			if (divides(pGb->front(), *pbegin))
				break;
		if (pGb == gb.end())
			result.push_back(std::move(*pbegin++));
		else {
			Mon q = div(*pbegin, pGb->front());
			Poly rel1 = mul(*pGb, q);
			Poly poly1;
			std::set_symmetric_difference(pbegin, pend, rel1.begin(), rel1.end(),
				std::back_inserter(poly1));
			poly = std::move(poly1);
			pbegin = poly.begin(); pend = poly.end();
		}
	}
	return result;
}

bool gcd_nonzero(const Mon& mon1, const Mon& mon2)
{
	MonInd k = mon1.begin(), l = mon2.begin();
	while (k != mon1.end() && l != mon2.end()) {
		if (k->gen < l->gen)
			k++;
		else if (k->gen > l->gen)
			l++;
		else
			return true;
	}
	return false;
}

/* Comsume relations from `heap` that is at most in degree `deg` while adding new relations to `heap` that is at most in degree `deg_max`.
** deg=-1 or deg_max=-1 means infinity.
*/
template <typename Fn>
void add_rels_from_heap(Poly1d& gb, RelHeap& heap, Fn _get_deg, int deg, int deg_max)
{
	while (!heap.empty() && (deg == -1 ? (deg_max == -1 || heap.top().t <= deg_max) : heap.top().t <= deg)) {
		PolyWithT heap_ele = MoveFromTop(heap);
		Poly rel = reduce(heap_ele.poly, gb);
		if (!rel.empty()) {
			for (Poly& g : gb) {
				if (gcd_nonzero(rel[0], g[0])) {
					if (divides(rel[0], g[0])) {
						Poly new_rel = rel * div(g[0], rel[0]) + g;
						if (!new_rel.empty())
							heap.push(PolyWithT{ std::move(new_rel), _get_deg(g[0]) });
						g.clear();
					}
					else {
						Mon mlcm = lcm(rel[0], g[0]);
						int deg_new_rel = _get_deg(mlcm);
						if (deg_max == -1 || deg_new_rel <= deg_max) {
							Poly new_rel = rel * div(mlcm, rel[0]) + g * div(mlcm, g[0]);
							if (!new_rel.empty())
								heap.push(PolyWithT{ std::move(new_rel), deg_new_rel });
						}
					}
				}
			}
			RemoveEmptyElements(gb);
			gb.push_back(std::move(rel));
		}
	}
}

template <typename Fn>
void add_rels(Poly1d& gb, const Poly1d& rels, Fn _get_deg, int deg_max)
{
	RelHeap heap;
	for (const Poly& rel : rels)
		if (!rel.empty())
			heap.push(PolyWithT{ rel, _get_deg(rel[0]) });  // rel is copied by value
	add_rels_from_heap(gb, heap, _get_deg, -1, deg_max);
}

void add_rels_from_heap(Poly1d& gb, RelHeap& heap, const array& gen_degs, int t, int deg_max)
{
	add_rels_from_heap(gb, heap, [&gen_degs](const Mon& m) {return get_deg(m, gen_degs); }, t, deg_max);
}

void add_rels_from_heap(Poly1d& gb, RelHeap& heap, const array& gen_degs, const array& gen_degs1, int t, int deg_max)
{
	add_rels_from_heap(gb, heap, FnGetDegV2{ gen_degs, gen_degs1 }, t, deg_max);
}

void add_rels(Poly1d& gb, const Poly1d& rels, const array& gen_degs, int deg_max)
{
	add_rels(gb, rels, FnGetDeg{ gen_degs }, deg_max);
}

void add_rels(Poly1d& gb, const Poly1d& rels, const array& gen_degs, const array& gen_degs1, int deg_max)
{
	add_rels(gb, rels, FnGetDegV2{ gen_degs, gen_degs1 }, deg_max);
}

template <typename Fn>
void add_rels_freemodule(Poly1d& gb, RelHeap& heap, Fn get_deg, int deg, int deg_max)
{
	while (!heap.empty() && (deg == -1 ? (deg_max == -1 || heap.top().t <= deg_max) : heap.top().t <= deg)) {
		PolyWithT heap_ele = MoveFromTop(heap);
		Poly rel = reduce(heap_ele.poly, gb);
		if (!rel.empty()) {
			for (Poly& g : gb) {
				if ((rel[0][0].gen >= 0 || g[0][0].gen >= 0 || rel[0][0].gen == g[0][0].gen) && gcd_nonzero(rel[0], g[0])) {
					if (divides(rel[0], g[0])) {
						Poly new_rel = rel * div(g[0], rel[0]) + g;
						if (!new_rel.empty())
							heap.push(PolyWithT{ std::move(new_rel), get_deg(g[0]) });
						g.clear();
					}
					else {
						Mon mlcm = lcm(rel[0], g[0]);
						int deg_new_rel = get_deg(mlcm);
						if (deg_max == -1 || deg_new_rel <= deg_max) {
							Poly new_rel = rel * div(mlcm, rel[0]) + g * div(mlcm, g[0]);
							if (!new_rel.empty())
								heap.push(PolyWithT{ std::move(new_rel), deg_new_rel });
						}
					}
				}
			}
			RemoveEmptyElements(gb);
			gb.push_back(std::move(rel));
		}
	}
}

/* Compute the generating set of `vectors` inplace */
Poly2d& indecomposables(const Poly1d& gb, Poly2d& vectors, const array& gen_degs, const array& basis_degs)
{
	if (vectors.empty())
		return vectors;
	Poly1d gb1 = gb;
	FnGetDegV2 get_deg_{ gen_degs, basis_degs };

	/* Convert each vector v into a relation \\sum vi x_{-i-1} */
	Poly1d rels;
	array degs;
	for (const Poly1d& v : vectors) {
		Poly rel;
		for (int i = 0; i < basis_degs.size(); ++i)
			if (!v[i].empty())
				rel += v[i] * Mon{ {-i - 1, 1} };
		degs.push_back(get_deg_(rel[0]));
		rels.push_back(std::move(rel));
	}
	array indices = range((int)vectors.size());
	std::sort(indices.begin(), indices.end(), [&degs](int i, int j) {return degs[i] < degs[j]; });

	/* Add relations ordered by degree to gb1 */
	RelHeap heap;
	int deg_max = degs[indices.back()];
	for (int i : indices) {
		add_rels_freemodule(gb1, heap, get_deg_, degs[i], deg_max);
		Poly rel = reduce(rels[i], gb1);
		if (!rel.empty())
			heap.push(PolyWithT{ std::move(rel), degs[i] });
		else
			vectors[i].clear();
	}

	/* Keep only the indecomposables in `vectors` */
	RemoveEmptyElements(vectors);
	return vectors;
}

/* Compute the generating set of linear relations among `polys` */
Poly2d ann_seq(const Poly1d& gb, const Poly1d& polys, const array& gen_degs, int deg_max)
{
	Poly2d result;
	if (polys.empty())
		return result;
	Poly1d rels;
	array gen_degs1;
	int N = (int)polys.size();

	/* Add relations Xi=polys[i] to gb to obtain gb1 */
	for (int i = 0; i < N; ++i) {
		Poly p = polys[i];
		gen_degs1.push_back(get_deg(p, gen_degs));
		p.push_back({ {-i - 1, 1} });
		rels.push_back(std::move(p));
	}
	Poly1d gb1 = gb;
	add_rels(gb1, rels, gen_degs, gen_degs1, deg_max);

	/* Extract linear relations from gb1 */
	for (const Poly& g : gb1) {
		if (g[0][0].gen < 0) {
			Poly1d ann;
			ann.resize(N);
			for (const Mon& m : g) {
				MonInd p = m.begin();
				for (; p != m.end() && p->gen < 0; ++p);
				Mon m1(m.begin(), p), m2(p, m.end());
				ann[size_t(-m1[0].gen) - 1] += reduce(mul(evaluate({ div(m1, { {m1[0].gen, 1} }) }, [&polys](int i) {return polys[size_t(-i) - 1]; }, gb), m2), gb);
			}
			result.push_back(std::move(ann));
		}
	}

	/* Add commutators to linear relations */
	for (int i = 0; i < N; ++i) {
		for (int j = i + 1; j < N; ++j) {
			if (gen_degs1[i] + gen_degs1[j] <= deg_max) {
				Poly1d result_i;
				result_i.resize(N);
				result_i[i] = polys[j];
				result_i[j] = polys[i];
				result.push_back(std::move(result_i));
			}
		}
	}

	indecomposables(gb, result, gen_degs, gen_degs1);
	return result;
}

RelHeap GenerateHeap(const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t, int t_max)
{
	RelHeap heap;
	for (auto pg1 = gb.begin(); pg1 != gb.end(); ++pg1) {
		for (auto pg2 = pg1 + 1; pg2 != gb.end(); ++pg2) {
			if (gcd_nonzero(pg1->front(), pg2->front())) {
				Mon mlcm = lcm(pg1->front(), pg2->front());
				int deg_new_rel = get_deg(mlcm, gen_degs, gen_degs1);
				if (t <= deg_new_rel && deg_new_rel <= t_max) {
					Poly new_rel = (*pg1) * div(mlcm, pg1->front()) + (*pg2) * div(mlcm, pg2->front());
					heap.push(PolyWithT{ std::move(new_rel), deg_new_rel });
				}
			}
		}
	}
	return heap;
}