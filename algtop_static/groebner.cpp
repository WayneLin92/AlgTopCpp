#include "groebner.h"
//#include "myio.h"////
#include <algorithm>
#include <iterator>

#ifdef GROEBNER_MULTITHREAD
#include <future>
constexpr int NUM_THREADS = 5;
#endif

namespace grbn {

/**********************************************************
* Small functions
**********************************************************/

Poly pow(const Poly& poly, int n, const Poly1d& gb)
{
	Poly result = { {} };
	if (n == 0)
		return result;
	Poly power = poly;
	while (n) {
		if (n & 1)
			result = Reduce(mul(result, power), gb);
		n >>= 1;
		if (n)
			power = Reduce(mul(power, power), gb);
	}
	return result;
}

Poly Reduce(Poly poly, const Poly1d& gb)
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

#ifdef GROEBNER_MULTITHREAD
void BatchReduce(const Poly1d& gb, const Poly1d& rels_in, Poly1d& rels_out, int i_start)
{
	for (int i = i_start; i < rels_in.size(); i += NUM_THREADS)
		rels_out[i] = Reduce(rels_in[i], gb);
}
#endif

#ifdef GROEBNER_MULTITHREAD
void BatchReduceV2(const Poly1d& gb, const std::vector<std::unique_ptr<BaseEleBuffer>>& rels_in, Poly1d& rels_out, int i_start)
{
	for (int i = i_start; i < rels_in.size(); i += NUM_THREADS)
		rels_out[i] = Reduce(rels_in[i]->GetPoly(gb), gb);
}
#endif

/**********************************************************
* Groebner basis
**********************************************************/
/* Comsume relations from `buffer` in degree <= `deg`
while adding new relations to `heap` in degree <= `deg_max`.
* deg=-1 or deg_max=-1 means infinity */
template <typename Fn>
void AddRelsB(Poly1d& gb, GbBuffer& buffer, Fn _get_deg, int deg, int deg_max)
{
	auto p_buffer = buffer.begin();
	for (; p_buffer != buffer.end() && (deg == -1 || p_buffer->first <= deg); ++p_buffer) {
		/* Reduce relations from buffer in degree `p_buffer->first` */
#ifndef GROEBNER_MULTITHREAD
		Poly1d rels;
		for (auto& poly : p_buffer->second) {
			Poly rel = Reduce(std::move(poly), gb);
			for (Poly& rel1 : rels)
				if (std::binary_search(rel.begin(), rel.end(), rel1[0]))
					rel += rel1;
			if (!rel.empty())
				rels.push_back(std::move(rel));
		}
#else
		Poly1d rels_tmp(p_buffer->second.size());
		std::vector<std::future<void>> futures;
		for (int i = 0; i < NUM_THREADS; ++i)
			futures.push_back(std::async(std::launch::async, BatchReduce, std::ref(gb), std::ref(p_buffer->second), std::ref(rels_tmp), i));
		for (int i = 0; i < NUM_THREADS; ++i)
			futures[i].wait();
		Poly1d rels;
		for (auto& rel : rels_tmp) {
			for (Poly& rel1 : rels)
				if (std::binary_search(rel.begin(), rel.end(), rel1[0]))
					rel += rel1;
			if (!rel.empty())
				rels.push_back(std::move(rel));
		}
#endif

		/* Add these relations */
		for (auto& rel : rels) {
			if (!rel.empty()) {
				for (Poly& g : gb) {
					if (gcd_nonzero(rel[0], g[0])) {
						Mon lcm = LCM(rel[0], g[0]);
						int deg_new_rel = _get_deg(lcm);
						if (deg_max == -1 || deg_new_rel <= deg_max) {
							Poly new_rel = rel * div(lcm, rel[0]) + g * div(lcm, g[0]);
							if (!new_rel.empty())
								buffer[deg_new_rel].push_back(std::move(new_rel));
						}
					}
				}
				gb.push_back(std::move(rel));
			}
		}
	}
	buffer.erase(buffer.begin(), p_buffer);
}
void AddRelsB(Poly1d& gb, GbBuffer& buffer, const array& gen_degs, int t, int deg_max) { AddRelsB(gb, buffer, FnGetDeg{ gen_degs }, t, deg_max); }
void AddRelsB(Poly1d& gb, GbBuffer& buffer, const array& gen_degs, const array& gen_degs1, int t, int deg_max) { AddRelsB(gb, buffer, FnGetDegV2{ gen_degs, gen_degs1 }, t, deg_max); }

template <typename Fn>
void AddRels(Poly1d& gb, Poly1d rels, Fn _get_deg, int deg_max)
{
	GbBuffer buffer;
	for (Poly& rel : rels) {
		if (!rel.empty()) {
			int deg = _get_deg(rel[0]);
			buffer[deg].push_back(std::move(rel));
		}
	}
	AddRelsB(gb, buffer, _get_deg, -1, deg_max);
}
void AddRels(Poly1d& gb, Poly1d rels, const array& gen_degs, int deg_max) { AddRels(gb, std::move(rels), FnGetDeg{ gen_degs }, deg_max); }
void AddRels(Poly1d& gb, Poly1d rels, const array& gen_degs, const array& gen_degs1, int deg_max) { AddRels(gb, std::move(rels), FnGetDegV2{ gen_degs, gen_degs1 }, deg_max); }

void AddRelsM(Poly1d& gb, GbBuffer& buffer, const array& gen_degs, const array& gen_degs1, int deg, int deg_max)
{
	auto p_buffer = buffer.begin();
	for (; p_buffer != buffer.end() && (deg == -1 || p_buffer->first <= deg); ++p_buffer) {
		/* Reduce relations from buffer in degree `p_buffer->first` */
#ifndef GROEBNER_MULTITHREAD
		Poly1d rels;
		for (auto& poly : p_buffer->second) {
			Poly rel = Reduce(std::move(poly), gb);
			for (Poly& rel1 : rels)
				if (std::binary_search(rel.begin(), rel.end(), rel1[0]))
					rel += rel1;
			if (!rel.empty())
				rels.push_back(std::move(rel));
		}
#else
		Poly1d rels_tmp(p_buffer->second.size());
		std::vector<std::future<void>> futures;
		for (int i = 0; i < NUM_THREADS; ++i)
			futures.push_back(std::async(std::launch::async, BatchReduce, std::ref(gb), std::ref(p_buffer->second), std::ref(rels_tmp), i));
		for (int i = 0; i < NUM_THREADS; ++i)
			futures[i].wait();
		Poly1d rels;
		for (auto& rel : rels_tmp) {
			for (Poly& rel1 : rels)
				if (std::binary_search(rel.begin(), rel.end(), rel1[0]))
					rel += rel1;
			if (!rel.empty())
				rels.push_back(std::move(rel));
		}
#endif

		/* Add these relations */
		for (auto& rel : rels) {
			if (!rel.empty()) {
				for (Poly& g : gb) {
					if (gcd_nonzero(rel[0], g[0]) && !(rel[0][0].gen < 0 && g[0][0].gen < 0 && rel[0][0].gen != g[0][0].gen)) { /* Ignore x_ix_j when i, j<0 */
						Mon lcm = LCM(rel[0], g[0]);
						int deg_new_rel = get_deg(lcm, gen_degs, gen_degs1);
						if (deg_max == -1 || deg_new_rel <= deg_max) {
							Poly new_rel = rel * div(lcm, rel[0]) + g * div(lcm, g[0]);
							if (!new_rel.empty())
								buffer[deg_new_rel].push_back(std::move(new_rel));
						}
					}
				}
				gb.push_back(std::move(rel));
			}
		}
	}
	buffer.erase(buffer.begin(), p_buffer);
}

GbBuffer GenerateBuffer(const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t, int t_max)
{
	GbBuffer buffer;
	for (auto pg1 = gb.begin(); pg1 != gb.end(); ++pg1) {
		for (auto pg2 = pg1 + 1; pg2 != gb.end(); ++pg2) {
			if (gcd_nonzero(pg1->front(), pg2->front())) {
				Mon lcm = LCM(pg1->front(), pg2->front());
				int deg_new_rel = get_deg(lcm, gen_degs, gen_degs1);
				if (t <= deg_new_rel && deg_new_rel <= t_max) {
					Poly new_rel = (*pg1) * div(lcm, pg1->front()) + (*pg2) * div(lcm, pg2->front());
					buffer[deg_new_rel].push_back(std::move(new_rel));
				}
			}
		}
	}
	return buffer;
}

/**********************************************************
* Groebner basis version 2
* This implementation uses polymorphism to reduce the use of memory
* with little cost.
**********************************************************/
template <typename Fn>
void AddRelsBV2(Poly1d& gb, GbBufferV2& buffer, Fn _get_deg, int deg, int deg_max)
{
	auto p_buffer = buffer.begin();
	for (; p_buffer != buffer.end() && (deg == -1 || p_buffer->first <= deg); ++p_buffer) {
		/* Reduce relations from buffer in degree `p_buffer->first` */
#ifndef GROEBNER_MULTITHREAD
		Poly1d rels;
		for (auto& poly : p_buffer->second) {
			Poly rel = Reduce(std::move(poly->GetPoly(gb)), gb);
			for (Poly& rel1 : rels)
				if (std::binary_search(rel.begin(), rel.end(), rel1[0]))
					rel += rel1;
			if (!rel.empty())
				rels.push_back(std::move(rel));
		}
#else
		Poly1d rels_tmp(p_buffer->second.size());
		std::vector<std::future<void>> futures;
		for (int i = 0; i < NUM_THREADS; ++i)
			futures.push_back(std::async(std::launch::async, BatchReduceV2, std::ref(gb), std::ref(p_buffer->second), std::ref(rels_tmp), i));
		for (int i = 0; i < NUM_THREADS; ++i)
			futures[i].wait();
		Poly1d rels;
		for (auto& rel : rels_tmp) {
			for (Poly& rel1 : rels)
				if (std::binary_search(rel.begin(), rel.end(), rel1[0]))
					rel += rel1;
			if (!rel.empty())
				rels.push_back(std::move(rel));
		}
#endif

		/* Add these relations */
		for (auto& rel : rels) {
			if (!rel.empty()) {
				for (int i = 0; i < (int)gb.size(); ++i) {
					if (gcd_nonzero(rel[0], gb[i][0])) {
						Mon gcd = GCD(rel[0], gb[i][0]);
						int deg_new_rel = p_buffer->first + _get_deg(gb[i][0]) - _get_deg(gcd);
						if (deg_max == -1 || deg_new_rel <= deg_max)
							buffer[deg_new_rel].push_back(std::make_unique<GcdEleBuffer>(std::move(gcd), i, (int)gb.size()));
					}
				}
				gb.push_back(std::move(rel));
			}
		}
	}
	buffer.erase(buffer.begin(), p_buffer);
}
void AddRelsBV2(Poly1d& gb, GbBufferV2& buffer, const array& gen_degs, int t, int deg_max) { AddRelsBV2(gb, buffer, FnGetDeg{ gen_degs }, t, deg_max); }
void AddRelsBV2(Poly1d& gb, GbBufferV2& buffer, const array& gen_degs, const array& gen_degs1, int t, int deg_max) { AddRelsBV2(gb, buffer, FnGetDegV2{ gen_degs, gen_degs1 }, t, deg_max); }

template <typename Fn>
void AddRelsV2(Poly1d& gb, Poly1d rels, Fn _get_deg, int deg_max)
{
	GbBufferV2 buffer;
	for (Poly& rel : rels) {
		if (!rel.empty()) {
			int deg = _get_deg(rel[0]);
			buffer[deg].push_back(std::make_unique<PolyEleBuffer>(std::move(rel)));
		}
	}
	AddRelsBV2(gb, buffer, _get_deg, -1, deg_max);
}
void AddRelsV2(Poly1d& gb, Poly1d rels, const array& gen_degs, int deg_max) { AddRelsV2(gb, std::move(rels), FnGetDeg{ gen_degs }, deg_max); }
void AddRelsV2(Poly1d& gb, Poly1d rels, const array& gen_degs, const array& gen_degs1, int deg_max) { AddRelsV2(gb, std::move(rels), FnGetDegV2{ gen_degs, gen_degs1 }, deg_max); }

void AddRelsMV2(Poly1d& gb, GbBufferV2& buffer, const array& gen_degs, const array& gen_degs1, int deg, int deg_max)
{
	auto p_buffer = buffer.begin();
	for (; p_buffer != buffer.end() && (deg == -1 || p_buffer->first <= deg); ++p_buffer) {
		/* Reduce relations from buffer in degree `p_buffer->first` */
#ifndef GROEBNER_MULTITHREAD
		Poly1d rels;
		for (auto& poly : p_buffer->second) {
			Poly rel = Reduce(std::move(poly->GetPoly(gb)), gb);
			for (Poly& rel1 : rels)
				if (std::binary_search(rel.begin(), rel.end(), rel1[0]))
					rel += rel1;
			if (!rel.empty())
				rels.push_back(std::move(rel));
		}
#else
		Poly1d rels_tmp(p_buffer->second.size());
		std::vector<std::future<void>> futures;
		for (int i = 0; i < NUM_THREADS; ++i)
			futures.push_back(std::async(std::launch::async, BatchReduceV2, std::ref(gb), std::ref(p_buffer->second), std::ref(rels_tmp), i));
		for (int i = 0; i < NUM_THREADS; ++i)
			futures[i].wait();
		Poly1d rels;
		for (auto& rel : rels_tmp) {
			for (Poly& rel1 : rels)
				if (std::binary_search(rel.begin(), rel.end(), rel1[0]))
					rel += rel1;
			if (!rel.empty())
				rels.push_back(std::move(rel));
		}
#endif

		/* Add these relations */
		for (auto& rel : rels) {
			if (!rel.empty()) {
				for (int i = 0; i < (int)gb.size(); ++i) {
					if (gcd_nonzero(rel[0], gb[i][0])) {
						Mon gcd = GCD(rel[0], gb[i][0]);
						int deg_new_rel = p_buffer->first + get_deg(gb[i][0], gen_degs, gen_degs1) - get_deg(gcd, gen_degs, gen_degs1);
						if (deg_max == -1 || deg_new_rel <= deg_max)
							buffer[deg_new_rel].push_back(std::make_unique<GcdEleBuffer>(std::move(gcd), i, (int)gb.size()));
					}
				}
				gb.push_back(std::move(rel));
			}
		}
	}
	buffer.erase(buffer.begin(), p_buffer);
}

GbBufferV2 GenerateBufferV2(const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t, int t_max)
{
	GbBufferV2 buffer;
	for (auto pg1 = gb.begin(); pg1 != gb.end(); ++pg1) {
		for (auto pg2 = pg1 + 1; pg2 != gb.end(); ++pg2) {
			Mon gcd = GCD(pg1->front(), pg2->front());
			if (!gcd.empty()) {
				int deg_new_rel = get_deg(pg1->front(), gen_degs, gen_degs1) + get_deg(pg2->front(), gen_degs, gen_degs1) - get_deg(gcd, gen_degs, gen_degs1);
				if (t <= deg_new_rel && deg_new_rel <= t_max)
					buffer[deg_new_rel].push_back(std::make_unique<GcdEleBuffer>(std::move(gcd), (int)(pg1 - gb.begin()), (int)(pg2 - gb.begin())));
			}
		}
	}
	return buffer;
}


/**********************************************************
* Algorithms that use Groebner basis
**********************************************************/

/* Compute the generating set of `vectors` inplace */
Poly2d& indecomposables(const Poly1d& gb, Poly2d& vectors, const array& gen_degs, const array& basis_degs)
{
	if (vectors.empty())
		return vectors;
	Poly1d gb1 = gb;

	/* Convert each vector v into a relation \\sum vi x_{-i-1} */
	Poly1d rels;
	array degs;
	for (const Poly1d& v : vectors) {
		Poly rel;
		for (int i = 0; i < basis_degs.size(); ++i)
			if (!v[i].empty())
				rel += v[i] * Mon{ {-i - 1, 1} };
		degs.push_back(get_deg(rel[0], gen_degs, basis_degs));
		rels.push_back(std::move(rel));
	}
	array indices = range((int)vectors.size());
	std::sort(indices.begin(), indices.end(), [&degs](int i, int j) {return degs[i] < degs[j]; });

	/* Add relations ordered by degree to gb1 */
	GbBuffer buffer;
	int deg_max = degs[indices.back()];
	for (int i : indices) {
		AddRelsM(gb1, buffer, gen_degs, basis_degs, degs[i], deg_max);
		Poly rel = Reduce(rels[i], gb1);
		if (!rel.empty())
			buffer[degs[i]].push_back(std::move(rel));
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
	AddRels(gb1, rels, gen_degs, gen_degs1, deg_max);

	/* Extract linear relations from gb1 */
	for (const Poly& g : gb1) {
		if (g[0][0].gen < 0) {
			Poly1d ann;
			ann.resize(N);
			for (const Mon& m : g) {
				MonInd p = m.begin();
				for (; p != m.end() && p->gen < 0; ++p);
				Mon m1(m.begin(), p), m2(p, m.end());
				ann[size_t(-m1[0].gen) - 1] += Reduce(mul(evaluate({ div(m1, { {m1[0].gen, 1} }) }, [&polys](int i) {return polys[size_t(-i) - 1]; }, gb), m2), gb);
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


} /* namespace grbn */