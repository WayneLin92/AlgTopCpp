#include "mymath.h"

/******** Monomials and Polynomials ********/

Mon mul(const Mon& mon1, const Mon& mon2)
{
	Mon result;
	MonInd k = mon1.begin(), l = mon2.begin();
	while (k != mon1.end() && l != mon2.end()) {
		if (k->gen < l->gen)
			result.push_back(*k++);
		else if (k->gen > l->gen)
			result.push_back(*l++);
		else {
			result.emplace_back(k->gen, k->exp + l->exp);
			k++; l++;
		}
	}
	if (k != mon1.end())
		result.insert(result.end(), k, mon1.end());
	else
		result.insert(result.end(), l, mon2.end());
	return result;
}

Mon div(const Mon& mon1, const Mon& mon2)
{
	Mon result;
	MonInd k = mon1.begin(), l = mon2.begin();
	while (k != mon1.end() && l != mon2.end()) {
		if (k->gen < l->gen)
			result.push_back(*k++);
#ifdef _DEBUG
		else if (k->gen > l->gen){
			std::cout << "mon1/mon2 not divisible!\n";
			throw "1227de8e";
		}
#endif
		else if (k->exp > l->exp) {
			result.emplace_back(k->gen, k->exp - l->exp);
			k++; l++;
		}
		else if (k->exp == l->exp) {
			k++;
			l++;
		}
#ifdef _DEBUG
		else {
			std::cout << "mon1/mon2 not divisible!\n";
			throw "a9c74ef9";
		}
#endif
	}
#ifdef _DEBUG
	if (l != mon2.end()) {
		std::cout << "mon1/mon2 not divisible!\n";
		throw "6cdd66bd";
	}
	else
#endif
		result.insert(result.end(), k, mon1.end());
	return result;
}

Poly add(const Poly& poly1, const Poly& poly2) {
	Poly result;
	std::set_symmetric_difference(poly1.begin(), poly1.end(), poly2.begin(), poly2.end(),
		std::back_inserter(result));
	return result;
}

Poly mul(const Poly& poly, const Mon& mon) {
	Poly result;
	for (const Mon& m : poly)
		result.push_back(mul(m, mon));
	return result;
}

Poly mul(const Poly& poly1, const Poly& poly2) {
	Poly result;
	for (const Mon& mon2 : poly2)
		result = add(result, mul(poly1, mon2));
	return result;
}

Poly pow(const Poly& poly, int n)
{
	Poly result = { {} };
	if (n == 0)
		return result;
	Poly power = poly;
	while (n) {
		if (n & 1)
			result = mul(result, power);
		n >>= 1;
		if (n)
			power = mul(power, power);
	}
	return result;
}

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

bool divides(const Mon& mon1, const Mon& mon2)
{
	MonInd k = mon1.begin(), l = mon2.begin();
	while (k != mon1.end() && l != mon2.end()) {
		if (k->gen < l->gen)
			return false;
		else if (k->gen > l->gen)
			l++;
		else if (k->exp > l->exp)
			return false;
		else {
			k++;
			l++;
		}
	}
	if (k != mon1.end())
		return false;
	return true;
}

Mon pow(const Mon& m, int e)
{
	if (e == 0)
		return {};
	else if (e == 1)
		return m;
	Mon result;
	for (MonInd p = m.begin(); p != m.end(); ++p)
		result.emplace_back(p->gen, p->exp * e);
	return result;
}

int log(const Mon& mon1, const Mon& mon2)
{
	if (mon2.empty()) {
		std::cout << "log with 0 base!\n";
		throw "f50d7f56";
	}
	int q = -1;

	/* Compute q */
	MonInd k = mon1.begin(), l = mon2.begin();
	while (k != mon1.end() && l != mon2.end()) {
		if (k->gen < l->gen)
			k++;
		else if (k->gen > l->gen) {
			q = 0;
			break;
		}
		else if (k->exp < l->exp) {
			q = 0;
			break;
		}
		else {
			int q1 = k->exp / l->exp;
			if (q == -1 || q > q1)
				q = q1;
			k++;
			l++;
		}
	}
	if (l != mon2.end())
		q = 0;
	return q;
}

Mon gcd(const Mon& mon1, const Mon& mon2)
{
	Mon result;
	MonInd k = mon1.begin(), l = mon2.begin();
	while (k != mon1.end() && l != mon2.end()) {
		if (k->gen < l->gen)
			k++;
		else if (k->gen > l->gen)
			l++;
		else {
			result.emplace_back(k->gen, std::min(k->exp, l->exp));
			k++; l++;
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

Mon lcm(const Mon& mon1, const Mon& mon2)
{
	Mon result;
	MonInd k = mon1.begin(), l = mon2.begin();
	while (k != mon1.end() && l != mon2.end()) {
		if (k->gen < l->gen)
			result.push_back(*k++);
		else if (k->gen > l->gen)
			result.push_back(*l++);
		else {
			result.emplace_back(k->gen, std::max(k->exp, l->exp));
			k++; l++;
		}
	}
	if (k != mon1.end())
		result.insert(result.end(), k, mon1.end());
	else
		result.insert(result.end(), l, mon2.end());
	return result;
}

Poly get_diff(const Mon& mon, const Poly1d& diffs)
{
	Poly result;
	for (MonInd k = mon.begin(); k != mon.end(); ++k) {
		if (k->exp % 2)
			result = add(result, mul(diffs[k->gen], div(mon, { { k->gen, 1 } })));
	}
	return result;
}

Poly get_diff(const Poly& poly, const Poly1d& diffs)
{
	Poly result;
	for (const Mon& mon : poly) {
		for (MonInd k = mon.begin(); k != mon.end(); ++k) {
			if (k->exp % 2)
				result = add(result, mul(diffs[k->gen], div(mon, { { k->gen, 1 } })));
		}
	}
	return result;
}

/******** Linear Algebra ********/

array2d& simplify_space(array2d& spaceV)
{
	for (size_t i = spaceV.size() - 1; i != -1; i--)
		for (size_t j = 0; j < i; j++)
			if (std::binary_search(spaceV[j].begin(), spaceV[j].end(), spaceV[i][0]))
				spaceV[j] = add_vectors(spaceV[j], spaceV[i]);
	return spaceV;
}

array residue(const array2d& spaceV, const array& v)
{
	array result(v);
	for (size_t i = 0; i < spaceV.size(); ++i)
		if (std::binary_search(result.begin(), result.end(), spaceV[i][0]))
			result = add_vectors(result, spaceV[i]);
	return result;
}

array residue(const array2d& spaceV, array&& v)//
{
	array result(v);
	for (size_t i = 0; i < spaceV.size(); ++i)
		if (std::binary_search(result.begin(), result.end(), spaceV[i][0]))
			result = add_vectors(result, spaceV[i]);
	return result;
}

inline void add_map(array2d& spaceV, const array& v)
{
	array v1 = residue(spaceV, v);
	if (!v1.empty())
		spaceV.push_back(std::move(v1));
}

void set_linear_map(const array2d& fx, array2d& image, array2d& kernel, array2d& g)
{
	/* f(g[i]) = image[i] */
	for (size_t i = 0; i < fx.size(); ++i) {
		array src = { int(i) };
		array tgt(fx[i]);
		for (size_t j = 0; j < image.size(); j++) {
			if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
				tgt = add_vectors(tgt, image[j]);
				src = add_vectors(src, g[j]);
			}
		}
		if (tgt.empty())
			add_map(kernel, src);
		else {
			image.push_back(std::move(tgt));
			g.push_back(std::move(src));
		}
	}
}

void set_linear_map(const array& x, const array2d& fx, array2d& image, array2d& kernel, array2d& g)
{
	/* f(g[i]) = image[i] */
	for (size_t i = 0; i < fx.size(); ++i) {
		array src = { x[i] };
		array tgt(fx[i]);
		for (size_t j = 0; j < image.size(); j++) {
			if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
				tgt = add_vectors(tgt, image[j]);
				src = add_vectors(src, g[j]);
			}
		}
		if (tgt.empty())
			add_map(kernel, src);
		else {
			image.push_back(std::move(tgt));
			g.push_back(std::move(src));
		}
	}
}

array get_image(const array2d& spaceV, const array2d& f, const array& v)
{
	array result;
	array v1(v);
	for (size_t i = 0; i < spaceV.size() && !v1.empty(); ++i)
		if (std::binary_search(v1.begin(), v1.end(), spaceV[i][0])) {
			v1 = add_vectors(v1, spaceV[i]);
			result = add_vectors(result, f[i]);
		}
#if _DEBUG
	if (!v1.empty()) {
		std::cerr << "v is not in space V!\n";
		throw "6a4fe8a1";
	}
#endif
	return result;
}

array2d quotient_space(const array2d& spaceV, const array2d& spaceW)
{
	array2d quotient;
	size_t dimQuo = spaceV.size() - spaceW.size();
#if _DEBUG
	for (size_t i = 0; i < spaceV.size(); i++) {
#else
	for (size_t i = 0; i < spaceV.size() && quotient.size() < dimQuo; i++) {
#endif
		auto v1 = residue(quotient, residue(spaceW, spaceV[i]));
		if (!v1.empty())
			quotient.push_back(std::move(v1));
	}
#if _DEBUG
	if (quotient.size() != dimQuo) {
		std::cerr << "W is not a subspace of V!\n";
		throw "cec7f701";
	}
#endif
	return quotient;
}

/******** Groebner Basis ********/

Poly reduce(Poly poly, const Poly1d& gb)
{
	Poly result;
	auto pbegin = poly.begin(); auto pend = poly.end();
	while (pbegin != pend) {
		auto pGb = gb.begin();
		for (; pGb != gb.end(); ++pGb)
			if (divides((*pGb)[0], *pbegin))
				break;
		if (pGb == gb.end())
			result.push_back(std::move(*pbegin++));
		else {
			Mon q = div(*pbegin, (*pGb)[0]);
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

/* Comsume relations from heap that is at most in degree `deg` while adding new relations to heap that is at most in degree `deg_max`. */
template <typename Fn>
void add_rels(Poly1d& gb, std::vector<rel_heap_t>& heap, Fn _get_deg, int deg, int deg_max)
{
	while (!heap.empty() && (deg == -1 ? (deg_max == -1 || heap.front().deg <= deg_max) : heap.front().deg <= deg)) {
		std::pop_heap(heap.begin(), heap.end(), cmp_heap_rels);
		rel_heap_t heap_ele = std::move(heap.back());
		heap.pop_back();

		Poly rel = reduce(heap_ele.poly, gb);
		if (!rel.empty()) {
			for (Poly& g : gb) {
				if (gcd_nonzero(rel[0], g[0])) {
					if (divides(rel[0], g[0])) {
						Mon q = div(g[0], rel[0]);
						Poly new_rel = add(mul(rel, q), g);
						if (!new_rel.empty()) {
							heap.push_back(rel_heap_t{ std::move(new_rel), _get_deg(g[0]) });
							std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
						}
						g.clear();
					}
					else {
						Mon mlcm = lcm(rel[0], g[0]);
						int deg_new_rel = _get_deg(mlcm);
						if (deg_max == -1 || deg_new_rel <= deg_max) {
							Mon q_r = div(mlcm, rel[0]);
							Mon q_g = div(mlcm, g[0]);
							Poly new_rel = add(mul(rel, q_r), mul(g, q_g));

							heap.push_back(rel_heap_t{ std::move(new_rel), deg_new_rel });
							std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
						}
					}
				}
			}
			gb.erase(std::remove_if(gb.begin(), gb.end(), [](const Poly& g) {return g.empty(); }), gb.end());
			gb.push_back(std::move(rel));
		}
	}
}

template <typename Fn>
void add_rels(Poly1d& gb, const Poly1d& rels, Fn _get_deg, int deg_max)
{
	std::vector<rel_heap_t> heap;
	for (const Poly& rel : rels) {
		if (!rel.empty()) {
			heap.push_back(rel_heap_t{ rel, _get_deg(rel[0]) });
			std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
		}
	}
	add_rels(gb, heap, _get_deg, -1, deg_max);
}

void add_rels(Poly1d& gb, std::vector<rel_heap_t>& heap, const array& gen_degs, int t, int deg_max)
{
	add_rels(gb, heap, [&gen_degs](const Mon& m) {return get_deg(m, gen_degs); }, t, deg_max);
}

void add_rels(Poly1d& gb, const Poly1d& rels, const array& gen_degs, int deg_max)
{
	std::vector<rel_heap_t> heap;
	for (const Poly& rel : rels) {
		if (!rel.empty()) {
			heap.push_back(rel_heap_t{ rel, get_deg(rel[0], gen_degs) });
			std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
		}
	}
	add_rels(gb, heap, [&gen_degs](const Mon& m) {return get_deg(m, gen_degs); }, -1, deg_max);
}

template <typename Fn>
void add_rels_freemodule(Poly1d& gb, std::vector<rel_heap_t>& heap, Fn get_deg, int deg, int deg_max)
{
	while (!heap.empty() && (deg == -1 ? (deg_max == -1 || heap.front().deg <= deg_max) : heap.front().deg <= deg)) {
		std::pop_heap(heap.begin(), heap.end(), cmp_heap_rels);
		rel_heap_t heap_ele = std::move(heap.back());
		heap.pop_back();

		Poly rel = reduce(heap_ele.poly, gb);
		if (!rel.empty()) {
			for (Poly& g : gb) {
				if ((rel[0][0].gen >= 0 || g[0][0].gen >= 0 || rel[0][0].gen == g[0][0].gen) && gcd_nonzero(rel[0], g[0])) {
					if (divides(rel[0], g[0])) {
						Mon q = div(g[0], rel[0]);
						Poly new_rel = add(mul(rel, q), g);
						if (!new_rel.empty()) {
							heap.push_back(rel_heap_t{ std::move(new_rel), get_deg(g[0]) });
							std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
						}
						g.clear();
					}
					else {
						Mon mlcm = lcm(rel[0], g[0]);
						int deg_new_rel = get_deg(mlcm);
						if (deg_max == -1 || deg_new_rel <= deg_max) {
							Mon q_r = div(mlcm, rel[0]);
							Mon q_g = div(mlcm, g[0]);
							Poly new_rel = add(mul(rel, q_r), mul(g, q_g));

							heap.push_back(rel_heap_t{ std::move(new_rel), deg_new_rel });
							std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
						}
					}
				}
			}
			gb.erase(std::remove_if(gb.begin(), gb.end(), [](const Poly& g) {return g.empty(); }), gb.end());
			gb.push_back(std::move(rel));
		}
	}
}

Poly2d& indecomposables(const Poly1d& gb, Poly2d& vectors, const array& gen_degs, const array& basis_degs)
{
	if (vectors.empty())
		return vectors;
	Poly1d gb1 = gb;
	int N = (int)basis_degs.size();
	auto get_deg = [&gen_degs, &basis_degs, &N](const Mon& mon) {
		int result = 0;
		for (MonInd p = mon.begin(); p != mon.end(); ++p)
			result += (p->gen >= 0 ? gen_degs[p->gen] : basis_degs[p->gen + size_t(N)]) * p->exp;
		return result;
	};

	Poly1d rels;
	array degs;
	for (const Poly1d& v : vectors) {
		Poly rel;
		for (int i = 0; i < N; ++i) {
			if (!v[i].empty())
				rel = add(rel, mul(v[i], { {i - N, 1} }));
		}
		degs.push_back(get_deg(rel[0]));
		rels.push_back(std::move(rel));
	}
	array indices = range((int)vectors.size());
	std::sort(indices.begin(), indices.end(), [&degs](int i, int j) {return degs[i] < degs[j]; });

	std::vector<rel_heap_t> heap;
	int deg_max = get_deg(rels[indices.back()][0]);
	for (int i : indices) {
		int deg = get_deg(rels[i][0]);
		add_rels_freemodule(gb1, heap, get_deg, deg, deg_max);
		Poly rel = reduce(rels[i], gb1);
		if (!rel.empty()) {
			heap.push_back(rel_heap_t{ std::move(rel), deg });
			std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
		}
		else
			vectors[i].clear();
	}

	vectors.erase(std::remove_if(vectors.begin(), vectors.end(), [](const Poly1d& v) {return v.empty(); }), vectors.end());
	return vectors;
}

Poly2d ann_seq(const Poly1d& gb, const Poly1d& polys, const array& gen_degs, int deg_max)
{
	Poly2d result;
	if (polys.empty())
		return result;
	Poly1d rels;
	array gen_degs1;
	int N = (int)polys.size();
	for (int i = 0; i < N; ++i) {
		Poly p = polys[i];
		gen_degs1.push_back(get_deg(p, gen_degs));
		p.push_back({ {i - N, 1} });
		rels.push_back(std::move(p));
	}
	Poly1d gb1 = gb;
	add_rels(gb1, rels, [&gen_degs, &gen_degs1, &N](const Mon& mon) {
		int result_ = 0;
		for (MonInd p = mon.begin(); p != mon.end(); ++p)
			result_ += (p->gen >= 0 ? gen_degs[p->gen] : gen_degs1[p->gen + size_t(N)]) * p->exp;
		return result_;
		}, deg_max);
	for (const Poly& g : gb1) {
		if (g[0][0].gen < 0) {
			Poly1d result_i;
			result_i.resize(N);
			for (const Mon& m : g) {
				MonInd p = m.begin();
				for (; p != m.end() && p->gen < 0; ++p);
				Mon m1(m.begin(), p), m2(p, m.end());
				result_i[m1[0].gen + size_t(N)] = add(result_i[m1[0].gen + size_t(N)], reduce(mul(evaluate({ div(m1, { {m1[0].gen, 1} }) }, [&polys, &N](int i) {return polys[i + size_t(N)]; }, gb), m2), gb));
			}
			result.push_back(std::move(result_i));
		}
	}
	for (int i = 0; i < N; ++i) {
		for (int j = i + 1; j < N; j++) {
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