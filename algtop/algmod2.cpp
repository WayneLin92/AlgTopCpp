#include "mymath.h"


array mul(const array& mon1, const array& mon2)
{
	array result;
	size_t k = 0, l = 0;
	while (k < mon1.size() && l < mon2.size()) {
		if (mon1[k] < mon2[l]) {
			result.push_back(mon1[k]);
			result.push_back(mon1[k + 1]);
			k += 2;
		}
		else if (mon1[k] > mon2[l]) {
			result.push_back(mon2[l]);
			result.push_back(mon2[l + 1]);
			l += 2;
		}
		else {
			result.push_back(mon1[k]);
			result.push_back(mon1[k + 1] + mon2[l + 1]);
			k += 2;
			l += 2;
		}
	}
	if (k < mon1.size())
		result.insert(result.end(), mon1.begin() + k, mon1.end());
	else
		result.insert(result.end(), mon2.begin() + l, mon2.end());
	return result;
}

array div(const array& mon1, const array& mon2)
{
	array result;
	size_t k = 0, l = 0;
	while (k < mon1.size() && l < mon2.size()) {
		if (mon1[k] < mon2[l]) {
			result.push_back(mon1[k]);
			result.push_back(mon1[k + 1]);
			k += 2;
		}
#ifdef _DEBUG
		else if (mon1[k] > mon2[l]) {
			std::cout << "mon1/mon2 not divisible!\n";
			throw "1227de8e-31a8-40ae-9ad3-85b1cd6872cf";
		}
#endif
		else if (mon1[k + 1] > mon2[l + 1]) {
			result.push_back(mon1[k]);
			result.push_back(mon1[k + 1] - mon2[l + 1]);
			k += 2;
			l += 2;
		}
		else if (mon1[k + 1] == mon2[l + 1]) {
			k += 2;
			l += 2;
		}
#ifdef _DEBUG
		else {
			std::cout << "mon1/mon2 not divisible!\n";
			throw "a9c74ef9-c5ef-484a-8b53-261be63349e3";
		}
#endif
	}
#ifdef _DEBUG
	if (l < mon2.size()) {
		std::cout << "mon1/mon2 not divisible!\n";
		throw "6cdd66bd-597e-4c3b-ab1d-1bf8014d84a0";
	}
	else
#endif
		result.insert(result.end(), mon1.begin() + k, mon1.end());
	return result;
}

array2d add(const array2d& poly1, const array2d& poly2) {
	array2d result;
	std::set_symmetric_difference(poly1.begin(), poly1.end(), poly2.begin(), poly2.end(),
		std::back_inserter(result), cmp_mons);
	return result;
}

array2d mul(const array2d& poly, const array& mon) {
	array2d result;
	for (const array& m : poly)
		result.push_back(mul(m, mon));
	return result;
}

array2d mul(const array2d& poly1, const array2d& poly2) {
	array2d result;
	for (const array& mon2 : poly2)
		result = add(result, mul(poly1, mon2));
	return result;
}

array2d pow(const array2d& poly, int n)
{
	array2d result = { {} };
	if (n == 0)
		return result;
	array2d power = poly;
	while (n) {
		if (n & 1)
			result = mul(result, power);
		n >>= 1;
		if (n)
			power = mul(power, power);
	}
	return result;
}

array2d pow(const array2d& poly, int n, const array3d& gb)
{
	array2d result = { {} };
	if (n == 0)
		return result;
	array2d power = poly;
	while (n) {
		if (n & 1)
			result = reduce(mul(result, power), gb);
		n >>= 1;
		if (n)
			power = reduce(mul(power, power), gb);
	}
	return result;
}

bool cmp_mons(const array& m1, const array& m2)
{
	if (m1.empty())  /* 1 is the biggest monomial */
		return false;
	else if (m2.empty())
		return true;
	size_t i;
	for (i = 0; i < m1.size() && i < m2.size(); i += 2) {
		if (m1[i] > m2[i])
			return true;
		else if (m1[i] < m2[i])
			return false;
		else if (m1[i + 1] < m2[i + 1])
			return true;
		else if (m1[i + 1] > m2[i + 1])
			return false;
	}
	if (i < m2.size())
		return true;
	else
		return false;
}

bool divides(const array& m1, const array& m2)
{
	size_t k = 0, l = 0;
	while (k < m1.size() && l < m2.size()) {
		if (m1[k] < m2[l])
			return false;
		else if (m1[k] > m2[l])
			l += 2;
		else if (m1[k + 1] > m2[l + 1])
			return false;
		else {
			k += 2;
			l += 2;
		}
	}
	if (k < m1.size())
		return false;
	return true;
}

array pow(const array& m, int e)
{
	if (e == 0)
		return array();
	else if (e == 1)
		return m;
	array result;
	for (size_t i = 0; i < m.size(); i += 2) {
		result.push_back(m[i]);
		result.push_back(m[i + 1] * e);
	}
	return result;
}

int log(const array& m1, const array& m2)
{
	if (m2.empty()) {
		std::cout << "log with 0 base!\n";
		throw "f50d7f56-8ca7-4efe-9b23-3c9cde40d069";
	}
	int q = -1;

	/* Compute q */
	size_t k = 0, l = 0;
	while (k < m1.size() && l < m2.size()) {
		if (m1[k] < m2[l])
			k += 2;
		else if (m1[k] > m2[l]) {
			q = 0;
			break;
		}
		else if (m1[k + 1] < m2[l + 1]) {
			q = 0;
			break;
		}
		else {
			int q1 = m1[k + 1] / m2[l + 1];
			if (q == -1 || q > q1)
				q = q1;
			k += 2;
			l += 2;
		}
	}
	if (l < m2.size())
		q = 0;
	return q;
}

array gcd(const array& m1, const array& m2)
{
	array result;
	size_t k = 0, l = 0;
	while (k < m1.size() && l < m2.size()) {
		if (m1[k] < m2[l])
			k += 2;
		else if (m1[k] > m2[l])
			l += 2;
		else {
			result.push_back(m1[k]);
			result.push_back(std::min(m1[k + 1], m2[l + 1]));
			k += 2;
			l += 2;
		}
	}
	return result;
}

bool gcd_nonzero(const array& mon1, const array& mon2)
{
	size_t k = 0, l = 0;
	while (k < mon1.size() && l < mon2.size()) {
		if (mon1[k] < mon2[l])
			k += 2;
		else if (mon1[k] > mon2[l])
			l += 2;
		else
			return true;
	}
	return false;
}

array lcm(const array& mon1, const array& mon2)
{
	array result;
	size_t k = 0, l = 0;
	while (k < mon1.size() && l < mon2.size()) {
		if (mon1[k] < mon2[l]) {
			result.push_back(mon1[k]);
			result.push_back(mon1[k + 1]);
			k += 2;
		}
		else if (mon1[k] > mon2[l]) {
			result.push_back(mon2[l]);
			result.push_back(mon2[l + 1]);
			l += 2;
		}
		else {
			result.push_back(mon1[k]);
			result.push_back(std::max(mon1[k + 1], mon2[l + 1]));
			k += 2;
			l += 2;
		}
	}
	if (k < mon1.size())
		result.insert(result.end(), mon1.begin() + k, mon1.end());
	else
		result.insert(result.end(), mon2.begin() + l, mon2.end());
	return result;
}

array2d get_diff(const array& mon, const array3d& diffs)
{
	array2d result;
	for (size_t k = 0; k < mon.size(); k += 2) {
		if (mon[k + 1] % 2) {
			array m1 = div(mon, { mon[k], 1 });
			result = add(result, mul({ std::move(m1) }, diffs[mon[k]]));
		}
	}
	return result;
}

array2d get_diff(const array2d& poly, const array3d& diffs)
{
	array2d result;
	for (const array& m : poly) {
		for (size_t k = 0; k < m.size(); k += 2) {
			if (m[k + 1] % 2) {
				array m1 = div(m, { m[k], m[k + 1] });
				result = add(result, mul({ std::move(m1) }, diffs[m[k]]));
			}
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

array residue(const array2d& spaceV, array&& v)
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
		throw "6a4fe8a1-608e-466c-ab5e-5fae459ce1b9";
	}
#endif
	return result;
}

array2d quotient_space(const array2d& spaceV, const array2d& spaceW)
{
	array2d quotient;
	size_t dimQuo = spaceV.size() - spaceW.size();
	for (size_t i = 0; i < spaceV.size() && quotient.size() < dimQuo; i++) {
		auto v1 = residue(quotient, residue(spaceW, spaceV[i]));
		if (!v1.empty())
			quotient.push_back(std::move(v1));
	}
#if _DEBUG
	if (quotient.size() != dimQuo) {
		std::cerr << "W is not a subspace of V!\n";
		throw "cec7f701-0911-482a-a63f-1caaa646591b";
	}
#endif
	return quotient;
}

/******** Groebner Basis ********/

array2d reduce(array2d poly, const array3d& gb)
{
	array2d result;
	auto pbegin = poly.begin(); auto pend = poly.end();
	while (pbegin != pend) {
		array3d::const_iterator pGb = gb.begin();
		for (; pGb != gb.end(); ++pGb)
			if (divides((*pGb)[0], *pbegin))
				break;
		if (pGb == gb.end())
			result.push_back(std::move(*pbegin++));
		else {
			array q = div(*pbegin, (*pGb)[0]);
			array2d rel1 = mul(*pGb, q);
			array2d poly1;
			std::set_symmetric_difference(pbegin, pend, rel1.begin(), rel1.end(),
				std::back_inserter(poly1), cmp_mons);
			poly = std::move(poly1);
			pbegin = poly.begin(); pend = poly.end();
		}
	}
	return result;
}

/* Comsume relations from heap that is at most in degree `deg` while adding new relations to heap that is at most in degree `deg_max`. */
template <typename Fn>
void add_rels(array3d& gb, std::vector<rel_heap_t>& heap, Fn _get_deg, int deg, int deg_max)
{
	while (!heap.empty() && (deg == -1 ? (deg_max == -1 || heap.front().deg <= deg_max) : heap.front().deg <= deg)) {
		std::pop_heap(heap.begin(), heap.end(), cmp_heap_rels);
		rel_heap_t heap_ele = std::move(heap.back());
		heap.pop_back();

		array2d rel = reduce(heap_ele.poly, gb);
		if (!rel.empty()) {
			for (array2d& g : gb) {
				if (gcd_nonzero(rel[0], g[0])) {
					if (divides(rel[0], g[0])) {
						array q = div(g[0], rel[0]);
						array2d new_rel = add(mul(rel, q), g);
						if (!new_rel.empty()) {
							heap.push_back(rel_heap_t{ std::move(new_rel), _get_deg(g[0]) });
							std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
						}
						g.clear();
					}
					else {
						array mlcm = lcm(rel[0], g[0]);
						int deg_new_rel = _get_deg(mlcm);
						if (deg_max == -1 || deg_new_rel <= deg_max) {
							array q_r = div(mlcm, rel[0]);
							array q_g = div(mlcm, g[0]);
							array2d new_rel = add(mul(rel, q_r), mul(g, q_g));

							heap.push_back(rel_heap_t{ std::move(new_rel), deg_new_rel });
							std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
						}
					}
				}
			}
			gb.erase(std::remove_if(gb.begin(), gb.end(), [](const array2d& g) {return g.empty(); }), gb.end());
			gb.push_back(std::move(rel));
		}
	}
}

template <typename Fn>
void add_rels(array3d& gb, const array3d& rels, Fn _get_deg, int deg_max)
{
	std::vector<rel_heap_t> heap;
	for (const array2d& rel : rels) {
		if (!rel.empty()) {
			heap.push_back(rel_heap_t{ rel, _get_deg(rel[0]) });
			std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
		}
	}
	add_rels(gb, heap, _get_deg, -1, deg_max);
}

void add_rels(array3d& gb, std::vector<rel_heap_t>& heap, const array& gen_degs, int t, int deg_max)
{
	add_rels(gb, heap, [&gen_degs](const array& m) {return get_deg(m, gen_degs); }, t, deg_max);
}

void add_rels(array3d& gb, const array3d& rels, const array& gen_degs, int deg_max)
{
	std::vector<rel_heap_t> heap;
	for (const array2d& rel : rels) {
		if (!rel.empty()) {
			heap.push_back(rel_heap_t{ rel, get_deg(rel[0], gen_degs) });
			std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
		}
	}
	add_rels(gb, heap, [&gen_degs](const array& m) {return get_deg(m, gen_degs); }, -1, deg_max);
}

template <typename Fn>
void add_rels_freemodule(array3d& gb, std::vector<rel_heap_t>& heap, Fn get_deg, int deg, int deg_max)
{
	while (!heap.empty() && (deg == -1 ? (deg_max == -1 || heap.front().deg <= deg_max) : heap.front().deg <= deg)) {
		std::pop_heap(heap.begin(), heap.end(), cmp_heap_rels);
		rel_heap_t heap_ele = std::move(heap.back());
		heap.pop_back();

		array2d rel = reduce(heap_ele.poly, gb);
		if (!rel.empty()) {
			for (array2d& g : gb) {
				if ((rel[0][0] >= 0 || g[0][0] >= 0 || rel[0][0] == g[0][0]) && gcd_nonzero(rel[0], g[0])) {
					if (divides(rel[0], g[0])) {
						array q = div(g[0], rel[0]);
						array2d new_rel = add(mul(rel, q), g);
						if (!new_rel.empty()) {
							heap.push_back(rel_heap_t{ std::move(new_rel), get_deg(g[0]) });
							std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
						}
						g.clear();
					}
					else {
						array mlcm = lcm(rel[0], g[0]);
						int deg_new_rel = get_deg(mlcm);
						if (deg_max == -1 || deg_new_rel <= deg_max) {
							array q_r = div(mlcm, rel[0]);
							array q_g = div(mlcm, g[0]);
							array2d new_rel = add(mul(rel, q_r), mul(g, q_g));

							heap.push_back(rel_heap_t{ std::move(new_rel), deg_new_rel });
							std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
						}
					}
				}
			}
			gb.erase(std::remove_if(gb.begin(), gb.end(), [](const array2d& g) {return g.empty(); }), gb.end());
			gb.push_back(std::move(rel));
		}
	}
}

array4d& indecomposables(const array3d& gb, array4d& vectors, const array& gen_degs, const array& basis_degs)
{
	if (vectors.empty())
		return vectors;
	array3d gb1 = gb;
	int N = (int)basis_degs.size();
	auto get_deg = [&gen_degs, &basis_degs, &N](const array& mon) {
		int result = 0;
		for (size_t i = 0; i < mon.size(); i += 2)
			result += (mon[i] >= 0 ? gen_degs[mon[i]] : basis_degs[mon[i] + size_t(N)]) * mon[i + 1];
		return result;
	};

	array3d rels;
	array degs;
	for (const array3d& v : vectors) {
		array2d rel;
		for (int i = 0; i < N; ++i) {
			if (!v[i].empty())
				rel = add(rel, mul(v[i], { i - N, 1 }));
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
		array2d rel = reduce(rels[i], gb1);
		if (!rel.empty()) {
			heap.push_back(rel_heap_t{ std::move(rel), deg });
			std::push_heap(heap.begin(), heap.end(), cmp_heap_rels);
		}
		else
			vectors[i].clear();
	}

	vectors.erase(std::remove_if(vectors.begin(), vectors.end(), [](const array3d& v) {return v.empty(); }), vectors.end());
	return vectors;
}

array4d ann_seq(const array3d& gb, const array3d& polys, const array& gen_degs, int deg_max)
{
	array4d result;
	if (polys.empty())
		return result;
	array3d rels;
	array gen_degs1;
	int N = (int)polys.size();
	for (int i = 0; i < N; ++i) {
		array2d p = polys[i];
		gen_degs1.push_back(get_deg(p, gen_degs));
		p.push_back({ i - N, 1 });
		rels.push_back(std::move(p));
	}
	array3d gb1 = gb;
	add_rels(gb1, rels, [&gen_degs, &gen_degs1, &N](const array& mon) {
		int result_ = 0;
		for (size_t i = 0; i < mon.size(); i += 2)
			result_ += (mon[i] >= 0 ? gen_degs[mon[i]] : gen_degs1[mon[i] + size_t(N)]) * mon[i + 1];
		return result_;
		}, deg_max);
	for (const array2d& g : gb1) {
		if (g[0][0] < 0) {
			array3d result_i;
			result_i.resize(N);
			for (const array& m : g) {
				auto p = m.begin();
				for (; p < m.end() && *p < 0; p += 2);
				array m1(m.begin(), p), m2(p, m.end());
				result_i[m1[0] + size_t(N)] = add(result_i[m1[0] + size_t(N)], reduce(mul(evaluate({ div(m1, { m1[0], 1 }) }, [&polys, &N](int i) {return polys[i + size_t(N)]; }, gb), m2), gb));
			}
			result.push_back(std::move(result_i));
		}
	}
	for (int i = 0; i < N; ++i) {
		for (int j = i + 1; j < N; j++) {
			if (gen_degs1[i] + gen_degs1[j] <= deg_max) {
				array3d result_i;
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