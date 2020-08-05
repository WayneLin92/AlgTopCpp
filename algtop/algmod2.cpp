#include "mymath.h"

// operator overloading
array mul(const array& m1, const array& m2)
{
	array result;
	size_t k = 0, l = 0;
	while (k < m1.size() && l < m2.size()) {
		if (m1[k] < m2[l]) {
			result.push_back(m1[k]);
			result.push_back(m1[k + 1]);
			k += 2;
		}
		else if (m1[k] > m2[l]) {
			result.push_back(m2[l]);
			result.push_back(m2[l + 1]);
			l += 2;
		}
		else {
			result.push_back(m1[k]);
			result.push_back(m1[k + 1] + m2[l + 1]);
			k += 2;
			l += 2;
		}
	}
	if (k < m1.size())
		result.insert(result.end(), m1.begin() + k, m1.end());
	else
		result.insert(result.end(), m2.begin() + l, m2.end());
	return result;
}

array div(const array& m1, const array& m2)
{
	array result;
	size_t k = 0, l = 0;
	while (k < m1.size() && l < m2.size()) {
		if (m1[k] < m2[l]) {
			result.push_back(m1[k]);
			result.push_back(m1[k + 1]);
			k += 2;
		}
		else if (m1[k] > m2[l]) {
			std::cout << "m1/m2 not divisible!\n";
			throw "1227de8e-31a8-40ae-9ad3-85b1cd6872cf";
		}
		else if (m1[k + 1] > m2[l + 1]) {
			result.push_back(m1[k]);
			result.push_back(m1[k + 1] - m2[l + 1]);
			k += 2;
			l += 2;
		}
		else if (m1[k + 1] == m2[l + 1]) {
			k += 2;
			l += 2;
		}
		else {
			std::cout << "m1/m2 not divisible!\n";
			throw "a9c74ef9-c5ef-484a-8b53-261be63349e3";
		}
	}
	if (l < m2.size()) {
		std::cout << "m1/m2 not divisible!\n";
		throw "6cdd66bd-597e-4c3b-ab1d-1bf8014d84a0";
	}
	else
		result.insert(result.end(), m1.begin() + k, m1.end());
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

bool cmp_mons(const array& m1, const array& m2)
{
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
	array result;
	for (size_t i = 0; i < m.size(); i += 2) {
		result.push_back(m[i]);
		result.push_back(m[i + 1] * e);
	}
	return result;
}

log_result log(const array& m1, const array& m2)
{
	if (m2.empty()) {
		std::cout << "log with 0 base!\n";
		throw "f50d7f56-8ca7-4efe-9b23-3c9cde40d069";
	}
	array r;
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

	return q == 0 ? log_result({ 0, m1 }) : log_result({ q, div(m1, pow(m2, q)) });
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

bool gcd_nonzero(const array& m1, const array& m2)
{
	size_t k = 0, l = 0;
	while (k < m1.size() && l < m2.size()) {
		if (m1[k] < m2[l])
			k += 2;
		else if (m1[k] > m2[l])
			l += 2;
		else
			return true;
	}
	return false;
}

array lcm(const array& m1, const array& m2)
{
	array result;
	size_t k = 0, l = 0;
	while (k < m1.size() && l < m2.size()) {
		if (m1[k] < m2[l]) {
			result.push_back(m1[k]);
			result.push_back(m1[k + 1]);
			k += 2;
		}
		else if (m1[k] > m2[l]) {
			result.push_back(m2[l]);
			result.push_back(m2[l + 1]);
			l += 2;
		}
		else {
			result.push_back(m1[k]);
			result.push_back(std::max(m1[k + 1], m2[l + 1]));
			k += 2;
			l += 2;
		}
	}
	return result;
}

struct rel_heap_t {
	array2d poly;
	int deg;
};

array2d reduce(const array2d& poly, const array3d& rels)
{
	array2d result;
	array2d poly1(poly);
	auto pbegin = poly1.begin(); auto pend = poly1.end();
	while (pbegin != pend) {
		array3d::const_iterator p_rel = rels.begin();
		for (; p_rel != rels.end(); ++p_rel)
			if (divides((*p_rel)[0], *pbegin))
				break;
		if (p_rel == rels.end())
			result.push_back(std::move(*pbegin++));
		else {
			const array2d& rel = *p_rel;
			array q = div(*pbegin, rel[0]);
			array2d rel1 = mul(rel, q); // TODO: improve
			array2d poly2;
			std::set_symmetric_difference(pbegin, pend, rel1.begin(), rel1.end(),
				std::back_inserter(poly2), cmp_mons);
			poly1 = std::move(poly2);
			pbegin = poly1.begin(); pend = poly1.end();
		}
	}
	return result;
}

inline bool compare(const rel_heap_t& s1, const rel_heap_t& s2) { return s1.deg > s2.deg; }

void add_rel(array3d& rels, const array2d& rel, const array& gen_degs)
{
	if (rel.empty())
		return;
	std::vector<rel_heap_t> heap = { rel_heap_t{rel, deg(rel, gen_degs)} };
	while (!heap.empty()) {
		std::pop_heap(heap.begin(), heap.end(), compare);
		rel_heap_t heap_ele = std::move(heap.back());
		heap.pop_back();

		array2d r = reduce(heap_ele.poly, rels);
		if (!r.empty()) {
			for (array2d& r1 : rels) {
				if (gcd_nonzero(r[0], r1[0])) {
					if (divides(r[0], r1[0]))
						r1.clear();
					else {
						array mlcm = lcm(r[0], r1[0]);
						array q = div(mlcm, r[0]);
						array q1 = div(mlcm, r1[0]);
						array2d new_rel = add(mul(r, q), mul(r1, q1));

						heap.push_back(rel_heap_t{ new_rel, deg(mlcm, gen_degs) });
						std::push_heap(heap.begin(), heap.end(), [](const rel_heap_t& s1, const rel_heap_t& s2) { return s1.deg > s2.deg; });
					}
				}
			} 
			rels.erase(std::remove_if(rels.begin(), rels.end(), [](const array2d& p) {return p.empty(); }), rels.end());
			rels.push_back(std::move(r));
		}
	}
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
	for (size_t i = 0; i < spaceV.size(); i++)
		if (std::binary_search(result.begin(), result.end(), spaceV[i][0]))
			result = add_vectors(result, spaceV[i]);
	return result;
}

array residue(const array2d& spaceV, array&& v)
{
	array result(v);
	for (size_t i = 0; i < spaceV.size(); i++)
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

void get_image_kernel(const array2d& fx, array2d& image, array2d& kernel)
{
	/* f(g[i]) = image[i] */
	array2d g;
	for (size_t i = 0; i < fx.size(); i++) {
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

void get_image_kernel(const array& x, const array2d& fx, array2d& image, array2d& kernel)
{
	/* f(g[i]) = image[i] */
	array2d g;
	for (size_t i = 0; i < fx.size(); i++) {
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

array2d quotient_space(const array2d& spaceV, const array2d& spaceW)
{
	array2d quotient;
	size_t dimQuo = spaceV.size() - spaceW.size();
	for (const auto& v : spaceV) {
		auto v1 = residue(quotient, residue(spaceW, v));
		if (!v1.empty()) {
			quotient.push_back(std::move(v1));
#if !_DEBUG
			if (quotient.size() == dimQuo)
				return quotient;
#endif
		}
	}
#if _DEBUG
	if (quotient.size() != dimQuo) {
		std::cerr << "W is not a subspace of V!\n";
		throw "cec7f701-0911-482a-a63f-1caaa646591b";
	}
	return quotient;
#endif
}