#include "mymath.h"
#include <sstream> ////
#include <string> ////

// operator overloading
MonRaw mul(const MonRaw& m1, const MonRaw& m2)
{
	MonRaw result;
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

MonRaw div(const MonRaw& m1, const MonRaw& m2)
{
	MonRaw result;
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

PolyRaw add(const PolyRaw& poly1, const PolyRaw& poly2) {
	PolyRaw result;
	std::set_symmetric_difference(poly1.begin(), poly1.end(), poly2.begin(), poly2.end(),
		std::back_inserter(result), cmp_mons);
	return result;
}

PolyRaw mul(const PolyRaw& poly, const MonRaw& mon) {
	PolyRaw result;
	for (const MonRaw& m : poly)
		result.push_back(mul(m, mon));
	return result;
}

PolyRaw mul(const PolyRaw& poly1, const PolyRaw& poly2) {
	PolyRaw result;
	for (const MonRaw& mon2 : poly2)
		result = add(result, mul(poly1, mon2));
	return result;
}

bool cmp_mons(const MonRaw& m1, const MonRaw& m2)
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

bool divides(const MonRaw& m1, const MonRaw& m2)
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

MonRaw pow(const MonRaw& m, int e)
{
	MonRaw result;
	for (size_t i = 0; i < m.size(); i += 2) {
		result.push_back(m[i]);
		result.push_back(m[i + 1] * e);
	}
	return result;
}

log_result log(const MonRaw& m1, const MonRaw& m2)
{
	if (m2.empty()) {
		std::cout << "log with 0 base!\n";
		throw "f50d7f56-8ca7-4efe-9b23-3c9cde40d069";
	}
	MonRaw r;
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

PolyRaw reduce(const PolyRaw& poly, const PolysRaw& rels)
{
	PolyRaw result;
	PolyRaw poly1(poly);
	auto pbegin = poly1.begin(); auto pend = poly1.end();
	while (pbegin != pend) {
		PolysRaw::const_iterator p_rel = rels.begin();
		for (; p_rel != rels.end(); ++p_rel)
			if (divides((*p_rel)[0], *pbegin))
				break;
		if (p_rel == rels.end())
			result.push_back(std::move(*pbegin++));
		else {
			const PolyRaw& rel = *p_rel;
			MonRaw q = div(*pbegin, rel[0]);
			PolyRaw rel1 = mul(rel, q); // TODO: improve
			PolyRaw poly2;
			std::set_symmetric_difference(pbegin, pend, rel1.begin(), rel1.end(),
				std::back_inserter(poly2), cmp_mons);
			poly1 = std::move(poly2);
			pbegin = poly1.begin(); pend = poly1.end();
		}
	}
	return result;
}

MonRaw gcd(const MonRaw& m1, const MonRaw& m2)
{
	MonRaw result;
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

bool gcd_nonzero(const MonRaw& m1, const MonRaw& m2)
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

MonRaw lcm(const MonRaw& m1, const MonRaw& m2)
{
	MonRaw result;
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
	PolyRaw poly;
	int deg;
};

inline bool compare(const rel_heap_t& s1, const rel_heap_t& s2) { return s1.deg > s2.deg; }

void add_rel(PolysRaw& rels, const PolyRaw& rel, const std::vector<int> gen_degs)
{
	if (rel.empty())
		return;
	std::vector<rel_heap_t> heap = { rel_heap_t{rel, deg(rel, gen_degs)} };
	while (!heap.empty()) {
		std::pop_heap(heap.begin(), heap.end(), compare);
		rel_heap_t heap_ele = std::move(heap.back());
		heap.pop_back();

		PolyRaw r = reduce(heap_ele.poly, rels);
		if (!r.empty()) {
			for (PolyRaw& r1 : rels) {
				if (gcd_nonzero(r[0], r1[0])) {
					if (divides(r[0], r1[0]))
						r1.clear();
					else {
						MonRaw mlcm = lcm(r[0], r1[0]);
						MonRaw q = div(mlcm, r[0]);
						MonRaw q1 = div(mlcm, r1[0]);
						PolyRaw new_rel = add(mul(r, q), mul(r1, q1));

						heap.push_back(rel_heap_t{ new_rel, deg(mlcm, gen_degs) });
						std::push_heap(heap.begin(), heap.end(), [](const rel_heap_t& s1, const rel_heap_t& s2) { return s1.deg > s2.deg; });
					}
				}
			} 
			rels.erase(std::remove_if(rels.begin(), rels.end(), [](const PolyRaw& p) {return p.empty(); }), rels.end());
			rels.push_back(std::move(r));
		}
	}
}