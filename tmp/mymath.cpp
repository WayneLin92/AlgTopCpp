#include "mymath.h"
#include <iterator>
#include <numeric>
#include <algorithm>

std::ostream& operator<<(std::ostream& sout, Mon& s)
{
	Mon::iterator i;
	sout << "(";
	for (i = s.begin(); i != s.end(); ++i)
	{
		if (i != s.begin())
			sout << ", ";
		sout << *i;
	}
	sout << ")";
	return sout;
}

std::ostream& operator<<(std::ostream& sout, Poly& s)
{
	Poly::iterator i;
	sout << "{";
	for (i = s.begin(); i != s.end(); ++i)
	{
		if (i != s.begin())
			sout << ", ";
		sout << *i;
	}
	sout << "}";
	return sout;
}

std::ostream& operator<<(std::ostream& sout, Relations& s)
{
	Relations::iterator i;
	sout << "[";
	for (i = s.begin(); i != s.end(); ++i)
	{
		if (i != s.begin())
			sout << ", ";
		sout << *i;
	}
	sout << "]";
	return sout;
}

Mon operator*(const Mon& m1, const Mon& m2)
{
	Mon result;
	unsigned int i;
	if (m1.size() <= m2.size()) {
		for (i = 0; i < m1.size(); i++)
			result.push_back(m1[i] + m2[i]);
		result.insert(result.end(), m2.begin() + m1.size(), m2.end());
	}
	else {
		for (i = 0; i < m2.size(); i++)
			result.push_back(m1[i] + m2[i]);
		result.insert(result.end(), m1.begin() + m2.size(), m1.end());
	}
	return result;
}

Mon operator/(const Mon& m1, const Mon& m2) // assume m2 divides m1
{
	Mon result = m1;
	for (unsigned int i = 0; i < m2.size(); i++) {
		result[i] -= m2[i];
	}
	while (result.size() && !result.back())
		result.pop_back();
	return result;
}

Mon& operator*=(Mon& m1, const Mon& m2)
{
	unsigned int i;
	if (m1.size() <= m2.size()) {
		for (i = 0; i < m1.size(); i++)
			m1[i] += m2[i];
		m1.insert(m1.end(), m2.begin() + m1.size(), m2.end());
	}
	else {
		for (i = 0; i < m2.size(); i++)
			m1[i] += m2[i];
	}
	return m1;
}



Poly operator^(Poly& lhs, Poly& rhs)
{
	Poly result;
	std::set_symmetric_difference(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::inserter(result, result.begin()));
	return result;
}

int deg(Mon& m) { return std::accumulate(m.begin(), m.end(), 0); }

int deg(Poly& p) { return p.size() ? deg(p[0]) : -1; }

bool divides(Mon& m1, Mon& m2)
{
	if (m1.size() <= m2.size()) {
		for (unsigned int i = 0; i < m1.size(); i++)
			if (m1[i] > m2[i])
				return false;
		return true;
	}
	return false;
}

Mon lcm_mon(Mon& m1, Mon& m2)
{
	Mon result;
	unsigned int i;
	if (m1.size() <= m2.size()) {
		for (i = 0; i < m1.size(); i++)
			result.push_back(std::max(m1[i], m2[i]));
		result.insert(result.end(), m2.begin() + m1.size(), m2.end());
	}
	else {
		for (i = 0; i < m2.size(); i++)
			result.push_back(std::max(m1[i], m2[i]));
		result.insert(result.end(), m1.begin() + m2.size(), m1.end());
	}
	return result;
}

bool compare(Poly& poly1, Poly& poly2) { return deg(poly1) < deg(poly2); }

void simplify(Relations& rels, Poly& poly)
{

}

void add_rel(Relations& rels, Poly& rel)
{
	if (rel.empty())
		return;
	Relations heap = { rel };
	while (heap.size() > 0) {
		std::pop_heap(heap.begin(), heap.end(), compare);
		Poly r = std::move(heap.back());
		heap.pop_back();
		if (r.size() > 0) {
			for (Poly& r1 : rels) {
				int len = std::min(r[0].size(), r1[0].size());
				bool gcdnonzero = false;
				for (int i = 0; i < len; i++) {
					if (r[0][i] || r1[0][i]) {
						gcdnonzero = true;
						break;
					}
				}
				if (gcdnonzero) {
					if (divides(r[0], r1[0]))
						r1.clear();
					else {
						Mon lcm = lcm_mon(r[0], r1[0]);
						Mon q = lcm / r[0];
						Mon q1 = lcm / r1[0];
					}
				}
				else {

				}
			}
		}
	}
}

void test()
{
	Mon m1 = { 1, 9, 3, 8, 2, 4, 5 };
	Mon m2 = { 2, 9, 3, 8, 2, 4, 5, 1 };
	lprod_mon(m1, m2);
	std::cout << m1 << std::endl;
}

// 91