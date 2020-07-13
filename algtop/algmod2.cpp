#include "mymath.h"
#include <sstream> ////
#include <string> ////

// operator overloading
Poly operator+(const Poly& poly1, const Poly& poly2)
{
	Poly result;
	std::set_symmetric_difference(poly1.begin(), poly1.end(), poly2.begin(), poly2.end(),
		std::back_inserter(result), std::greater<Mon>());
	return result;
}

Poly& operator+=(Poly& poly1, const Poly& poly2)
{
	Poly result;
	std::set_symmetric_difference(poly1.begin(), poly1.end(), poly2.begin(), poly2.end(),
		std::back_inserter(result), std::greater<Mon>());
	poly1 = std::move(result);
	return poly1;
}

std::vector<int> mul_mons(const std::vector<int>& m1, const std::vector<int>& m2) //
{
	Mon result;
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

Mon operator/(const Mon& m1, const Mon& m2) // assume m2 divides m1
{
	Mon result = m1;
	for (size_t i = 0; i < m2.size(); ++i)
		result[i] -= m2[i];
	size_t i;
	for (i = result.size(); i > 0 && result[i-1] == 0; --i);
	result.resize(i);
	return result;
}

Mon& operator*=(Mon& m1, const Mon& m2)
{
	if (m1.size() <= m2.size()) {
		for (size_t i = 0; i < m1.size(); ++i)
			m1[i] += m2[i];
		m1.insert(m1.end(), m2.begin() + m1.size(), m2.end());
	}
	else {
		for (size_t i = 0; i < m2.size(); ++i)
			m1[i] += m2[i];
	}
	return m1;
}

Poly operator*(const Poly& poly, const Mon& mon)
{
	Poly result = poly;
	for (Mon& m : result)
		m *= mon;
	return result;
}

int deg(const Mon& mon, const std::vector<int>& gen_degs)
{
	int result = 0;
	for (size_t i = 0; i < mon.size(); ++i)
		result += mon[i] * gen_degs[i];
	return result;
}

// other functions
bool divides(const std::vector<int>& m1, const std::vector<int>& m2)
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

bool cmp_mons(const std::vector<int>& m1, const std::vector<int>& m2)
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

Mon gcd(const Mon& m1, const Mon& m2)
{
	Mon result;
	size_t len = std::min(m1.size(), m2.size());
	for (size_t i = 0; i < len; ++i)
		result.push_back(std::min(m1[i], m2[i]));
	return result;
}

Mon lcm(const Mon& m1, const Mon& m2)
{
	Mon result;
	unsigned int i;
	if (m1.size() <= m2.size()) {
		for (i = 0; i < m1.size(); ++i)
			result.push_back(std::max(m1[i], m2[i]));
		result.insert(result.end(), m2.begin() + m1.size(), m2.end());
	}
	else {
		for (i = 0; i < m2.size(); ++i)
			result.push_back(std::max(m1[i], m2[i]));
		result.insert(result.end(), m1.begin() + m2.size(), m1.end());
	}
	return result;
}

log_mon_t log_mon(const Mon& m1, const Mon& m2) // assume m2 divides m1 and m2 != ()
{
	Mon r = m1;
	int q = -1;
	for (unsigned int i = 0; i < m2.size(); ++i) {
		if (m2[i] > 0) {
			int q1 = r[i] / m2[i];
			if (q > q1 || q == -1) q = q1;
		}
	}
	for (unsigned int i = 0; i < m2.size(); ++i)
		r[i] -= m2[i] * q;
	size_t i;
	for (i = r.size(); i > 0 && r[i-1] == 0; --i);
	r.resize(i);
	return log_mon_t{ q, r };
}

Mon pow(const Mon& m, int e)
{
	Mon result;
	for (int i : m)
		result.push_back(i * e);
	return result;
}
// 255, 205, 147, 143