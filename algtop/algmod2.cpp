#include "mymath.h"

// operator overloading
std::ostream& operator<<(std::ostream& sout, const Mon& s)
{
	sout << "(";
	for (auto i = s.begin(); i != s.end(); ++i)
	{
		if (i != s.begin())
			sout << ", ";
		sout << *i;
	}
	sout << ")";
	return sout;
}

std::ostream& operator<<(std::ostream& sout, const Poly& s)
{
	sout << "{";
	for (auto i = s.begin(); i != s.end(); ++i)
	{
		if (i != s.begin())
			sout << ", ";
		sout << *i;
	}
	sout << "}";
	return sout;
}

std::ostream& operator<<(std::ostream& sout, const Relations& s)
{
	sout << "[";
	for (auto i = s.begin(); i != s.end(); ++i)
	{
		if (i != s.begin())
			sout << ", ";
		sout << *i;
	}
	sout << "]";
	return sout;
}

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
	int i;
	for (i = result.size() - 1; i >= 0 && result[i] == 0; i--);
	result.resize(i + 1);
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

Poly operator*(const Poly& poly, const Mon& mon)
{
	Poly result = poly;
	for (Mon& m : result)
		m *= mon;
	return result;
}

// other functions
bool divides(const Mon& m1, const Mon& m2)
{
	if (m1.size() <= m2.size()) {
		for (unsigned int i = 0; i < m1.size(); i++)
			if (m1[i] > m2[i])
				return false;
		return true;
	}
	return false;
}

Mon gcd(const Mon& m1, const Mon& m2)
{
	Mon result;
	unsigned int len = std::min(m1.size(), m2.size());
	for (unsigned int i = 0; i < len; i++)
		result.push_back(std::min(m1[i], m2[i]));
	return result;
}

Mon lcm(const Mon& m1, const Mon& m2)
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

log_mon_t log_mon(const Mon& m1, const Mon& m2) // assume m2 divides m1 and m2 != ()
{
	Mon r = m1;
	int q = -1;
	for (unsigned int i = 0; i < m2.size(); i++) {
		if (m2[i] > 0) {
			int q1 = r[i] / m2[i];
			if (q > q1 || q == -1) q = q1;
		}
	}
	for (unsigned int i = 0; i < m2.size(); i++)
		r[i] -= m2[i] * q;
	int i;
	for (i = r.size() - 1; i >= 0 && r[i] == 0; i--);
	r.resize(i + 1);
	log_mon_t result = { q, r };
	return result;
}

Mon pow(const Mon& m, int e)
{
	Mon result;
	for (int i : m)
		result.push_back(i * e);
	return result;
}
