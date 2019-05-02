#include "mymath.h"
#include <iostream>  // for test

std::ostream& operator<<(std::ostream& sout, Mon& s)
{
	Mon::iterator i;
	sout << "(";
	for (i = s.begin(); i != s.end(); ++i)
	{
		if (i != s.begin())
		{
			sout << ", ";
		}
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
		{
			sout << ", ";
		}
		sout << *i;
	}
	sout << "}";
	return sout;
}

std::ostream& operator<<(std::ostream& sout, Relations& s)
{
	Relations::iterator i;
	sout << "{";
	for (i = s.begin(); i != s.end(); ++i)
	{
		if (i != s.begin())
		{
			sout << ", ";
		}
		sout << *i;
	}
	sout << "}";
	return sout;
}

Poly operator^(Poly& lhs, Poly& rhs)
{
	Poly result;
	std::set_symmetric_difference(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::inserter(result, result.begin()));
	return result;
}

int deg(Mon &m)
{
	return std::accumulate(m.begin(), m.end(), 0);
}

int deg(Poly& p)
{
	if (p.size() == 0)
		return -1;
	else
		return deg(p[0]);
}

bool compare(Poly& poly1, Poly& poly2)
{
	return deg(poly1) < deg(poly2);
}

void Simplify(Relations& rels, Poly& poly)
{

}

void AddRel(Relations& rels, Poly& rel)
{
	if (rel.size() == 0)
		return;
	Relations heap = { rel };
	Relations::iterator h_end = heap.end();
	while (heap.size() > 0) {
		std::pop_heap(heap.begin(), h_end, compare);
		h_end--;
		if (h_end->size() > 0) {
			
		}
	}
}

void test()
{
	Poly s1 = { {1, 2}, {0, 2, 1}, {1, 1, 1} };
	Poly s2 = { {3, 2, 1}, {0, 0, 6}, {4, 2} };
	std::sort(s1.begin(), s1.end());
	std::sort(s2.begin(), s2.end());
	Relations rels = { s1, s2 };
	std::make_heap(rels.begin(), rels.end(), compare);
	std::cout << rels << std::endl;
}