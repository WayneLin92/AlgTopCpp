#include "mymath.h"

std::ostream& operator<<(std::ostream &sout, Mon s)
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

std::ostream& operator<<(std::ostream &sout, Poly s)
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

Poly operator^(Poly lhs, Poly rhs)
{
	Poly result;
	std::set_symmetric_difference(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::inserter(result, result.begin()));
	return result;
}