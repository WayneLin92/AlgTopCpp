#include "main.h"
#include <sstream>

/*********** FUNCTIONS **********/

void dump_MonPow(const GenPow& p, std::ostream& sout)
{
	sout << "x_";
	if (0 <= p.gen && p.gen < 10)
		sout << p.gen;
	else
		sout << '{' << p.gen << '}';
	if (p.exp > 1) {
		sout << '^';
		if (0 <= p.exp && p.exp < 10)
			sout << p.exp;
		else
			sout << '{' << p.exp << '}';
	}
};

std::string array_to_str(array::const_iterator pbegin, array::const_iterator pend)
{
	std::stringstream ss;
	for (auto p = pbegin; p < pend; p++) {
		ss << *p;
		if (p + 1 < pend)
			ss << ",";
	}
	return ss.str();
}

array str_to_array(const char* str_mon)
{
	array result;
	if (str_mon[0] == '\0')
		return array();
	std::stringstream ss(str_mon);
	while (ss.good()) {
		int i;
		ss >> i;
		result.push_back(i);
		if (ss.peek() == ',')
			ss.ignore();
	}
	return result;
}

std::string Mon_to_str(MonInd pbegin, MonInd pend)
{
	std::stringstream ss;
	for (auto p = pbegin; p < pend; p++) {
		ss << p->gen << ',' << p->exp;
		if (p + 1 < pend)
			ss << ",";
	}
	return ss.str();
}

Mon str_to_Mon(const char* str_mon)
{
	Mon result;
	if (str_mon[0] == '\0')
		return {};
	std::stringstream ss(str_mon);
	while (ss.good()) {
		int gen, exp;
		ss >> gen >> "," >> exp;
		result.emplace_back(gen, exp);
		if (ss.peek() == ',')
			ss.ignore();
	}
	return result;
}

std::string Poly_to_str(Poly::const_iterator pbegin, Poly::const_iterator pend) /* Warning: assume the algebra is connected */
{
	std::stringstream ss;
	for (auto pMon = pbegin; pMon < pend; ++pMon) {
		for (auto p = pMon->begin(); p != pMon->end(); ++p) {
			ss << p->gen << ',' << p->exp;
			if (p + 1 != pMon->end())
				ss << ",";
		}
		if (pMon + 1 != pend)
			ss << ";";
	}
	return ss.str();
}

Poly str_to_Poly(const char* str_poly) /* Warning: assume the algebra is connected */
{
	Poly result;
	if (str_poly[0] == '\0')
		return {};
	std::stringstream ss(str_poly);
	while (ss.good()) {
		int gen, exp;
		ss >> gen >> "," >> exp;
		if (result.empty())
			result.emplace_back();
		result.back().emplace_back(gen, exp);
		if (ss.peek() == ',')
			ss.ignore();
		else if (ss.peek() == ';') {
			ss.ignore();
			result.emplace_back();
		}
	}
	return result;
}

Poly indices_to_Poly(const array& indices, const Poly& basis)
{
	Poly result;
	for (int i : indices)
		result.push_back(basis[i]);
	return result;
}

array Poly_to_indices(const Poly& poly, const Poly& basis)
{
	array result;
	for (const Mon& mon : poly)
		result.push_back(get_index(basis, mon));
	return result;
}