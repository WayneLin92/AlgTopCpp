#include "database.h"
#include <sstream>
//#include <chrono>

/*********** FUNCTIONS **********/

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

std::string array2d_to_str(array2d::const_iterator pbegin, array2d::const_iterator pend) /* Warning: assume inner arrays are all nontrivial */
{
	std::stringstream ss;
	for (auto pMon = pbegin; pMon < pend; pMon++) {
		for (auto p = pMon->begin(); p < pMon->end(); p++) {
			ss << *p;
			if (p + 1 < pMon->end())
				ss << ",";
		}
		if (pMon + 1 < pend)
			ss << ";";
	}
	return ss.str();
}

array2d str_to_array2d(const char* str_poly) /* Warning: assume inner arrays are all nontrivial */
{
	array2d result;
	if (str_poly[0] == '\0')
		return array2d();
	std::stringstream ss(str_poly);
	while (ss.good()) {
		int i;
		ss >> i;
		if (result.empty())
			result.push_back(array());
		result.back().push_back(i);
		if (ss.peek() == ',')
			ss.ignore();
		else if (ss.peek() == ';') {
			ss.ignore();
			result.push_back(array());
		}
	}
	return result;
}