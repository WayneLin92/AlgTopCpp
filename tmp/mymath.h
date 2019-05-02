#pragma once

#include <vector>
#include <iterator>
#include <algorithm>

typedef std::vector<int> Mon;
typedef std::vector<Mon> Poly;
struct Relation {
	Mon lead;
	Poly rest;
};
typedef std::vector<Relation> Relations;

std::ostream& operator<<(std::ostream &sout, Mon s);
std::ostream& operator<<(std::ostream &sout, Poly s);
Poly operator^(Poly lhs, Poly rhs);

void add_rel(Relations& rels, Relation& rel)
{
	
}
