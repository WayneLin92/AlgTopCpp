#pragma once
#include <vector>
#include <iostream>

typedef std::vector<int> Mon;  // ordered
typedef std::vector<Mon> Poly;  // ordered
typedef std::vector<Poly> Relations;

std::ostream& operator<<(std::ostream& sout, Mon& s);
std::ostream& operator<<(std::ostream& sout, Poly& s);
Poly operator^(Poly& lhs, Poly& rhs);

int deg(Mon& m);
int deg(Poly& p);
void simplify(Relations& rels, Poly& poly);
void add_rel(Relations& rels, Poly& rel);

void test();
