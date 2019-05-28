#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>

typedef std::vector<int> Mon;  // ordered
typedef std::vector<Mon> Poly;  // ordered
typedef std::vector<Poly> Relations;
struct log_mon_t {
	int q;
	Mon r;
};

// algmod2.cpp
inline void sort(Poly& p) { std::sort(p.begin(), p.end(), std::greater<Mon>()); };

std::ostream& operator<<(std::ostream& sout, const Mon& s);
std::ostream& operator<<(std::ostream& sout, const Poly& s);
Mon operator*(const Mon& m1, const Mon& m2);
Mon operator/(const Mon& m1, const Mon& m2);
Mon& operator*=(Mon& m1, const Mon& m2);
Poly operator+(const Poly& poly1, const Poly& poly2);
Poly operator+(const Poly& poly1, const Poly& poly2);
Poly& operator+=(Poly& poly1, const Poly& poly2);
Poly operator*(const Poly& poly, const Mon& mon);

inline int deg(const Mon& m) { return std::accumulate(m.begin(), m.end(), 0); };
inline int deg(const Poly& p) { return p.size() ? deg(p[0]) : -1; };
bool divides(const Mon& m1, const Mon& m2);
Mon gcd(const Mon& m1, const Mon& m2);
Mon lcm(const Mon& m1, const Mon& m2);
Mon pow(const Mon& m, int e);
log_mon_t log_mon(const Mon& m1, const Mon& m2);

// groebner.cpp
Poly& simplify(const Relations& rels, Poly& poly);
void add_rel(Relations& rels, const Poly& rel);
void test();
