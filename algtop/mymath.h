#pragma once
#include "myparser.h"
#include <algorithm>
#include <numeric>

typedef std::vector<int> Mon;  // ordered
typedef std::vector<Mon> Poly;  // ordered
typedef std::vector<Poly> Relations;

inline void sort(Poly& p) { std::sort(p.begin(), p.end(), std::greater<Mon>()); }
inline void load_mon(Mon& mon, std::istream& sin)         { load_vector(mon, sin, "(", ",", ")"); }
inline void load_poly(Poly& poly, std::istream& sin)      { load_vector(poly, sin, "{", ",", "}", load_mon); }
inline void load_rels(Relations& rels, std::istream& sin) { load_vector(rels, sin, "[", ",", "]", load_poly); }
inline void dump_mon(const Mon& mon, std::ostream& sout)              { dump_vector(mon, sout, "(", ", ", ")"); }
inline void dump_poly(const Poly& poly, std::ostream& sout)           { dump_vector(poly, sout, "{", ", ", "}", dump_mon); }
inline void dump_relations(const Relations& rels, std::ostream& sout) { dump_vector(rels, sout, "[", ", ", "]", dump_poly); }

inline std::istream& operator>>(std::istream& sin, Mon& mon)        { load_mon(mon, sin); return sin; }
inline std::istream& operator>>(std::istream& sin, Poly& poly)      { load_poly(poly, sin); return sin; }
inline std::istream& operator>>(std::istream& sin, Relations& rels) { load_rels(rels, sin); return sin; }
inline std::ostream& operator<<(std::ostream& sout, const Mon& mon)        { dump_mon(mon, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Poly& poly)      { dump_poly(poly, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Relations& rels) { dump_relations(rels, sout); return sout; }

// algmod2.cpp
Mon operator*(const Mon& m1, const Mon& m2);
Mon& operator*=(Mon& m1, const Mon& m2);
Mon operator/(const Mon& m1, const Mon& m2);
Poly operator*(const Poly& poly, const Mon& mon);
Poly operator+(const Poly& poly1, const Poly& poly2);
Poly operator+(const Poly& poly1, const Poly& poly2);
Poly& operator+=(Poly& poly1, const Poly& poly2);

inline int deg(const Mon& m) { return std::accumulate(m.begin(), m.end(), 0); };
inline int deg(const Poly& p) { return p.size() ? deg(p[0]) : -1; };
int deg(const Mon& mon, const std::vector<int>& gen_degs);
inline int deg(const Poly& p, const std::vector<int>& gen_degs) { return p.size() ? deg(p[0], gen_degs) : -1; };
bool divides(const Mon& m1, const Mon& m2);
Mon gcd(const Mon& m1, const Mon& m2);
Mon lcm(const Mon& m1, const Mon& m2);
Mon pow(const Mon& m, int e);
struct log_mon_t { int q; Mon r; };
log_mon_t log_mon(const Mon& m1, const Mon& m2);

// groebner.cpp
Poly& simplify(Poly& poly, const Relations& rels);
void add_rel(const Poly& rel, Relations& rels, std::vector<int>* pgen_degs = NULL);


void test();
