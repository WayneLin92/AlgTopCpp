#ifndef MYMATH_H
#define MYMATH_H

#include "myparser.h"
#include <algorithm>
#include <numeric>

/* sparse monomial i1, e1, i2, e2, ... */
typedef std::vector<int> MonRaw;
typedef std::vector<MonRaw> PolyRaw;
typedef std::vector<PolyRaw> PolysRaw;

struct Mon {
	MonRaw data;
	Mon(MonRaw data_) : data(data_) {};
};

struct Poly {
	PolyRaw data;
	Poly(PolyRaw data_) : data(data_) {};
	Poly(Mon mon) : data({ mon.data }) {};
};

inline void sort(PolyRaw& p) { std::sort(p.begin(), p.end(), std::greater<MonRaw>()); }
inline void load_mon(MonRaw& mon, std::istream& sin)         { load_vector(mon, sin, "(", ",", ")"); } // TODO: check if empty string is properly handled
inline void load_poly(PolyRaw& poly, std::istream& sin)      { load_vector(poly, sin, "{", ",", "}", load_mon); }
inline void load_polys(PolysRaw& polys, std::istream& sin)   { load_vector(polys, sin, "[", ",", "]", load_poly); }
inline void dump_mon(const MonRaw& mon, std::ostream& sout)              { dump_vector(mon, sout, "(", ", ", ")"); }
inline void dump_poly(const PolyRaw& poly, std::ostream& sout)           { dump_vector(poly, sout, "{", ", ", "}", dump_mon); }
inline void dump_polys(const PolysRaw& polys, std::ostream& sout)        { dump_vector(polys, sout, "[", ", ", "]", dump_poly); }

inline std::istream& operator>>(std::istream& sin, MonRaw& mon)        { load_mon(mon, sin); return sin; }
inline std::istream& operator>>(std::istream& sin, PolyRaw& poly)      { load_poly(poly, sin); return sin; }
inline std::istream& operator>>(std::istream& sin, PolysRaw& polys)    { load_polys(polys, sin); return sin; }
inline std::ostream& operator<<(std::ostream& sout, const MonRaw& mon)        { dump_mon(mon, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const PolyRaw& poly)      { dump_poly(poly, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const PolysRaw& polys)    { dump_polys(polys, sout); return sout; }

inline std::ostream& operator<<(std::ostream& sout, const Poly& poly) { dump_poly(poly.data, sout); return sout; }

/* algmod2.cpp */
MonRaw mul(const MonRaw& mon1, const MonRaw& mon2);
MonRaw div(const MonRaw& mon1, const MonRaw& mon2);
PolyRaw add(const PolyRaw& poly1, const PolyRaw& poly2);
PolyRaw mul(const PolyRaw& poly, const MonRaw& mon);
PolyRaw mul(const PolyRaw& poly1, const PolyRaw& poly2);
inline Mon operator*(const Mon& m1, const Mon& m2) { return mul(m1.data, m2.data); }
inline Poly operator+(const Poly& p1, const Poly& p2) { return add(p1.data, p2.data); }
inline Poly operator+=(Poly& p1, const Poly& p2) { return p1 = add(p1.data, p2.data); }
inline Poly operator*(const Poly& p, const Mon& m) { return mul(p.data, m.data); }
inline Poly operator*(const Mon& m, const Poly& p) { return mul(p.data, m.data); }
inline Poly operator*(const Poly& p1, const Poly& p2) { return mul(p1.data, p2.data); }

bool divides(const MonRaw& m1, const MonRaw& m2);
MonRaw pow(const MonRaw& m, int e);
struct log_result { int q; MonRaw r; };
/* m1 = m2^q * r */
log_result log(const MonRaw& m1, const MonRaw& m2);

/* reduce `poly` by groebner basis `rels` */
PolyRaw reduce(const PolyRaw& poly, const PolysRaw& rels);
/* Add new relation `rel` to groebner basis `rels` */
void add_rel(PolysRaw& rels, const PolyRaw& rel, const std::vector<int> gen_degs);

inline int deg(const MonRaw& mon) { int result = 0; for (size_t i = 0; i < mon.size(); i += 2) result += mon[i + 1]; return result; };
inline int deg(const PolyRaw& p) { return p.size() ? deg(p[0]) : -1; };
inline int deg(const MonRaw& mon, const std::vector<int>& gen_degs) {
	int result = 0; for (size_t i = 0; i < mon.size(); i += 2) result += gen_degs[mon[i]] * mon[i + 1]; return result;
};
inline int deg(const PolyRaw& p, const std::vector<int>& gen_degs) { return p.size() ? deg(p[0], gen_degs) : -1; }
/* cmp_mons(m1, m2) returns true if e1 < e2 in lexicographical order where e1, e2 are arrays of exponents. */
bool cmp_mons(const MonRaw& m1, const MonRaw& m2);
MonRaw gcd(const MonRaw& m1, const MonRaw& m2);
MonRaw lcm(const MonRaw& m1, const MonRaw& m2);


int main_test();

#endif /* MYMATH_H */
