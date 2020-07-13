#ifndef MYMATH_H
#define MYMATH_H

#include "myparser.h"
#include <algorithm>
#include <numeric>

/* sparse monomial i1, e1, i2, e2, ... */
typedef std::vector<int> Mon;
typedef std::vector<Mon> Poly;  // ordered
typedef std::vector<Poly> Polys;

inline void sort(Poly& p) { std::sort(p.begin(), p.end(), std::greater<Mon>()); }
inline void load_mon(Mon& mon, std::istream& sin)         { load_vector(mon, sin, "(", ",", ")"); } // TODO: check if empty string is properly handled
inline void load_poly(Poly& poly, std::istream& sin)      { load_vector(poly, sin, "{", ",", "}", load_mon); }
inline void load_polys(Polys& polys, std::istream& sin)   { load_vector(polys, sin, "[", ",", "]", load_poly); }
inline void dump_mon(const Mon& mon, std::ostream& sout)              { dump_vector(mon, sout, "(", ", ", ")"); }
inline void dump_poly(const Poly& poly, std::ostream& sout)           { dump_vector(poly, sout, "{", ", ", "}", dump_mon); }
inline void dump_polys(const Polys& polys, std::ostream& sout)        { dump_vector(polys, sout, "[", ", ", "]", dump_poly); }

inline std::istream& operator>>(std::istream& sin, Mon& mon)        { load_mon(mon, sin); return sin; }
inline std::istream& operator>>(std::istream& sin, Poly& poly)      { load_poly(poly, sin); return sin; }
inline std::istream& operator>>(std::istream& sin, Polys& polys)    { load_polys(polys, sin); return sin; }
inline std::ostream& operator<<(std::ostream& sout, const Mon& mon)        { dump_mon(mon, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Poly& poly)      { dump_poly(poly, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Polys& polys)    { dump_polys(polys, sout); return sout; }

/* algmod2.cpp */
std::vector<int> mul_mons(const std::vector<int>& mon1, const std::vector<int>& mon2);
inline Mon operator*(const Mon& m1, const Mon& m2) { return mul_mons(m1, m2); };
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
bool divides(const std::vector<int>& m1, const std::vector<int>& m2);
/* cmp_mons(m1, m2) returns true if e1 < e2 in lexicographical order where e1, e2 are arrays of exponents. */
bool cmp_mons(const std::vector<int>& m1, const std::vector<int>& m2);
Mon gcd(const Mon& m1, const Mon& m2);
Mon lcm(const Mon& m1, const Mon& m2);
Mon pow(const Mon& m, int e);
struct log_mon_t { int q; Mon r; };
log_mon_t log_mon(const Mon& m1, const Mon& m2);

// groebner.cpp
class GroebnerBasis
{
public:
	GroebnerBasis() {};
	~GroebnerBasis() {};
private:
	Polys::const_iterator find_leading_divisor(const Mon mon) const {
		// Return the pointer of r in m_rels if r[0] | mon otherwise return m_rels::end();
		Polys::const_iterator p_rel = m_rels.begin();
		for (; p_rel != m_rels.end(); ++p_rel)
			if (divides((*p_rel)[0], mon))
				return p_rel;
		return p_rel;
	}
public:
	Poly& simplify(Poly& poly) const;
	void add_rel(const Poly& rel, const std::vector<int>* pgen_degs = nullptr);
public:
	Polys m_rels;
};


void test();

#endif /* MYMATH_H */
