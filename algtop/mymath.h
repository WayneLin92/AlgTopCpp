#ifndef MYMATH_H
#define MYMATH_H
// TODO: greobner basis array4d version
// TODO: monomial std::vector<MonPow> version

#include "myparser.h"
#include <queue>
#include <algorithm>
#include <numeric>

/********** STRUCTS AND CLASSES **********/
using array = std::vector<int>;
using array2d = std::vector<array>;
using array3d = std::vector<array2d>;
using array4d = std::vector<array3d>;
using arrayInd = array::const_iterator;

struct MonPow {
	int gen, exp;
	MonPow(int g, int e) : gen(g), exp(e) {}
	bool operator<(const MonPow& rhs) const { return gen > rhs.gen || (gen == rhs.gen && exp < rhs.exp); }
};
using Mon = std::vector<MonPow>; // only monomials in the same degree are comparable
using Poly = std::vector<Mon>;
using Poly1d = std::vector<Poly>;
using Poly2d = std::vector<Poly1d>;
using Mon1d = std::vector<Mon>;
using Mon2d = std::vector<Mon1d>;
using MonInd = Mon::const_iterator;

struct Deg
{
	int s, t, v;
	Deg(int s_, int t_, int v_) : s(s_), t(t_), v(v_) {};
	bool operator<(const Deg& rhs) const {
		if (t < rhs.t)
			return true;
		else if (t == rhs.t)
			if (s < rhs.s)
				return true;
			else if (s == rhs.s)
				if (v < rhs.v)
					return true;
		return false;
	};
	Deg operator+(const Deg& rhs) const {
		return Deg({ s + rhs.s, t + rhs.t, v + rhs.v });
	};
	Deg operator-(const Deg& rhs) const {
		return Deg({ s - rhs.s, t - rhs.t, v - rhs.v });
	};
	Deg& operator+=(const Deg& rhs) {
		s += rhs.s; t += rhs.t; v += rhs.v;
		return *this;
	};
	Deg& operator-=(const Deg& rhs) {
		s -= rhs.s; t -= rhs.t; v -= rhs.v;
		return *this;
	};
	bool operator==(const Deg& rhs) const {
		return v == rhs.v && s == rhs.s && t == rhs.t;
	};
	bool operator!=(const Deg& rhs) const {
		return t != rhs.t || s != rhs.s || v != rhs.v;
	};
	Deg operator*(int rhs) const {
		return Deg({ s * rhs, t * rhs, v * rhs });
	};
};

struct PolyWithT {
	Poly poly;
	int deg;
	bool operator<(const PolyWithT& rhs) const { return deg > rhs.deg; }
};

using RelHeap = std::priority_queue<PolyWithT>;

/********** FUNCTIONS **********/
// Move to a header for python in the future

inline void load_array(array& mon, std::istream& sin) { load_vector(mon, sin, "(", ",", ")"); } // TODO: check if empty string is properly handled
inline void load_array2d(array2d& poly, std::istream& sin) { load_vector(poly, sin, "{", ",", "}", load_array); }
inline void load_array3d(array3d& polys, std::istream& sin) { load_vector(polys, sin, "[", ",", "]", load_array2d); }
inline void dump_array(const array& mon, std::ostream& sout) { dump_vector(mon, sout, "(", ", ", ")"); }
inline void dump_array2d(const array2d& poly, std::ostream& sout) { dump_vector(poly, sout, "{", ", ", "}", dump_array); }
inline void dump_array3d(const array3d& polys, std::ostream& sout) { dump_vector(polys, sout, "[", ", ", "]", dump_array2d); }
inline void dump_array4d(const array4d& a, std::ostream& sout) { dump_vector(a, sout, "[", ", ", "]", dump_array3d); }

void dump_MonPow(const MonPow& p, std::ostream& sout);
inline void dump_Mon(const Mon& mon, std::ostream& sout) { dump_vector(mon, sout, "", "", "", dump_MonPow); }
inline void dump_Poly(const Poly& poly, std::ostream& sout) { if (poly.empty()) sout << '0'; else dump_vector(poly, sout, "", "+", "", dump_Mon); }
inline void dump_Poly1d(const Poly1d& polys, std::ostream& sout) { dump_vector(polys, sout, "(", ", ", ")", dump_Poly); }
inline void dump_Poly2d(const Poly2d& polys, std::ostream& sout) { dump_vector(polys, sout, "[", ", ", "]", dump_Poly1d); }

inline std::istream& operator>>(std::istream& sin, array& mon) { load_array(mon, sin); return sin; }
inline std::istream& operator>>(std::istream& sin, array2d& poly) { load_array2d(poly, sin); return sin; }
inline std::istream& operator>>(std::istream& sin, array3d& polys) { load_array3d(polys, sin); return sin; }

inline std::ostream& operator<<(std::ostream& sout, const array& mon) { dump_array(mon, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const array2d& poly) { dump_array2d(poly, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const array3d& polys) { dump_array3d(polys, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const array4d& a) { dump_array4d(a, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Deg& d) { sout << '(' << d.s << ',' << d.t << ',' << d.v << ')'; return sout; }
inline std::ostream& operator<<(std::ostream& sout, const MonPow& p) { dump_MonPow(p, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Mon& mon) { dump_Mon(mon, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Poly& poly) { dump_Poly(poly, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Poly1d& polys) { dump_Poly1d(polys, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const Poly2d& polys) { dump_Poly2d(polys, sout); return sout; }

inline void hash_combine(std::size_t& seed, int v) { seed ^= std::hash<int>{}(v)+0x9e3779b9 + (seed << 6) + (seed >> 2); }
inline size_t hash_range(const array& a)
{
	size_t seed = 0;
	for (int i : a)
		hash_combine(seed, i);
	return seed;
}

/*--------- algmod2.cpp ---------*/
/* Functions for Monomials and Polynomials
**
** A monomial m is an array of i1, e1, i2, e2, ...
** Monomials are ordered lexicographically by exponents after filling in zero exponents.
** A polynomials is an increasing sequence of monomials.
*/

Mon mul(const Mon& mon1, const Mon& mon2);
Mon div(const Mon& mon1, const Mon& mon2);
Poly add(const Poly& poly1, const Poly& poly2);
Poly mul(const Poly& poly, const Mon& mon);
inline Poly mul(const Mon& mon, const Poly& poly) { return mul(poly, mon); }
Poly mul(const Poly& poly1, const Poly& poly2);
Mon pow(const Mon& m, int e);
Poly pow(const Poly& poly, int n);
Poly pow(const Poly& poly, int n, const Poly1d& gb);
bool divides(const Mon& m1, const Mon& m2);
/* return the largest integer e where m1 = m2^e * r */
int log(const Mon& m1, const Mon& m2);

inline int get_deg(const Mon& mon) { int result = 0; for (MonInd p = mon.begin(); p != mon.end(); ++p) result += p->exp; return result; };
inline int get_deg(const Poly& poly) { return poly.size() ? get_deg(poly[0]) : -1; };
inline int get_deg(const Mon& mon, const array& gen_degs) { int result = 0; for (MonInd p = mon.begin(); p != mon.end(); ++p) result += gen_degs[p->gen] * p->exp; return result; };
inline Deg get_deg(const Mon& mon, const std::vector<Deg>& gen_degs) { Deg result({ 0, 0, 0 }); for (MonInd p = mon.begin(); p != mon.end(); ++p) result += gen_degs[p->gen] * p->exp; return result; };
inline int get_deg_t(const Mon& mon, const std::vector<Deg>& gen_degs) { int result = 0; for (MonInd p = mon.begin(); p != mon.end(); ++p) result += gen_degs[p->gen].t * p->exp; return result; };
inline int get_deg(const Poly& p, const array& gen_degs) { return p.size() ? get_deg(p[0], gen_degs) : -1; }
inline Deg get_deg(const Poly& p, const std::vector<Deg>& gen_degs) { return p.size() ? get_deg(p[0], gen_degs) : Deg({ -1, -1, -1 }); }
Mon gcd(const Mon& m1, const Mon& m2);
Mon lcm(const Mon& m1, const Mon& m2);

Poly get_diff(const Mon& mon, const Poly1d& diffs);
Poly get_diff(const Poly& poly, const Poly1d& diffs);

/* Linear Algebra Mod 2
**
** Vector spaces are sparse triangular matrices
*/

inline array range(int n) {
	array result;
	for (int i = 0; i < n; ++i)
		result.push_back(i);
	return result;
};

inline array range(int a, int b) {
	array result;
	for (int i = a; i < b; ++i)
		result.push_back(i);
	return result;
};

inline array add_vectors(const array& v1, const array& v2) {
	array result;
	std::set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(result));
	return result;
};
/* Reduce the space to the rref form */
array2d& simplify_space(array2d& spaceV);

array residue(const array2d& spaceV, const array& v);
array residue(const array2d& spaceV, array&& v);
/* Setup a linear map */
void set_linear_map(const array2d& fx, array2d& image, array2d& kernel, array2d& g);
/* Setup a linear map */
void set_linear_map(const array& x, const array2d& fx, array2d& image, array2d& kernel, array2d& g);
/* Compute the image of v under a linear map */
array get_image(const array2d& spaceV, const array2d& f, const array& v);
/* Compute the quotient of linear spaces V/W assuming that W is a subspace of V */
array2d quotient_space(const array2d& spaceV, const array2d& spaceW);

/* Groebner Basis
**
** The default monomial ordering is revlex
*/
/* Move from and remove top() of heap */
inline PolyWithT MoveFromTop(RelHeap& heap) {
	PolyWithT result = std::move(const_cast<PolyWithT&>(heap.top()));
	heap.pop();
	return result;
}
/* Reduce `poly` by groebner basis `rels` */
Poly reduce(Poly poly, const Poly1d& gb);
/* Comsume relations from heap that is at most in degree `deg` while adding new relations to heap that is at most in degree `deg_max`. */
void add_rels_from_heap(Poly1d& gb, RelHeap& heap, const array& gen_degs, int t, int deg_max);
/* Add new relations `rels` to groebner basis `gb` */
void add_rels(Poly1d& gb, const Poly1d& rels, const array& gen_degs, int deg_max);
/* return a_{ij} such that a_{i1}p_1+...+a_{in}p_n=0 */
Poly2d ann_seq(const Poly1d& gb, const Poly1d& polys, const array& gen_degs, int deg_max);

template <typename Fn>
Poly evaluate(const Poly& poly, Fn map, const Poly1d& gb)
{
	Poly result;
	for (const Mon& m : poly) {
		Poly fm = { {} };
		for (MonInd p = m.begin(); p != m.end(); ++p)
			fm = reduce(mul(fm, pow(map(p->gen), p->exp, gb)), gb);
		result = add(result, fm);
	}
	return result;
}

#endif /* MYMATH_H */
