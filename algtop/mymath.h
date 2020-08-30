#ifndef MYMATH_H
#define MYMATH_H
// TODO: greobner basis array4d version
// TODO: monomial std::vector<MonPow> version

#include "myparser.h"
#include <algorithm>
#include <numeric>

/********** STRUCTS AND CLASSES **********/
/* sparse monomial i1, e1, i2, e2, ... */
using array = std::vector<int>;
using array2d = std::vector<array>;
using array3d = std::vector<array2d>;
using array4d = std::vector<array3d>;
using arrayInd = array::const_iterator;

struct MonPow { int gen, exp; };
using Mon = std::vector<MonPow>;
using Poly = std::vector<Mon>;
using Poly1d = std::vector<Poly>;
using Poly2d = std::vector<Poly1d>;

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

struct rel_heap_t {
	array2d poly;
	int deg;
};

/********** FUNCTIONS **********/
// Move to a header for python in the future

inline void load_array(array& mon, std::istream& sin)         { load_vector(mon, sin, "(", ",", ")"); } // TODO: check if empty string is properly handled
inline void load_array2d(array2d& poly, std::istream& sin)      { load_vector(poly, sin, "{", ",", "}", load_array); }
inline void load_array3d(array3d& polys, std::istream& sin)   { load_vector(polys, sin, "[", ",", "]", load_array2d); }
inline void dump_array(const array& mon, std::ostream& sout)              { dump_vector(mon, sout, "(", ", ", ")"); }
inline void dump_array2d(const array2d& poly, std::ostream& sout)           { dump_vector(poly, sout, "{", ", ", "}", dump_array); }
inline void dump_array3d(const array3d& polys, std::ostream& sout)        { dump_vector(polys, sout, "[", ", ", "]", dump_array2d); }
inline void dump_array4d(const array4d& a, std::ostream& sout) { dump_vector(a, sout, "[", ", ", "]", dump_array3d); }

inline std::istream& operator>>(std::istream& sin, array& mon)        { load_array(mon, sin); return sin; }
inline std::istream& operator>>(std::istream& sin, array2d& poly)      { load_array2d(poly, sin); return sin; }
inline std::istream& operator>>(std::istream& sin, array3d& polys)    { load_array3d(polys, sin); return sin; }
inline std::ostream& operator<<(std::ostream& sout, const array& mon)        { dump_array(mon, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const array2d& poly)      { dump_array2d(poly, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const array3d& polys)    { dump_array3d(polys, sout); return sout; }
inline std::ostream& operator<<(std::ostream& sout, const array4d& a) { dump_array4d(a, sout); return sout; }

inline std::ostream& operator<<(std::ostream& sout, const Deg& d) { sout << '(' << d.s << ',' << d.t << ',' << d.v << ')'; return sout; }

inline void hash_combine(std::size_t& seed, int v) { seed ^= std::hash<int>{}(v)+0x9e3779b9 + (seed << 6) + (seed >> 2); }
inline size_t hash_range(const array& mon)
{
	size_t seed = 0;
	for (int i : mon)
		hash_combine(seed, i);
	return seed;
}

/*--------- algmod2.cpp ---------*/
/* Functions for Monomials and Polynomials
**
** A monomial m is an array of i1, e1, i2, e2, ...
** Monomials are ordered lexicographically by exponents after filling in zero exponents.
** A polynomials is an increasing squence of monomials.
*/

array mul(const array& mon1, const array& mon2);
array div(const array& mon1, const array& mon2);
array2d add(const array2d& poly1, const array2d& poly2);
array2d mul(const array2d& poly, const array& mon); // TODO: add begin-end version of these functions
array2d mul(const array2d& poly1, const array2d& poly2);
array pow(const array& m, int e);
array2d pow(const array2d& poly, int n);
array2d pow(const array2d& poly, int n, const array3d& gb);
bool divides(const array& m1, const array& m2);
bool divides(arrayInd pFirst1, arrayInd pLast1, arrayInd pFirst2, arrayInd pLast2);
/* return the largest integer e where m1 = m2^e * r */
int log(const array& m1, const array& m2);

inline int get_deg(const array& mon) { int result = 0; for (size_t i = 0; i < mon.size(); i += 2) result += mon[i + 1]; return result; };
inline int get_deg(const array2d& poly) { return poly.size() ? get_deg(poly[0]) : -1; };
inline int get_deg(const array& mon, const array& gen_degs) { int result = 0; for (size_t i = 0; i < mon.size(); i += 2) result += gen_degs[mon[i]] * mon[i + 1]; return result; };
inline Deg get_deg(const array& mon, const std::vector<Deg>& gen_degs) { Deg result({ 0, 0, 0 }); for (size_t i = 0; i < mon.size(); i += 2) result += gen_degs[mon[i]] * mon[i + 1]; return result; };
inline int get_deg(const array2d& p, const array& gen_degs) { return p.size() ? get_deg(p[0], gen_degs) : -1; }
inline Deg get_deg(const array2d& p, const std::vector<Deg>& gen_degs) { return p.size() ? get_deg(p[0], gen_degs) : Deg({ -1, -1, -1 }); }
/* cmp_mons(m1, m2) returns true if e1 < e2 in lexicographical order where e1, e2 are arrays of exponents. */
bool cmp_mons(const array& m1, const array& m2);
array gcd(const array& m1, const array& m2);
array lcm(const array& m1, const array& m2);

array2d get_diff(const array& mon, const array3d& diffs);
array2d get_diff(const array2d& poly, const array3d& diffs);

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

/* Reduce `poly` by groebner basis `rels` */
array2d reduce(array2d poly, const array3d& gb);
inline bool cmp_heap_rels(const rel_heap_t& s1, const rel_heap_t& s2) { return s1.deg > s2.deg; }
/* Comsume relations from heap that is at most in degree `deg` while adding new relations to heap that is at most in degree `deg_max`. */
void add_rels(array3d& gb, std::vector<rel_heap_t>& heap, const array& gen_degs, int t, int deg_max);
/* Add new relations `rels` to groebner basis `gb` */
void add_rels(array3d& gb, const array3d& rels, const array& gen_degs, int deg_max);
/* return a_{ij} such that a_{i1}p_1+...+a_{in}p_n=0 */
array4d ann_seq(const array3d& gb, const array3d& polys, const array& gen_degs, int deg_max);

template <typename Fn>
array2d evaluate(const array2d& poly, Fn map, const array3d& gb)
{
	array2d result;
	for (const array& m : poly) {
		array2d fm = { {} };
		for (size_t i = 0; i < m.size(); i += 2)
			fm = reduce(mul(fm, pow(map(m[i]), m[i + 1], gb)), gb);
		result = add(result, fm);
	}
	return result;
}

#endif /* MYMATH_H */
