/*
** Groebner Basis
** The default monomial ordering is revlex
*/

#ifndef GROEBNER_H
#define GROEBNER_H

#include "algebras.h"
#include <map>
#include <memory>

//#define GROEBNER_MULTITHREAD

namespace grbn {

/**********************************************************
* Small functions
**********************************************************/
/* Reduce `poly` by groebner basis `gb` */
Poly Reduce(Poly poly, const Poly1d& gb);

template <typename Container1d>
inline void RemoveEmptyElements(Container1d& cont)
{
	cont.erase(std::remove_if(cont.begin(), cont.end(), [](const Container1d::value_type& g) {return g.empty(); }), cont.end());
}

inline array range(int n) {
	array result;
	for (int i = 0; i < n; ++i)
		result.push_back(i);
	return result;
};

Poly pow(const Poly& poly, int n, const Poly1d& gb);

template <typename Fn>
Poly evaluate(const Poly& poly, Fn map, const Poly1d& gb)
{
	Poly result;
	for (const Mon& m : poly) {
		Poly fm = { {} };
		for (MonInd p = m.begin(); p != m.end(); ++p)
			fm = Reduce(mul(fm, pow(map(p->gen), p->exp, gb)), gb);
		result = add(result, fm);
	}
	return result;
}

/**********************************************************
* Groebner basis
**********************************************************/

using GbBuffer = std::map<int, Poly1d>;

void AddRelsB(Poly1d& gb, GbBuffer& buffer, const array& gen_degs, int t, int deg_max);
void AddRelsB(Poly1d& gb, GbBuffer& buffer, const array& gen_degs, const array& gen_degs1, int t, int deg_max);
void AddRels(Poly1d& gb, Poly1d rels, const array& gen_degs, int deg_max);
void AddRels(Poly1d& gb, Poly1d rels, const array& gen_degs, const array& gen_degs1, int deg_max);
void AddRelsM(Poly1d& gb, GbBuffer& buffer, const array& gen_degs, const array& gen_degs1, int deg, int deg_max);

/* Assume gb is truncated in degree < t. Prepare the heap in degrees [t, t_max] */
GbBuffer GenerateBuffer(const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t, int t_max);


/**********************************************************
* Groebner basis version 2
* This implementation uses polymorphism to reduce the use of memory
* with little cost.
**********************************************************/

struct BaseEleBuffer
{
	virtual Poly GetPoly(const Poly1d& gb) = 0;
};

struct GcdEleBuffer : BaseEleBuffer
{
	Mon gcd_; int i1_, i2_;
	GcdEleBuffer(Mon gcd, int i1, int i2) : gcd_(std::move(gcd)), i1_(i1), i2_(i2) {};
	Poly GetPoly(const Poly1d& gb) override {
		Poly result = gb[i1_] * div(gb[i2_][0], gcd_) + gb[i2_] * div(gb[i1_][0], gcd_);
		Mon{}.swap(gcd_); /* Deallocate */
		return result;
	}
};

struct PolyEleBuffer : BaseEleBuffer
{
	Poly p_;
	PolyEleBuffer(Poly p) : p_(std::move(p)) {};
	Poly GetPoly(const Poly1d& gb) override {
		return std::move(p_);
	}
};

using GbBufferV2 = std::map<int, std::vector<std::unique_ptr<BaseEleBuffer>>>;

void AddRelsBV2(Poly1d& gb, GbBufferV2& buffer, const array& gen_degs, int t, int deg_max);
void AddRelsBV2(Poly1d& gb, GbBufferV2& buffer, const array& gen_degs, const array& gen_degs1, int t, int deg_max);
void AddRelsV2(Poly1d& gb, Poly1d rels, const array& gen_degs, int deg_max);
void AddRelsV2(Poly1d& gb, Poly1d rels, const array& gen_degs, const array& gen_degs1, int deg_max);
void AddRelsMV2(Poly1d& gb, GbBufferV2& buffer, const array& gen_degs, const array& gen_degs1, int deg, int deg_max);

/* Assume gb is truncated in degree < t. Prepare the buffer in degrees [t, t_max] */
GbBufferV2 GenerateBufferV2(const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t, int t_max);

/**********************************************************
* Algorithms that use Groebner basis
**********************************************************/

/* return a_{ij} such that a_{i1}p_1+...+a_{in}p_n=0 */
Poly2d ann_seq(const Poly1d& gb, const Poly1d& polys, const array& gen_degs, int deg_max);

} /* namespace grbn */

#endif /* GROEBNER_H */
