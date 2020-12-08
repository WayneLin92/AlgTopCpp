/* 
** Groebner Basis
** The default monomial ordering is revlex
*/

#ifndef GROEBNER_H
#define GROEBNER_H

#include "algebras.h"

/********** STRUCTS AND CLASSES **********/

struct PolyWithT {
	Poly poly;
	int t;
	bool operator<(const PolyWithT& rhs) const { return t < rhs.t; }
	bool operator>(const PolyWithT& rhs) const { return t > rhs.t; }
};

using RelHeap = std::priority_queue<PolyWithT, std::vector<PolyWithT>, std::greater<PolyWithT>>;

/********** FUNCTIONS **********/

/* Move from and remove top() of heap */
inline PolyWithT MoveFromTop(RelHeap& heap) {
	PolyWithT result = std::move(const_cast<PolyWithT&>(heap.top()));
	heap.pop();
	return result;
};

template <typename Container1d>
inline void RemoveEmptyElements(Container1d& gb)
{
	gb.erase(std::remove_if(gb.begin(), gb.end(), [](const Container1d::value_type& g) {return g.empty(); }), gb.end());
}

inline array range(int n) {
	array result;
	for (int i = 0; i < n; ++i)
		result.push_back(i);
	return result;
};

Poly pow(const Poly& poly, int n, const Poly1d& gb);

/* Reduce `poly` by groebner basis `rels` */
Poly reduce(Poly poly, const Poly1d& gb);

/* Comsume relations from heap that is at most in degree `deg` while adding new relations to `heap` that is at most in degree `deg_max`. */
void add_rels_from_heap(Poly1d& gb, RelHeap& heap, const array& gen_degs, int t, int deg_max);
void add_rels_from_heap(Poly1d& gb, RelHeap& heap, const array& gen_degs, const array& gen_degs1, int t, int deg_max);

/* Add new relations `rels` to groebner basis `gb` */
void add_rels(Poly1d& gb, const Poly1d& rels, const array& gen_degs, int deg_max);
void add_rels(Poly1d& gb, const Poly1d& rels, const array& gen_degs, const array& gen_degs1, int deg_max);

/* return a_{ij} such that a_{i1}p_1+...+a_{in}p_n=0 */
Poly2d ann_seq(const Poly1d& gb, const Poly1d& polys, const array& gen_degs, int deg_max);

/* Assume gb is truncated in degree < t. Prepare the heap in degrees [t, t_max] */
RelHeap GenerateHeap(const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t, int t_max);

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

#endif /* GROEBNER_H */
