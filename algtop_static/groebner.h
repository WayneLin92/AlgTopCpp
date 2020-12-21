/* 
** Groebner Basis
** The default monomial ordering is revlex
*/

#ifndef GROEBNER_H
#define GROEBNER_H

#include "algebras.h"
#include <map>
#include <queue>

namespace grbn {

/********** STRUCTS AND CLASSES **********/

struct PolyWithT {
	Poly poly;
	int t;
	bool operator<(const PolyWithT& rhs) const { return t < rhs.t; }
	bool operator>(const PolyWithT& rhs) const { return t > rhs.t; }
};

using RelHeap = std::priority_queue<PolyWithT, std::vector<PolyWithT>, std::greater<PolyWithT>>;

struct MonWithIndices {
	Mon gcd;
	int index1;
	int index2;
};

/* Reduce `poly` by groebner basis `gb` */
Poly Reduce(Poly poly, const Poly1d& gb);

class RelBuffer {
public:
	void push(int deg, const Mon& gcd, int index1, int index2) { gcds[deg].push_back({ gcd, index1, index2 }); }
	void push(int deg, Mon&& gcd, int index1, int index2) { gcds[deg].push_back({ gcd, index1, index2 }); }
	MonWithIndices pop() {
		MonWithIndices m = std::move(gcds.begin()->second.back());
		if (gcds.begin()->second.size() == 1)
			gcds.erase(gcds.begin());
		else
			gcds.begin()->second.pop_back();
		return m;
	}
	bool empty() const { return gcds.empty(); }
	int next_deg() const { return gcds.begin()->first; }
private:
	std::map<int, std::vector<MonWithIndices>> gcds;
};

/********** FUNCTIONS **********/

/* Move from and remove top() of heap */
inline PolyWithT MoveFromTop(RelHeap& heap) {
	PolyWithT result = std::move(const_cast<PolyWithT&>(heap.top()));
	heap.pop();
	return result;
};

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

/* Comsume relations from heap that is at most in degree `deg` while adding new relations to `heap` that is at most in degree `deg_max`. */
void add_rels_from_heap(Poly1d& gb, RelHeap& heap, const array& gen_degs, int t, int deg_max);
void add_rels_from_heap(Poly1d& gb, RelHeap& heap, const array& gen_degs, const array& gen_degs1, int t, int deg_max);

/* Add new relations `rels` to groebner basis `gb` */
void add_rels(Poly1d& gb, const Poly1d& rels, const array& gen_degs, int deg_max);
void add_rels(Poly1d& gb, const Poly1d& rels, const array& gen_degs, const array& gen_degs1, int deg_max);

void add_rels_freemodule(Poly1d& gb, RelHeap& heap, const array& gen_degs, const array& gen_degs1, int deg, int deg_max);

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
			fm = Reduce(mul(fm, pow(map(p->gen), p->exp, gb)), gb);
		result = add(result, fm);
	}
	return result;
}



/* Add new relations `rels` to groebner basis `gb` */
void AddRels(Poly1d& gb, const Poly1d& rels, const array& gen_degs, int deg_max);
void AddRels(Poly1d& gb, const Poly1d& rels, const array& gen_degs, const array& gen_degs1, int deg_max);
void AddRels(Poly1d& gb, const Poly1d& rels, RelBuffer& buffer, const array& gen_degs, int t, int deg_max);
void AddRels(Poly1d& gb, const Poly1d& rels, RelBuffer& buffer, const array& gen_degs, const array& gen_degs1, int t, int deg_max);

/* Assume gb is truncated in degree < t. Prepare the buffer in degrees [t, t_max] */
RelBuffer GenerateBuffer(const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t, int t_max);

} /* namespace grbn */

#endif /* GROEBNER_H */
