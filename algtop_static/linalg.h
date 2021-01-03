/*****************************************************************
** Linear Algebra Mod 2
** 
** All row vectors are compressed
** Spaces are upper triangular matrices with some column ordering
*****************************************************************/

#ifndef LINALG_H
#define LINALG_H

#include <vector>
#include <algorithm>
#include <iterator>

using array = std::vector<int>;
using array2d = std::vector<array>;

namespace lina {

/* Add two compressed vectors */
inline array AddVectors(const array& v1, const array& v2) {
	array result;
	std::set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(result));
	return result;
};

/* Return the space spanned by `vectors` */
array2d GetSpace(const array2d& vectors);

/* Reduce the space to rref form */
array2d& SimplifySpace(array2d& spaceV);

/* Return a newer v such that (spaceV\\ v) is triangular */
array Residue(const array2d& spaceV, const array& v);
array Residue(const array2d& spaceV, array&& v);

/* Setup a linear map */
void GetInvMap(const array2d& fx, array2d& image, array2d& g);
void SetLinearMap(const array2d& fx, array2d& image, array2d& kernel, array2d& g);
void SetLinearMapV2(const array& x, const array2d& fx, array2d& image, array2d& kernel, array2d& g);
void SetLinearMapV3(const array2d& x, const array2d& fx, array2d& domain, array2d& f, array2d& image, array2d& g, array2d& kernel);

/* Return f(v) for v\\in V. fi=f(vi) */
array GetImage(const array2d& spaceV, const array2d& f, const array& v);

/* Compute the quotient of linear spaces V/W assuming that W is a subspace of V */
array2d QuotientSpace(const array2d& spaceV, const array2d& spaceW);

} /* namespace lina */

#endif /* LINALG_H */