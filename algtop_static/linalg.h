/* 
** Linear Algebra Mod 2
**
** Vector spaces are sparse triangular matrices
*/

#ifndef LINALG_H
#define LINALG_H

#include <vector>

using array = std::vector<int>;
using array2d = std::vector<array>;

array add_vectors(const array& v1, const array& v2);

/* Reduce the space to the rref form */
array2d& simplify_space(array2d& spaceV);

array residue(const array2d& spaceV, const array& v);
array residue(const array2d& spaceV, array&& v);

/* Setup a linear map */
void set_linear_map(const array2d& fx, array2d& image, array2d& kernel, array2d& g);

/* Setup a linear map */
void set_linear_map_v2(const array& x, const array2d& fx, array2d& image, array2d& kernel, array2d& g);

/* Compute the image of v under a linear map */
array get_image(const array2d& spaceV, const array2d& f, const array& v);

/* Compute the quotient of linear spaces V/W assuming that W is a subspace of V */
array2d quotient_space(const array2d& spaceV, const array2d& spaceW);

#endif /* LINALG_H */