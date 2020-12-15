#ifndef MYCUDA_PUBLIC_H
#define MYCUDA_PUBLIC_H

#include <vector>

namespace cuda {
    using array = std::vector<int>;
    using array2d = std::vector<array>;

	array2d EchelonCuda(const array2d& matrix_csr);
}


#endif /* MYCUDA_PUBLIC_H */
