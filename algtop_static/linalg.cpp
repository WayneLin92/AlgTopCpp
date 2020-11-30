#include "linalg.h"
#include <algorithm>
#include <iterator>

array add_vectors(const array& v1, const array& v2) {
	array result;
	std::set_symmetric_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(result));
	return result;
};

array2d& simplify_space(array2d& spaceV)
{
	for (size_t i = spaceV.size() - 1; i != -1; i--)
		for (size_t j = 0; j < i; j++)
			if (std::binary_search(spaceV[j].begin(), spaceV[j].end(), spaceV[i][0]))
				spaceV[j] = add_vectors(spaceV[j], spaceV[i]);
	return spaceV;
}

array residue(const array2d& spaceV, const array& v)
{
	array result(v);
	for (size_t i = 0; i < spaceV.size(); ++i)
		if (std::binary_search(result.begin(), result.end(), spaceV[i][0]))
			result = add_vectors(result, spaceV[i]);
	return result;
}

array residue(const array2d& spaceV, array&& v)//
{
	array result(v);
	for (size_t i = 0; i < spaceV.size(); ++i)
		if (std::binary_search(result.begin(), result.end(), spaceV[i][0]))
			result = add_vectors(result, spaceV[i]);
	return result;
}

inline void add_map(array2d& spaceV, const array& v)
{
	array v1 = residue(spaceV, v);
	if (!v1.empty())
		spaceV.push_back(std::move(v1));
}

void set_linear_map(const array2d& fx, array2d& image, array2d& kernel, array2d& g)
{
	/* f(g[i]) = image[i] */
	for (size_t i = 0; i < fx.size(); ++i) {
		array src = { int(i) };
		array tgt(fx[i]);
		for (size_t j = 0; j < image.size(); j++) {
			if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
				tgt = add_vectors(tgt, image[j]);
				src = add_vectors(src, g[j]);
			}
		}
		if (tgt.empty())
			add_map(kernel, src);
		else {
			image.push_back(std::move(tgt));
			g.push_back(std::move(src));
		}
	}
}

void set_linear_map(const array& x, const array2d& fx, array2d& image, array2d& kernel, array2d& g)
{
	/* f(g[i]) = image[i] */
	for (size_t i = 0; i < fx.size(); ++i) {
		array src = { x[i] };
		array tgt(fx[i]);
		for (size_t j = 0; j < image.size(); j++) {
			if (std::binary_search(tgt.begin(), tgt.end(), image[j][0])) {
				tgt = add_vectors(tgt, image[j]);
				src = add_vectors(src, g[j]);
			}
		}
		if (tgt.empty())
			add_map(kernel, src);
		else {
			image.push_back(std::move(tgt));
			g.push_back(std::move(src));
		}
	}
}

/* Return f(v) for v\\in V. fi=f(vi) */
array get_image(const array2d& spaceV, const array2d& f, const array& v)
{
	array result;
	array v1(v);
	for (size_t i = 0; i < spaceV.size() && !v1.empty(); ++i)
		if (std::binary_search(v1.begin(), v1.end(), spaceV[i][0])) {
			v1 = add_vectors(v1, spaceV[i]);
			result = add_vectors(result, f[i]);
		}
#if _DEBUG
	if (!v1.empty()) {
		std::cerr << "v is not in space V!\n";
		throw "6a4fe8a1";
	}
#endif
	return result;
}

array2d quotient_space(const array2d& spaceV, const array2d& spaceW)
{
	array2d quotient;
	size_t dimQuo = spaceV.size() - spaceW.size();
#ifdef _DEBUG
	for (size_t i = 0; i < spaceV.size(); i++)
#else
	for (size_t i = 0; i < spaceV.size() && quotient.size() < dimQuo; i++)
#endif
	{
		auto v1 = residue(quotient, residue(spaceW, spaceV[i]));
		if (!v1.empty())
			quotient.push_back(std::move(v1));
	}
#ifdef _DEBUG
	if (quotient.size() != dimQuo) {
		std::cerr << "W is not a subspace of V!\n";
		throw "cec7f701";
	}
#endif
	return quotient;
}
