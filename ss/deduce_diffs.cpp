#include "main.h"
#include "groebner.h"
#include "benchmark.h"
#include "linalg.h"
#include "myexception.h"
#include <iostream>

constexpr auto kLevelMin = 4;

/* A custom Exception class */
class SSException : public MyException
{
public:
	SSException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
	SSException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

/* Add d_r(base)=diff. Return empty if nothing has changed. */
BasisSSV2 add_diff(const BasisSSV2& basis_d, const array& base, const array& diff, int r)
{
	array2d cycles_r; /* Ker(d_r) */
	array2d cycles_rm2_with_diff; /* Ker(d_{r-2}) */
	array2d cycles_rm2_without_diff; /* Ker(d_{r-2}) */
	array2d diffs_of_rm2;

	/* Initialize cycles_r and cycles_rm2 */
	int level_cycles_rm2 = r > kLevelMin ? kLevelMax - r + 2 : kLevelMax;
	for (int i = 0; i < (int)basis_d.levels.size(); ++i) {
		if (basis_d.levels[i] <= kLevelMax - r)
			cycles_r.push_back(basis_d.basis_ind[i]);
		else if (basis_d.levels[i] <= level_cycles_rm2) {
			if (basis_d.diffs_ind[i] != array{ -1 }) {
				cycles_rm2_with_diff.push_back(basis_d.basis_ind[i]);
				diffs_of_rm2.push_back(basis_d.diffs_ind[i]);
			}
			else
				array2d cycles_rm2_without_diff; /* Ker(d_{r-2}) */

		}
		else
			break;
	}

	/* Check if base is in the vector space cycles_rm2_without_diff mod cycles_r + cycles_rm2_with_diff */
	array base1 = lina::Residue(cycles_rm2_with_diff, lina::Residue(cycles_r, base));
	if (!lina::Residue(cycles_rm2_without_diff, base1).empty())
		throw SSException(0xa34147cbU, "The source is not in the E_r page.");
	if (base1.empty())
		;

}

/* Add new diff d_r(basis_ss[deg][i])=diff */
void set_diff(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, std::vector<std::map<Deg, BasisSSV2>>& basis_ss, const Deg& deg, const array& base, const array& diff, int r)
{
	for (auto& [d, basis_d] : basis_ss.front()) {
		;
	}
}

/* generate the table of the spectral sequence */
void deduce_diffs(const Database& db, const std::string& table_prefix, int r, int t_max)
{
	const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", t_max);
	const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", t_max);
	std::vector<std::map<Deg, BasisSSV2>> basis_ss = { db.load_basis_ss(table_prefix + "_ss", t_max) };


}