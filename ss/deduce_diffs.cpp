#include "main.h"
#include "benchmark.h"
#include "linalg.h"
#include "myio.h"
#include <iostream>
#include <sstream>

/* If n = 2^k1 + ... + 2^kn,
 * return the array k1, ..., kn. */
array two_expansion(int n)
{
	array result;
	int k = 0;
	while (n > 0) {
		if (n & 1)
			result.push_back(k);
		n >>= 1;
		++k;
	}
	return result;
}

int DeduceDiffsForEt(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, std::vector<std::map<Deg, Staircase>>& basis_ss)
{
	int count_new_diffs = 0;
	/* Look for pairs of basis that can be only killed in one way */
	for (auto p = basis_ss.front().begin(); p != basis_ss.front().end(); ++p)
	{
		const Deg deg = p->first;
		const Staircase& sc = GetRecentStaircase(basis_ss, deg);

		/* Check if there is only one base with undermined differential */
		auto pair = CountUndetermined(sc, 0);
		int count = pair.first, index = pair.second;
		int count_tgt = 0, index_tgt = -1;
		int count_src = 0, index_src = -1;
		Deg deg_tgt{ -1, -1, -1 }, deg_src{ -1, -1, -1 };
		if (count == 1) { /* If there is only one base, check if there is only one way to kill it */
			int level = sc.levels[index];
			if (level > kLevelMax / 2) {
				/* Count how many chances it can be a source of a differential */
				for (int r = kLevelMax - level; r <= deg.v; r += 2) {
					Deg d_tgt = deg + Deg{ 1, 0, -r };
					if (basis_ss.front().find(d_tgt) != basis_ss.front().end()) {
						auto pair_tgt = CountUndetermined(GetRecentStaircase(basis_ss, d_tgt), 0);
						count_tgt += pair_tgt.first;
						if (pair_tgt.first > 0) {
							index_tgt = pair_tgt.second;
							deg_tgt = d_tgt;
						}
					}
				}
			}
			/* Count how many chances it can be a target of a differential */
			int r_max = std::min(level, 2 * (deg.t - deg.s + 1) - deg.v);
			for (int r = kLevelMin; r <= r_max; r += 2) {
				Deg d_src = deg - Deg{ 1, 0, -r };
				if (basis_ss.front().find(d_src) != basis_ss.front().end()) {
					auto pair_src = CountUndetermined(GetRecentStaircase(basis_ss, d_src), kLevelMax - r);
					count_src += pair_src.first;
					if (pair_src.first > 0) {
						index_src = pair_src.second;
						deg_src = d_src;
					}
				}
			}

			if (count_tgt == 0 && count_src == 1) { /* It is a target */
				if (IsSingle(sc, index)) {
					std::cout << "type1: (" << deg_src.t - deg_src.s << ", " << deg_src.v / 2 << ") --> (" << deg.t - deg.s << ", " << deg.v / 2 << ")\n";
					++count_new_diffs;
					SetDiff(gb, basis, basis_ss, deg_src, GetRecentStaircase(basis_ss, deg_src).basis_ind[index_src], sc.basis_ind[index], deg_src.v - deg.v, true);
				}
			}
			else if (count_tgt == 1 && count_src == 0) { /* It is a source */
				Staircase sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
				if (IsSingle(sc_tgt, index_tgt)) {
					std::cout << "type2: (" << deg.t - deg.s << ", " << deg.v / 2 << ") --> (" << deg_tgt.t - deg_tgt.s << ", " << deg_tgt.v / 2 << ")\n";
					++count_new_diffs;
					SetDiff(gb, basis, basis_ss, deg, sc.basis_ind[index], sc_tgt.basis_ind[index_tgt], deg.v - deg_tgt.v, true);
				}
			}
			else if (count_tgt == 0 && count_src > 1) { /* It is a permanent cycle */
				if (sc.levels[index] > kLevelMax / 4 * 3) {
					std::cout << "type3: (" << deg.t - deg.s << ", " << deg.v / 2 << ") --> 0\n";
					++count_new_diffs;
					SetDiff(gb, basis, basis_ss, deg, sc.basis_ind[index], array{}, kLevelMax / 4, true); //
				}
			}
			else if ((level < kLevelMax / 2 || count_tgt == 0) && count_src == 0) /* It is a nontrivial permanent cycle */
				throw SSException(0x9eacf10fU, "Nontrivial permanant cycle");
		}
	}
	return count_new_diffs;
}


/* Check if there exists at least one global differential structure */
bool FindGlobalStructureForEt(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, std::vector<std::map<Deg, Staircase>>& basis_ss)
{
	const size_t initial_history_length = basis_ss.size();
	CacheDeduction cache_history;
	cache_history.push_back();
	cache_history.InitNullDiffs(basis_ss);

	while (true) { /* Deep search for a global structure */
		/* Find the place to start in cached history */
		Deg deg = cache_history.degs_.back();
		if (deg.s == -1)
			return true;
		const NullDiff& nd = cache_history.GetRecentNullDiff(deg);
		if (nd.num_tgts > 10) {
			std::cout << "Too many choices!\n";
			basis_ss.resize(initial_history_length);
			return false;
		}
		int& i = cache_history.indices_.back();
		if (++i == (1 << nd.num_tgts)) { /* All targets fail. Roll back. */
			if (cache_history.degs_.size() > 1) {
				basis_ss.pop_back();
				cache_history.pop_back();
				continue;
			}
			else
				return false; /* No global structure */
		}

		const Staircase& sc = GetRecentStaircase(basis_ss, deg);
		const int index = nd.index;
		const int r = kLevelMax - sc.levels[index];
		const Deg deg_tgt = deg + Deg{ 1, 0, -r };
		const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
		
		/* Calculate the diff */
		array diff;
		for (int j : two_expansion(i))
			diff = lina::AddVectors(diff, sc_tgt.basis_ind[(size_t)nd.first_r + j]);

		try {
			array src = sc.basis_ind[index];
			std::cout << "Try: " << deg << " d_{" << r << '}' << src << '=' << diff << '\n';
			basis_ss.push_back({});
			cache_history.push_back();
			SetDiff(gb, basis, basis_ss, deg, src, diff, r, true);
			while (DeduceDiffsForEt(gb, basis, basis_ss));
			cache_history.InitNullDiffs(basis_ss);
		}
		catch (SSException&) { /* Try the next possibility */
			std::cout << "Failed!\n";
			basis_ss.pop_back();
			cache_history.pop_back();
		}
	}
}

int DeduceDiffsByTrying(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, std::vector<std::map<Deg, Staircase>>& basis_ss)
{
	return 0;
}

/* generate the table of the spectral sequence */
void DeduceDiffs(const Database& db, const std::string& table_prefix, int t_max)
{
	const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", t_max);
	const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", t_max);
	std::vector<std::map<Deg, Staircase>> basis_ss = { db.load_basis_ss(table_prefix + "_ss", t_max) };

	try {
		while (DeduceDiffsForEt(gb, basis, basis_ss));
	}
#ifdef _DEBUG
	catch (SSException& e)
	{
		std::cerr << "Error-" << std::hex << e.id() << ": " << e.what() << '\n';
	}
#else
	catch (NoException& e) {} /* This does not catch any exception */
#endif
	std::cout << "Find global differential structure:\n";
	if (FindGlobalStructureForEt(gb, basis, basis_ss)) {
		ApplyChanges(basis_ss);
		std::cout << "A global structure is found!\n";
	}
	else {
		std::cout << "no global structure\n";
	}
	
	
	/* insert into the database */
	db.begin_transaction();
	db.execute_cmd("DELETE FROM " + table_prefix + "_ss_tmp;");
	db.save_basis_ss(table_prefix + "_ss_tmp", basis_ss.front());
	db.end_transaction();
}
