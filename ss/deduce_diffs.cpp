#include "main.h"
#include "benchmark.h"
#include "linalg.h"
#include "myio.h"
#include <iostream>
#include <sstream>

#define DEDUCE_LOGGING

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

/* Deduce zero differentials for degree reason */
int DeduceZeroDiffs(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, bool bForEt = false)
{
	int count_new_diffs = 0;
	int prev_t = -1;
	for (auto p = basis_ss.front().begin(); p != basis_ss.front().end(); ++p) {
		const Staircase& sc = GetRecentStaircase(basis_ss, p->first);
		for (size_t i = 0; i < sc.levels.size(); ++i) {
			if (sc.diffs_ind[i] == array{ -1 } && sc.levels[i] > kLevelMax / 2) {
				int r = kLevelMax - sc.levels[i];
				auto [count_tgt, deg_tgt, index_tgt] = CountTgt(basis_ss, p->first, r);
				if (count_tgt > 0) {
					if (p->first.v - deg_tgt.v > r) {
						SetDiff(gb, basis, basis_ss, p->first, sc.basis_ind[i], {}, p->first.v - deg_tgt.v - 2);
						++count_new_diffs;
					}
				}
				else {
					SetDiff(gb, basis, basis_ss, p->first, sc.basis_ind[i], {}, kLevelMax / 4);
					if (bForEt)
						AddImage(basis_ss, p->first, sc.basis_ind[i], { -1 }, kLevelMax / 2 - 1);
				}
			}
		}
	}
	return count_new_diffs;
}

/* Use the fact everything in positive degree should be killed by differentials */
int DeduceDiffsForEt(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss)
{
	int count_new_diffs = 0;
	/* Look for pairs of basis that can be only killed in one way */
	for (auto p = basis_ss.front().begin(); p != basis_ss.front().end(); ++p) {
		const Deg deg = p->first;
		const Staircase& sc = GetRecentStaircase(basis_ss, deg);

		/* Check if there is only one base with undermined (inv)differential in its level */
		size_t index = sc.levels.size();
		while (index-- > 0) {
			if (sc.diffs_ind[index] == array{ -1 })
				break;
		}
		if (index < sc.levels.size() && (index == 0 || sc.levels[index] != sc.levels[index - 1])) {
			int level = sc.levels[index];
			auto [count_tgt, deg_tgt, index_tgt] = level > kLevelMax / 2 ? CountTgt(basis_ss, deg, kLevelMax - level) : std::make_tuple(0, Deg(), -1);
			auto [count_src, deg_src, index_src] = CountSrc(basis_ss, deg, level);
			
			if (count_tgt == 0 && count_src == 1) { /* It is a target */
				int r = deg_src.v - deg.v;
				if (index == 0 || sc.levels[index - 1] < r) {
#ifdef DEDUCE_LOGGING
					std::cout << deg << " <-- " << deg_src << '\n';
#endif
					++count_new_diffs;
					SetDiff(gb, basis, basis_ss, deg_src, GetRecentStaircase(basis_ss, deg_src).basis_ind[index_src], sc.basis_ind[index], deg_src.v - deg.v);
				}
			}
			else if (count_tgt == 1 && count_src == 0) { /* It is a source */
				Staircase sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
#ifdef DEDUCE_LOGGING
				std::cout << deg << " --> " << deg_tgt << '\n';
#endif
				++count_new_diffs;
				SetDiff(gb, basis, basis_ss, deg, sc.basis_ind[index], sc_tgt.basis_ind[index_tgt], deg.v - deg_tgt.v);
			}
			else if (count_tgt == 0 && count_src > 1) { /* It is a permanent cycle */
				if (sc.levels[index] > kLevelMax / 2) {
#ifdef DEDUCE_LOGGING
					std::cout << deg << " --> 0\n";
#endif
					++count_new_diffs;
					SetDiff(gb, basis, basis_ss, deg, sc.basis_ind[index], array{}, kLevelMax / 2 + 2);
					AddImage(basis_ss, deg, sc.basis_ind[index], { -1 }, kLevelMax / 2 - 1);
				}
			}
			else if ((level < kLevelMax / 2 || count_tgt == 0) && count_src == 0) /* It is a nontrivial permanent cycle */
				throw SSException(0x9eacf10fU, "Nontrivial permanant cycle");
		}
	}
	return count_new_diffs;
}

/* Check if there exists global differential structures
 * Return 0 if there is no global structure
 * Return 1 if there exists a global structure 
 * Return 2 if the search space is too big
 * Return 3 if the search depth is reached */
int FindGlobalStructure(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, int depth, bool bForEt = false)
{
	const size_t initial_history_length = basis_ss.size();
	CacheDeduction cache_history;
	cache_history.push_back();
	cache_history.InitNullDiffs(basis_ss);

	while (true) { /* Deep search for a global structure */
		/* Find the place to start in cached history */
		Deg deg = cache_history.degs_.back();
		if (deg.s == -1) {
			if (bForEt && !NoNullDiff(basis_ss)) { /* Invalid global structure. Roll back. */
				if (cache_history.degs_.size() > 1) {
					basis_ss.pop_back();
					cache_history.pop_back();
					continue;
				}
				else
					return 0; /* No global structure */
			}
			//basis_ss.resize(initial_history_length);
			return 1; /* Find one global structure */
		}
		if (depth != -1 && cache_history.degs_.size() > depth) {
			basis_ss.resize(initial_history_length);
			return 3; /* search depth is reached */
		}
		const NullDiff& nd = cache_history.GetRecentNullDiff(deg);
		if (nd.num_tgts > 10) {
			basis_ss.resize(initial_history_length);
			return 2; /* Too many choices */
		}
		int& i = cache_history.indices_.back();
		if (++i == (1 << nd.num_tgts)) { /* All targets fail. Roll back. */
			if (cache_history.degs_.size() > 1) {
				basis_ss.pop_back();
				cache_history.pop_back();
				continue;
			}
			else
				return 0; /* No global structure */
		}

		const Staircase& sc = GetRecentStaircase(basis_ss, deg);
		const int r = kLevelMax - sc.levels[nd.index];
		
		/* Calculate the diff */
		array src = sc.basis_ind[nd.index];
		array tgt;
		if (nd.num_tgts > 0) { /* tgt is zero when there is only one choice */
			const Deg deg_tgt = deg + Deg{ 1, 0, -r };
			const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
			for (int j : two_expansion(i))
				tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)nd.first_r + j]);
		}

		basis_ss.push_back({}); /* Warning: this invalidates references above */
		cache_history.push_back();
		try {
			SetDiff(gb, basis, basis_ss, deg, src, tgt, r);
			DeduceZeroDiffs(gb, basis, basis_ss, bForEt);
			if (bForEt)
				DeduceDiffsForEt(gb, basis, basis_ss);
			cache_history.InitNullDiffs(basis_ss);
		}
		catch (SSException&) { /* Try the next possibility */
			basis_ss.pop_back();
			cache_history.pop_back();
		}
	}
}

int DeduceDiffsByTrying(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, int depth, bool bForEt = false)
{
	NullDiffs nullDiffs;
	nullDiffs.InitNullDiffs(basis_ss, true);

	/* List the degrees */
	std::vector<Deg> degs;
	array nums_tgts;
	for (auto p = nullDiffs.null_diffs_.begin(); p != nullDiffs.null_diffs_.end(); ++p) {
		int num_tgts = nullDiffs.null_diffs_.at(p->first).num_tgts;
		if (num_tgts <= 10 && num_tgts >= 0) {
			degs.push_back(p->first);
			nums_tgts.push_back(num_tgts);
		}
	}
	std::sort(degs.begin(), degs.end(), [&degs, &nums_tgts](const Deg& d1, const Deg& d2) {return nums_tgts[&d1 - &degs[0]] < nums_tgts[&d2 - &degs[0]]; });

	for (auto& deg : degs) {
		const Staircase& sc = GetRecentStaircase(basis_ss, deg);
		const NullDiff& nd = nullDiffs.null_diffs_.at(deg);
		const int r = kLevelMax - sc.levels[nd.index];
		const Deg deg_tgt = deg + Deg{ 1, 0, -r };
		array src = sc.basis_ind[nd.index];
		array tgt_pass;

		int count_pass = 0;
		for (int i = 0; i < (1 << nd.num_tgts); ++i) {
			/* Calculate the diff */
			array tgt;
			if (nd.num_tgts > 0) {
				const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
				for (int j : two_expansion(i))
					tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)nd.first_r + j]);
			}
#ifdef DEDUCE_LOGGING
			std::cout << "Try " << i << '/' << (1 << nd.num_tgts) << ": " << deg << " d_{" << r << '}' << src << '=' << tgt;
#endif
			basis_ss.push_back({}); /* Warning: invalidate references above */
			bool bException = false;
			try {
				SetDiff(gb, basis, basis_ss, deg, src, tgt, r);
				if (bForEt)
					while (DeduceDiffsForEt(gb, basis, basis_ss));
				if (!FindGlobalStructure(gb, basis, basis_ss, depth, bForEt))
					bException = true;
			}
			catch (SSException&) {
				bException = true;
			}
			basis_ss.pop_back();
			if (!bException) {
				tgt_pass = tgt;
				++count_pass;
#ifdef DEDUCE_LOGGING
				std::cout << " Success\n";
#endif
				if (count_pass > 1)
					break;
			}
#ifdef DEDUCE_LOGGING
			else
				std::cout << " Fail\n";
#endif
		}
		if (count_pass == 0)
			throw SSException(0x5b3d7e35U, "no compatible differentials");
		else if (count_pass == 1) {
			basis_ss.push_back({});
			SetDiff(gb, basis, basis_ss, deg, src, tgt_pass, r);
			DeduceZeroDiffs(gb, basis, basis_ss, bForEt);
			if (bForEt)
				DeduceDiffsForEt(gb, basis, basis_ss);
			nullDiffs.InitNullDiffs(basis_ss, false);
			ApplyRecentChanges(basis_ss);
		}
	}
	return 0;
}

/* generate the table of the spectral sequence */
void DeduceDiffs(const Database& db, const std::string& table_prefix, int t_max, bool bForEt /*= false*/)
{
	const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", t_max);
	const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", t_max);
	Staircases1d basis_ss = { db.load_basis_ss(table_prefix + "_ss", t_max) };

	try {
		while (true) {
			int count1 = DeduceZeroDiffs(gb, basis, basis_ss, true);
			std::cout << "DeduceZeroDiffs: " << count1 << '\n';

			int count2 = DeduceDiffsForEt(gb, basis, basis_ss);
			std::cout << "DeduceDiffsForEt: " << count2 << '\n';

			//int count3 = DeduceDiffsByTrying(gb, basis, basis_ss, 0, true);
			//std::cout << "DeduceDiffsByTrying: " << count3 << '\n';

			if (count1 + count2 /*+ count3*/ == 0)
				break;
		}
	}
	catch (SSException& e)
	{
		std::cerr << "Error-" << std::hex << e.id() << ": " << e.what() << '\n';
	}
	FindGlobalStructure(gb, basis, basis_ss, -1, true);
	ApplyAllChanges(basis_ss);
	
	/* insert into the database */
	db.begin_transaction();
	db.execute_cmd("DELETE FROM " + table_prefix + "_ss_tmp;");
	db.save_basis_ss(table_prefix + "_ss_tmp", basis_ss.front());
	db.end_transaction();
}
