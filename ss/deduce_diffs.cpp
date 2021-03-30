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

/* Deduce zero differentials for degree reason 
 * t_test is for DeduceDiffs() */
int DeduceZeroDiffs(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, DiffType d_type, int t_max, bool bForEt, bool bPrint, int t_test)
{
	auto [s_diff, t_diff] = GetST(d_type);
	if (t_max == -1)
		t_max = basis.rbegin()->first.t;

	int count_new_diffs = 0;
	int prev_t = -1;
	for (auto p = basis_ss.front().begin(); p != basis_ss.front().end(); ++p) {
		if (bPrint && p->first.t != prev_t) {
			std::cout << p->first.t << " " << count_new_diffs << "     \r";
			prev_t = p->first.t;
		}
		if (p->first.t > t_max + t_diff) //TODO: Consider other places for t_diff
			break;
		const Staircase& sc = GetRecentStaircase(basis_ss, p->first);
		for (size_t i = 0; i < sc.levels.size(); ++i) {
			if (sc.diffs_ind[i] == array{ -1 } && sc.levels[i] > kLevelPC) {
				int r = kLevelMax - sc.levels[i];
				auto [count_tgt, deg_tgt, index_tgt] = CountTgt(basis_ss, p->first, d_type, r); /* Find the first possible d_{r1} target for r1>=r */
				if (count_tgt > 0) {
					if (p->first.v - deg_tgt.v > r) {
						SetDiff(gb, basis, basis_ss, p->first, d_type, sc.basis_ind[i], {}, p->first.v - deg_tgt.v - 2, t_test);
						++count_new_diffs;
					}
				}
				else {
					SetDiff(gb, basis, basis_ss, p->first, d_type, sc.basis_ind[i], {}, kLevelMax / 4, t_test);
					if (bForEt)
						AddImage(basis_ss, p->first, d_type, sc.basis_ind[i], { -1 }, kLevelMax / 2 - 1);
				}
			}
		}
	}
	return count_new_diffs;
}

/* Use the fact everything in positive degree should be killed by differentials */
int DeduceDiffsForEt(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, DiffType d_type, int t_max, bool bPrint = false, int t_test = -1)
{
	auto [s_diff, t_diff] = GetST(d_type);
	if (t_max == -1)
		t_max = basis.rbegin()->first.t;
	int count_new_diffs = 0;
	int prev_t = -1;
	/* Look for pairs of basis that can be only killed in one way */
	for (auto p = basis_ss.front().begin(); p != basis_ss.front().end(); ++p) {
		if (bPrint && p->first.t != prev_t) {
			std::cout << p->first.t << " " << count_new_diffs << "     \r";
			prev_t = p->first.t;
		}
		if (p->first.t > t_max + t_diff) //TODO: Consider other places for t_diff
			break;
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
			auto [count_tgt, deg_tgt, index_tgt] = level > kLevelMax / 2 ? CountTgt(basis_ss, deg, d_type, kLevelMax - level) : std::make_tuple(0, Deg(), -1);
			auto [count_src, deg_src, index_src] = CountSrc(basis_ss, deg, d_type, level);
			
			if (count_tgt == 0 && count_src == 1) { /* It is a target */
				int r = deg_src.v - deg.v;
				if (index == 0 || sc.levels[index - 1] < r) {
					++count_new_diffs;
					SetDiff(gb, basis, basis_ss, deg_src, d_type, GetRecentStaircase(basis_ss, deg_src).basis_ind[index_src], sc.basis_ind[index], deg_src.v - deg.v, t_test);
				}
			}
			else if (count_tgt == 1 && count_src == 0) { /* It is a source */
				Staircase sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
				++count_new_diffs;
				SetDiff(gb, basis, basis_ss, deg, d_type, sc.basis_ind[index], sc_tgt.basis_ind[index_tgt], deg.v - deg_tgt.v, t_test);
			}
			else if (count_tgt == 0 && count_src > 1) { /* It is a permanent cycle */
				if (sc.levels[index] > kLevelMax / 2) {
					++count_new_diffs;
					SetDiff(gb, basis, basis_ss, deg, d_type, sc.basis_ind[index], {}, kLevelMax / 4, t_test); //TODO: Implement SetUnknownImage
					AddImage(basis_ss, deg, d_type, sc.basis_ind[index], { -1 }, kLevelMax / 2 - 1); //
				}
			}
			else if ((level < kLevelMax / 2 || count_tgt == 0) && count_src == 0) /* It is a nontrivial permanent cycle */
				throw SSException(0x9eacf10fU, "Nontrivial permanant cycle");
		}
	}
	return count_new_diffs;
}

/* Check if there exists a global differential structure
 * Return 0 if there is no global structure
 * Return 1 if there exists a global structure 
 * Return 2 if the search space is too big
 * Return 3 if maximum number of iterations is reached */
int FindGlobalStructure(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, DiffType d_type, int t_try, int t_test, int iter_max)
{
	if (iter_max == 0)
		return 3; /* maximum number of iterations is reached */
	auto [s_diff, t_diff] = GetST(d_type);
	bool bForEt = IsEt(d_type);

	const size_t initial_history_length = basis_ss.size();
	CacheDeduction cache_history;
	cache_history.push_back();
	cache_history.InitNullDiffs(basis_ss, d_type, t_try);

	int count_iterations = 0;
	while (true) { /* Deep search for a global structure */
		++count_iterations;
		/* Find the place to start in cached history */
		Deg deg = cache_history.degs_.back();
		if (bForEt && cache_history.degs_.size() > 1 && deg.t != cache_history.degs_[cache_history.degs_.size() - 2].t) {
			if (!NoNullDiff(basis_ss, cache_history.degs_[cache_history.degs_.size() - 2].t)) { /* Invalid global structure at degree t. Roll back. */
				basis_ss.pop_back();
				cache_history.pop_back();
				continue;
			}
		}
		if (deg.s == -1) {
			basis_ss.resize(initial_history_length);
			return 1; /* Find one global structure */
		}
		if (iter_max != -1 && count_iterations >= iter_max) {
			basis_ss.resize(initial_history_length);
			return 3; /* maximum number of iterations is reached */
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
			const Deg deg_tgt = deg + Deg{ s_diff, t_diff, -r };
			const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
			for (int j : two_expansion((1 << nd.num_tgts) - 1 - i))
				tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)nd.first + j]);
		}

		basis_ss.push_back({}); /* Warning: this invalidates references sc above */
		cache_history.push_back();
		try {
			std::cout << deg << " " << nd.num_tgts << " d_{" << r << '}' << src << '=' << tgt << '\n';
			SetDiff(gb, basis, basis_ss, deg, d_type, std::move(src), std::move(tgt), r, t_test);
			cache_history.InitNullDiffs(basis_ss, d_type, t_try);
		}
		catch (SSException&) { /* Try the next possibility */
			std::cout << "Rollback\n";
			basis_ss.pop_back();
			cache_history.pop_back();
		}
	}
}

/* Return true if it cannot find a contradiction */
bool CheckCompatibility(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, DiffType d_type, int t_try, int t_test, int rCheck, int nCheck)
{
	if (nCheck == 0)
		return true;
	auto [s_diff, t_diff] = GetST(d_type);
	bool bForEt = IsEt(d_type);
	int count_new_diffs = 0;
	NullDiffs nullDiffs;
	nullDiffs.InitNullDiffs(basis_ss, d_type, t_try, true);

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
	array indices_degs = grbn::range((int)degs.size());
	std::stable_sort(indices_degs.begin(), indices_degs.end(), [nums_tgts](int i1, int i2) {return nums_tgts[i1] < nums_tgts[i2]; });

	int count_try = 0;
	for (int index_degs : indices_degs) {
		const Deg& deg = degs[index_degs];
		const Staircase& sc = GetRecentStaircase(basis_ss, deg);
		const NullDiff& nd = nullDiffs.null_diffs_.at(deg);
		if (nd.index == -1 || nd.num_tgts > 10 || sc.levels[nd.index] != kLevelMax - rCheck) // skip if there are too many possible targets
			continue;
		if (++count_try > nCheck)
			return true;

		const int r = kLevelMax - sc.levels[nd.index];
		const Deg deg_tgt = deg + Deg{ s_diff, t_diff, -r };
		array src = sc.basis_ind[nd.index];
		array tgt_pass;

		int count_pass = 0;
		for (int i = 0; i < (1 << nd.num_tgts); ++i) {
			/* Calculate tgt */
			array tgt;
			if (nd.num_tgts > 0) {
				const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
				for (int j : two_expansion((1 << nd.num_tgts) - 1 - i))
					tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)nd.first + j]);
			}

			if (i == (1 << nd.num_tgts) - 1 && count_pass == 0) {
				tgt_pass = tgt;
				++count_pass;
				break;
			}

			basis_ss.push_back({}); /* Warning: invalidate references sc above */
			bool bException = false;
			try {
				SetDiff(gb, basis, basis_ss, deg, d_type, src, tgt, r, t_test);
			}
			catch (SSException&) {
				bException = true;
			}
			basis_ss.pop_back();
			if (!bException) {
				tgt_pass = std::move(tgt);
				++count_pass;
				break;
			}
		}
		if (count_pass == 0)
			return false;
	}
	return true;
}

int DeduceDiffs(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, DiffType d_type, int t_try, int t_test, int max_poss, int nCheck, double secconds)
{
	auto [s_diff, t_diff] = GetST(d_type);
	bool bForEt = IsEt(d_type);
	Timer timer;
	int count_new_diffs = 0;

	int index_count = 0;
	for (auto p = basis_ss.front().begin(); p != basis_ss.front().end(); ++p) {
		const Deg& deg = p->first;
		std::cout << deg << "      \r";
		if (deg.t > t_try)
			break;

		int index = 0;
		while (true) {
			const Staircase& sc = GetRecentStaircase(basis_ss, deg);
			NullDiff nd = GetNextNullDiff(basis_ss, deg, d_type, index);
			if (nd.index == -1)
				break;
			else
				index = nd.index + 1;
			if (nd.num_tgts > max_poss) // skip if there are too many possible targets
				continue;

			if (sc.levels[nd.index] < kLevelMax / 2) {
				const int r = sc.levels[nd.index];
				const Deg deg_src = deg - Deg{ s_diff, t_diff, -r };
				array tgt = sc.basis_ind[nd.index];
				array src_pass;

				int count_pass = 0;
				for (int i = 0; i < (1 << nd.num_tgts); ++i) {
					/* Calculate tgt */
					array src;
					if (nd.num_tgts > 0) {
						const Staircase& sc_src = GetRecentStaircase(basis_ss, deg_src);
						for (int j : two_expansion((1 << nd.num_tgts) - 1 - i))
							src = lina::AddVectors(tgt, sc_src.basis_ind[(size_t)nd.first + j]);
					}

					//std::cout << "Try " << i + 1 << '/' << (1 << nd.num_tgts) << ": " << deg << " d_{" << r << '}' << src << '=' << tgt << '\n';

					if (i == (1 << nd.num_tgts) - 1 && count_pass == 0) {
						//std::cout << " Accepted\n";
						src_pass = src;
						++count_pass;
						break;
					}

					basis_ss.push_back({}); /* Warning: invalidate the reference sc above */
					bool bException = false;
					try {
						SetDiff(gb, basis, basis_ss, deg, d_type, src, tgt, r, t_test);
						/*if (!CheckCompatibility(gb, basis, basis_ss, d_type, r, t_try, t_test, nCheck))
							throw SSException(0xb07764bfU, "contradiction");*/
						if (!FindGlobalStructure(gb, basis, basis_ss, d_type, t_try, t_test, nCheck))
							throw SSException(0x890b3a9fU, "no global structure");
					}
					catch (SSException&) {
						bException = true;
					}
					basis_ss.pop_back();

					if (!bException) {
						src_pass = std::move(src);
						++count_pass;
						if (count_pass > 1)
							break;
					}
				}
				if (count_pass == 0)
					throw SSException(0x5b3d7e35U, "no compatible differentials");
				else if (count_pass == 1) {
					SetDiff(gb, basis, basis_ss, deg, d_type, src_pass, tgt, r);
					std::cout << deg << " " << nd.num_tgts << " d_{" << r << "}*" << src_pass << '=' << tgt << '\n';
					index = 0;
					++count_new_diffs;
					if (timer.Elapsed() > secconds) {
						timer.SuppressPrint();
						break;
					}
				}
			}
			else {
				const int r = kLevelMax - sc.levels[nd.index];
				const Deg deg_tgt = deg + Deg{ s_diff, t_diff, -r };
				array src = sc.basis_ind[nd.index];
				array tgt_pass;

				int count_pass = 0;
				for (int i = 0; i < (1 << nd.num_tgts); ++i) {
					/* Calculate tgt */
					array tgt;
					if (nd.num_tgts > 0) {
						const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
						for (int j : two_expansion((1 << nd.num_tgts) - 1 - i))
							tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)nd.first + j]);
					}

					if (i == (1 << nd.num_tgts) - 1 && count_pass == 0) {
						tgt_pass = tgt;
						++count_pass;
						break;
					}

					basis_ss.push_back({}); /* Warning: invalidate the reference sc above */
					bool bException = false;
					try {
						SetDiff(gb, basis, basis_ss, deg, d_type, src, tgt, r, t_test);
						/*if (!CheckCompatibility(gb, basis, basis_ss, d_type, r, t_try, t_test, nCheck))
							throw SSException(0xb07764bfU, "contradiction");*/
						if (!FindGlobalStructure(gb, basis, basis_ss, d_type, t_try, t_test, nCheck))
							throw SSException(0xb89c864eU, "no global structure");
					}
					catch (SSException&) {
						bException = true;
					}
					basis_ss.pop_back();

					if (!bException) {
						tgt_pass = std::move(tgt);
						++count_pass;
						if (count_pass > 1)
							break;
					}
				}
				if (count_pass == 0)
					throw SSException(0x66d5bd9aU, "no compatible differentials");
				else if (count_pass == 1) {
					SetDiff(gb, basis, basis_ss, deg, d_type, src, tgt_pass, r);
					std::cout << deg << " " << nd.num_tgts << " d_{" << r << '}' << src << '=*' << tgt_pass << '\n';
					index = 0;
					++count_new_diffs;
					if (timer.Elapsed() > secconds) {
						timer.SuppressPrint();
						break;
					}
				}
			}
		}
	}
	return count_new_diffs;
}

int DeduceDiffsChoosingSources(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, DiffType d_type, int t_try, int t_test, int nCheck)
{
	auto [s_diff, t_diff] = GetST(d_type);
	bool bForEt = IsEt(d_type);
	Timer timer;
	int count_new_diffs = 0;
	NullDiffs nullDiffs;
	nullDiffs.InitNullDiffs(basis_ss, d_type, t_try, true);

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
	array indices_degs = grbn::range((int)degs.size());
	std::stable_sort(indices_degs.begin(), indices_degs.end(), [nums_tgts](int i1, int i2) {return nums_tgts[i1] < nums_tgts[i2]; });

	int index_count = 0;
	for (int index_degs : indices_degs) {
		std::cout << ++index_count << '/' << indices_degs.size() << " " << count_new_diffs << "      \r";
		const Deg& deg = degs[index_degs];
		const Staircase& sc = GetRecentStaircase(basis_ss, deg);
		const NullDiff& nd = nullDiffs.null_diffs_.at(deg);
		if (nd.index == -1 || nd.num_tgts > 10) // skip if there are too many possible targets
			continue;
		const int r = kLevelMax - sc.levels[nd.index];
		const Deg deg_tgt = deg + Deg{ s_diff, t_diff, -r };
		array src = sc.basis_ind[nd.index];
		array tgt_pass;

		int count_pass = 0;
		for (int i = 0; i < (1 << nd.num_tgts); ++i) {
			/* Calculate tgt */
			array tgt;
			if (nd.num_tgts > 0) {
				const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
				for (int j : two_expansion((1 << nd.num_tgts) - 1 - i))
					tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)nd.first + j]);
			}

			//std::cout << "Try " << i + 1 << '/' << (1 << nd.num_tgts) << ": " << deg << " d_{" << r << '}' << src << '=' << tgt << '\n';

			if (i == (1 << nd.num_tgts) - 1 && count_pass == 0) {
				//std::cout << " Accepted\n";
				tgt_pass = tgt;
				++count_pass;
				break;
			}

			basis_ss.push_back({}); /* Warning: invalidate references sc above */
			bool bException = false;
			try {
				SetDiff(gb, basis, basis_ss, deg, d_type, src, tgt, r, t_test);
				if (!CheckCompatibility(gb, basis, basis_ss, d_type, r, t_try, t_test, nCheck))
					throw SSException(0xb07764bfU, "contradiction");
			}
			catch (SSException&) {
				bException = true;
			}
			basis_ss.pop_back();
			if (!bException) {
				//std::cout << " Success\n";
				tgt_pass = std::move(tgt);
				++count_pass;
				if (count_pass > 1)
					break;
			}
			//else
				//std::cout << " Fail\n";
		}
		if (count_pass == 0)
			throw SSException(0x5b3d7e35U, "no compatible differentials");
		else if (count_pass == 1) {
			basis_ss.push_back({});
			SetDiff(gb, basis, basis_ss, deg, d_type, src, tgt_pass, r);
			std::cout << deg << " " << nd.num_tgts << " d_{" << r << '}' << src << '=' << tgt_pass << '\n';
			if (count_new_diffs % 50 == 49) {
				DeduceZeroDiffs(gb, basis, basis_ss, d_type, -1, bForEt, false, -1);
				if (bForEt)
					DeduceDiffsForEt(gb, basis, basis_ss, d_type, -1);
			}
			nullDiffs.InitNullDiffs(basis_ss, d_type, t_try, false);
			ApplyRecentChanges(basis_ss);
			++count_new_diffs;

			if (timer.Elapsed() > 3600.0 * 1) { //
				timer.SuppressPrint();
				break;
			}
		}
	}
	return count_new_diffs;
}

int DeduceDiffsByNat(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, const Poly1d& gen_images, DiffType d_type,
	const std::vector<Deg>& gen_degs1, const grbn::GbWithCache& gb1, const std::map<Deg, Mon1d>& basis1, Staircases1d& basis_ss1, DiffType d_type1,
	int t_try, int t_test, bool bForEt)
{
	auto [s_diff, t_diff] = GetST(d_type);
	auto [s_diff1, t_diff1] = GetST(d_type1);
	Timer timer;
	int count_new_diffs = 0;
	NullDiffs nullDiffs;
	nullDiffs.InitNullDiffs(basis_ss, d_type, t_try, true);

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
	array indices_degs = grbn::range((int)degs.size());
	std::stable_sort(indices_degs.begin(), indices_degs.end(), [nums_tgts](int i1, int i2) {return nums_tgts[i1] < nums_tgts[i2]; });

	int index_count = 0;
	for (int index_degs : indices_degs) {
		std::cout << ++index_count << '/' << indices_degs.size() << " " << count_new_diffs << "      \r";
		const Deg& deg = degs[index_degs];
		const Staircase& sc = GetRecentStaircase(basis_ss, deg);
		const NullDiff& nd = nullDiffs.null_diffs_.at(deg);
		if (nd.index == -1 || nd.num_tgts > 10) // skip if there are no null diffs or a null diff has too many possible targets
			continue;
		const int r = kLevelMax - sc.levels[nd.index];
		const Deg deg_tgt = deg + Deg{ s_diff, t_diff, -r };
		array src = sc.basis_ind[nd.index];
		array tgt_pass;

		Poly poly_src = src.empty() ? Poly{} : indices_to_Poly(src, basis.at(deg));
		Poly poly_src_image = poly_src.empty() ? Poly{} : get_image(poly_src, gen_images, gb1);
		Deg deg1 = get_deg(poly_src_image, gen_degs1);
		array src_image = poly_src_image.empty() ? array{} : Poly_to_indices(poly_src_image, basis1.at(deg1));

		int count_pass = 0;
		for (int i = 0; i < (1 << nd.num_tgts); ++i) {
			/* Calculate tgt */
			array tgt;
			if (nd.num_tgts > 0) {
				const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
				for (int j : two_expansion((1 << nd.num_tgts) - 1 - i))
					tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)nd.first + j]);
			}

			Poly poly_tgt = tgt.empty() ? Poly{} : indices_to_Poly(tgt, basis.at(deg_tgt));
			Poly poly_tgt_image = poly_tgt.empty() ? Poly{} : get_image(poly_tgt, gen_images, gb1);
			Deg deg1_tgt = get_deg(poly_tgt_image, gen_degs1);
			array tgt_image = poly_tgt_image.empty() ? array{} : Poly_to_indices(poly_tgt_image, basis1.at(deg1_tgt));

			//std::cout << "Try " << i + 1 << '/' << (1 << nd.num_tgts) << ": " << deg << " d_{" << r << '}' << src << '=' << tgt;
			//std::cout << "-> " << deg1 << " d_{" << r << '}' << src_image1 << '=' << tgt_image1;

			if (i == (1 << nd.num_tgts) - 1 && count_pass == 0) {
				//std::cout << " Accepted\n";
				tgt_pass = tgt;
				++count_pass;
				break;
			}

			basis_ss.push_back({}); /* Warning: invalidate references sc above */
			basis_ss1.push_back({});
			bool bException = false;
			try {
				if (!src_image.empty() && !poly_tgt_image.empty() && deg1_tgt != deg1 + Deg{ s_diff1, t_diff1, -r })
					throw MyException(0xf434b447U, "degrees not compatible in E4t");
				SetDiffV2(gb1, basis1, basis_ss1, deg1, deg1_tgt, d_type1, src_image, tgt_image, r, -1);

				SetDiff(gb, basis, basis_ss, deg, d_type, src, tgt, r, t_test);
			}
			catch (SSException&) {
				bException = true;
			}
			basis_ss.pop_back();
			basis_ss1.pop_back();
			if (!bException) {
				//std::cout << " Success\n";
				tgt_pass = std::move(tgt);
				++count_pass;
				if (count_pass > 1)
					break;
			}
			/*else
				std::cout << " Fail\n";*/
		}
		if (count_pass == 0)
			throw SSException(0x5b3d7e35U, "no compatible differentials");
		else if (count_pass == 1) {
			basis_ss.push_back({});
			SetDiff(gb, basis, basis_ss, deg, d_type, src, tgt_pass, r);
			std::cout << deg << " " << nd.num_tgts << " d_{" << r << '}' << src << '=' << tgt_pass << '\n';
			/*if (count_new_diffs % 50 == 0) {
				DeduceZeroDiffs(gb, basis, basis_ss, d_type, -1, bForEt, false, -1);
				if (bForEt)
					DeduceDiffsForEt(gb, basis, basis_ss, d_type, -1);
			}*/
			nullDiffs.InitNullDiffs(basis_ss, d_type, t_try, false);
			ApplyRecentChanges(basis_ss);
			++count_new_diffs;

			if (timer.Elapsed() > 3600.0 * 7) { //
				timer.SuppressPrint();
				break;
			}
		}
	}
	return count_new_diffs;
}

/* For a map SS1->SS2, add induced differentials from SS1 to SS2 */
int DeduceDiffsByNatV2(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, const Poly1d& gen_to_E4t, DiffType d_type,
	const std::vector<Deg>& gen_degs1, const grbn::GbWithCache& gb1, const std::map<Deg, Mon1d>& basis1, Staircases1d& basis_ss1, DiffType d_type1)
{
	auto [s_diff, t_diff] = GetST(d_type);
	auto [s_diff1, t_diff1] = GetST(d_type1);
	Timer timer;

	int index_count = 0;
	for (auto p = basis_ss.front().begin(); p != basis_ss.front().end(); ++p) {
		const Deg& deg = p->first;
		std::cout << deg << "      \r";
		const Staircase& sc = GetRecentStaircase(basis_ss, deg);
		for (size_t i = 0; i < sc.levels.size(); ++i) {
			int r = kLevelMax - sc.levels[i];
			const Deg deg_tgt = deg + Deg{ s_diff, t_diff, -r };
			const array& src = sc.basis_ind[i];
			array tgt;
			if (r < 5000) {
				if (sc.diffs_ind[i] == array{ -1 }) {
					tgt = array{};
					r -= 2;
					if (r < kLevelMin)
						continue;
				}
				else {
					tgt = sc.diffs_ind[i];
				}
			}
			else {
				if (sc.diffs_ind[i] == array{ -1 }) {
					tgt = array{};
					r = kLevelMax / 4;
				}
				else {
					continue;
				}
			}

			Poly poly_src = src.empty() ? Poly{} : indices_to_Poly(src, basis.at(deg));
			Poly poly_src_image = poly_src.empty() ? Poly{} : get_image(poly_src, gen_to_E4t, gb1);
			Deg deg1 = get_deg(poly_src_image, gen_degs1);
			array src_image = poly_src_image.empty() ? array{} : Poly_to_indices(poly_src_image, basis1.at(deg1));

			Poly poly_tgt = tgt.empty() ? Poly{} : indices_to_Poly(tgt, basis.at(deg_tgt));
			Poly poly_tgt_image = poly_tgt.empty() ? Poly{} : get_image(poly_tgt, gen_to_E4t, gb1);
			Deg deg1_tgt = get_deg(poly_tgt_image, gen_degs1);
			array tgt_image = poly_tgt_image.empty() ? array{} : Poly_to_indices(poly_tgt_image, basis1.at(deg1_tgt));
			
			int nNewDiffs = SetDiffV2(gb1, basis1, basis_ss1, deg1, deg1_tgt, d_type1, src_image, tgt_image, r, -1);
			index_count += nNewDiffs;
			if (nNewDiffs)
				std::cout << deg << ' ' << deg1 << ' ' << deg1_tgt << ' ' << " d_{" << r << '}' << src_image << '=' << tgt_image << '\n';
		}
	}
	return index_count;
}

int DeduceDiffsForE4(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, const Poly1d& gen_loc, const Poly1d& gen_loc_sq0, const Poly1d& gen_to_E4t, DiffType d_type,
	const std::vector<Deg>& gen_degs_E4h0, const grbn::GbWithCache& gb_E4h0, const std::map<Deg, Mon1d>& basis_E4h0, Staircases1d& basis_ss_E4h0, DiffType d_type_E4h0,
	const std::vector<Deg>& gen_degs_E4t, const grbn::GbWithCache& gb_E4t, const std::map<Deg, Mon1d>& basis_E4t, Staircases1d& basis_ss_E4t, DiffType d_type_E4t,
	int t_try, int t_test)
{
	auto [s_diff, t_diff] = GetST(d_type);
	bool bForEt = false;
	Timer timer;
	int count_new_diffs = 0;
	NullDiffs nullDiffs;
	nullDiffs.InitNullDiffs(basis_ss, d_type, t_try, true);

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
	array indices_degs = grbn::range((int)degs.size());
	std::stable_sort(indices_degs.begin(), indices_degs.end(), [nums_tgts](int i1, int i2) {return nums_tgts[i1] < nums_tgts[i2]; });

	int index_count = 0;
	for (int index_degs : indices_degs) {
		std::cout << ++index_count << '/' << indices_degs.size() << " " << count_new_diffs << "      \r";
		const Deg& deg = degs[index_degs];
		const Staircase& sc = GetRecentStaircase(basis_ss, deg);
		const NullDiff& nd = nullDiffs.null_diffs_.at(deg);
		if (nd.index == -1 || nd.num_tgts > 10) // skip if there are no null diffs or a null diff has too many possible targets
			continue;
		const int r = kLevelMax - sc.levels[nd.index];
		const Deg deg_tgt = deg + Deg{ s_diff, t_diff, -r };
		array src = sc.basis_ind[nd.index];
		array tgt_pass;

		Poly poly_src = src.empty() ? Poly{} : indices_to_Poly(src, basis.at(deg));
		Poly poly_src_image1 = poly_src.empty() ? Poly{} : get_image(poly_src, gen_loc, gb_E4h0);
		Deg deg1 = get_deg(poly_src_image1, gen_degs_E4h0);
		array src_image1 = poly_src_image1.empty() ? array{} : Poly_to_indices(poly_src_image1, basis_E4h0.at(deg1));

		Poly poly_src_image2 = poly_src.empty() ? Poly{} : get_image(poly_src, gen_loc_sq0, gb_E4t);
		Deg deg2 = get_deg(poly_src_image2, gen_degs_E4t);
		array src_image2 = poly_src_image2.empty() ? array{} : Poly_to_indices(poly_src_image2, basis_E4t.at(deg2));

		Poly poly_src_image3;
		Deg deg3;
		array src_image3;
		if (deg.t <= 189) {
			poly_src_image3 = poly_src.empty() ? Poly{} : get_image(poly_src, gen_to_E4t, gb_E4t);
			deg3 = get_deg(poly_src_image3, gen_degs_E4t);
			src_image3 = poly_src_image3.empty() ? array{} : Poly_to_indices(poly_src_image3, basis_E4t.at(deg3));
		}

		int count_pass = 0;
		for (int i = 0; i < (1 << nd.num_tgts); ++i) {
			/* Calculate tgt */
			array tgt;
			if (nd.num_tgts > 0) {
				const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
				for (int j : two_expansion((1 << nd.num_tgts) - 1 - i))
					tgt = lina::AddVectors(tgt, sc_tgt.basis_ind[(size_t)nd.first + j]);
			}

			Poly poly_tgt = tgt.empty() ? Poly{} : indices_to_Poly(tgt, basis.at(deg_tgt));
			Poly poly_tgt_image1 = poly_tgt.empty() ? Poly{} : get_image(poly_tgt, gen_loc, gb_E4h0);
			Deg deg_tgt1 = get_deg(poly_tgt_image1, gen_degs_E4h0);
			array tgt_image1 = poly_tgt_image1.empty() ? array{} : Poly_to_indices(poly_tgt_image1, basis_E4h0.at(deg_tgt1));

			Poly poly_tgt_image2 = poly_tgt.empty() ? Poly{} : get_image(poly_tgt, gen_loc_sq0, gb_E4t);
			Deg deg_tgt2 = get_deg(poly_tgt_image2, gen_degs_E4t);
			array tgt_image2 = poly_tgt_image2.empty() ? array{} : Poly_to_indices(poly_tgt_image2, basis_E4t.at(deg_tgt2));

			Poly poly_tgt_image3;
			Deg deg_tgt3;
			array tgt_image3;
			if (deg.t <= 189) {
				poly_tgt_image3 = poly_tgt.empty() ? Poly{} : get_image(poly_tgt, gen_to_E4t, gb_E4t);
				deg_tgt3 = get_deg(poly_tgt_image3, gen_degs_E4t);
				tgt_image3 = poly_tgt_image3.empty() ? array{} : Poly_to_indices(poly_tgt_image3, basis_E4t.at(deg_tgt3));
			}

			//std::cout << "Try " << i + 1 << '/' << (1 << nd.num_tgts) << ": " << deg << " d_{" << r << '}' << src << '=' << tgt;
			//std::cout << "-> " << deg1 << " d_{" << r << '}' << src_image1 << '=' << tgt_image1;

			if (i == (1 << nd.num_tgts) - 1 && count_pass == 0) {
				//std::cout << " Accepted\n";
				tgt_pass = tgt;
				++count_pass;
				break;
			}

			basis_ss.push_back({}); /* Warning: invalidate references sc above */
			basis_ss_E4h0.push_back({});
			basis_ss_E4t.push_back({});
			bool bException = false;
			try {
				//if (!src_image1.empty() && !poly_tgt_image1.empty() && deg_tgt1 != deg1 + Deg{ s_diff_E4h0, t_diff_E4h0, -r })
					//throw MyException(0x5f3f1bc7U, "degrees not compatible in E4t");
				SetDiffV2(gb_E4h0, basis_E4h0, basis_ss_E4h0, deg1, deg_tgt1, d_type_E4h0, src_image1, tgt_image1, r, -1);
				SetDiffV2(gb_E4t, basis_E4t, basis_ss_E4t, deg2, deg_tgt2, d_type_E4t, src_image2, tgt_image2, r, -1);
				if (deg.t <= 189)
					SetDiffV2(gb_E4t, basis_E4t, basis_ss_E4t, deg3, deg_tgt3, d_type_E4t, src_image3, tgt_image3, r, -1);

				SetDiff(gb, basis, basis_ss, deg, d_type, src, tgt, r, t_test);
			}
			catch (SSException&) {
				bException = true;
			}
			basis_ss.pop_back();
			basis_ss_E4h0.pop_back();
			basis_ss_E4t.pop_back();
			if (!bException) {
				//std::cout << " Success\n";
				tgt_pass = std::move(tgt);
				++count_pass;
				if (count_pass > 1)
					break;
			}
			/*else
				std::cout << " Fail\n";*/
		}
		if (count_pass == 0)
			throw SSException(0x5b3d7e35U, "no compatible differentials");
		else if (count_pass == 1) {
			basis_ss.push_back({});
			SetDiff(gb, basis, basis_ss, deg, d_type, src, tgt_pass, r);
			std::cout << deg << " " << nd.num_tgts << " d_{" << r << '}' << src << '=' << tgt_pass << '\n';
			/*if (count_new_diffs % 50 == 0) {
				DeduceZeroDiffs(gb, basis, basis_ss, d_type, -1, bForEt, false, -1);
				if (bForEt)
					DeduceDiffsForEt(gb, basis, basis_ss, d_type, -1);
			}*/
			nullDiffs.InitNullDiffs(basis_ss, d_type, t_try, false);
			ApplyRecentChanges(basis_ss);
			++count_new_diffs;

			if (timer.Elapsed() > 3600.0 * 1) { //
				timer.SuppressPrint();
				break;
			}
		}
	}
	return count_new_diffs;
}

void WrapDeduceZeroDiffs(const Database& db, const std::string& table_prefix, int t_max)
{
	DiffType d_type = GetDiffType(table_prefix);
	bool bForEt = IsEt(table_prefix);
	Timer timer;
	const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", -1);
	const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", -1);
	Staircases1d basis_ss = { db.load_basis_ss(table_prefix + "_ss", -1), {} };

	int count = 0;
	try {
		count = DeduceZeroDiffs(gb, basis, basis_ss, d_type, t_max, bForEt, true, -1);

		db.begin_transaction();
		db.update_basis_ss(table_prefix + "_ss", basis_ss[1]);
		db.end_transaction();
	}
	catch (SSException& e) {
		std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
	}
	catch (MyException& e) {
		std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
	}

	std::cout << "Changed differentials: " << count << '\n';
}

void WrapDeduceDiffsForEt(const Database& db, const std::string& table_prefix, int t_max)
{
	try {
		if (!IsEt(table_prefix))
			throw MyException(0x2f64a754U, table_prefix + " is not an Et type table");
		DiffType d_type = GetDiffType(table_prefix);
		Timer timer;
		const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", -1);
		const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", -1);
		Staircases1d basis_ss = { db.load_basis_ss(table_prefix + "_ss", -1), {} };

		int count = 0;
		count = DeduceDiffsForEt(gb, basis, basis_ss, d_type, t_max, true);

		db.begin_transaction();
		db.update_basis_ss(table_prefix + "_ss", basis_ss[1]);
		db.end_transaction();

		std::cout << "Changed differentials: " << count << '\n';
	}
	catch (SSException& e) {
		std::cerr << "Error " << std::hex << e.id() << ": " << e.what() << '\n';
	}
	catch (MyException& e) {
		std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
	}
}

void WrapDeduceDiffs(const Database& db, const std::string& table_prefix, int t_try, int t_test, int max_poss, int nCheck, double seconds)
{
	DiffType d_type = GetDiffType(table_prefix);
	bool bForEt = IsEt(table_prefix);
	Timer timer;
	const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", -1);
	const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", -1);
	Staircases1d basis_ss = { db.load_basis_ss(table_prefix + "_ss", -1), {} };

	int count = 0;
	try {
		count = DeduceDiffs(gb, basis, basis_ss, d_type, t_try, t_test, max_poss, nCheck, seconds);

		db.begin_transaction();
		db.update_basis_ss(table_prefix + "_ss", basis_ss[1]);
		db.end_transaction();
	}
	catch (SSException& e) {
		std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
	}
	catch (MyException& e) {
		std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
	}
	/*catch (NoException&) {
		;
	}*/

	std::cout << "Changed differentials: " << count << '\n';
}

void WrapDeduceDiffsByNat(const Database& db, const std::string& table_prefix, const std::string& table_prefix1, std::string& column,
	int t_try, int t_test)
{
	try {
		DiffType d_type = GetDiffType(table_prefix);
		DiffType d_type1 = GetDiffType(table_prefix1);
		bool bForEt = IsEt(table_prefix);
		Timer timer;
		const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", -1);
		const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", -1);
		Staircases1d basis_ss = { db.load_basis_ss(table_prefix + "_ss", -1), {} };

		Poly1d gen_to_E4t = db.load_gen_images(table_prefix + "_generators", column, column == "to_E4t" ? 189 : 230);

		const std::vector<Deg> gen_degs1 = db.load_gen_degs(table_prefix1 + "_generators");
		const grbn::GbWithCache gb1 = db.load_gb(table_prefix1 + "_relations", t_test + 1);
		const std::map<Deg, Mon1d> basis1 = db.load_basis(table_prefix1 + "_basis", t_test + 1);
		Staircases1d basis_ss1 = { db.load_basis_ss(table_prefix1 + "_ss", t_test + 1), {} };

		int count = 0;
		count = DeduceDiffsByNat(gb, basis, basis_ss, gen_to_E4t, d_type, gen_degs1, gb1, basis1, basis_ss1, d_type1, t_try, t_test, bForEt);

		db.begin_transaction();
		db.update_basis_ss(table_prefix + "_ss", basis_ss[1]);
		db.end_transaction();

		std::cout << "Changed differentials: " << count << '\n';
	}
	catch (SSException& e) {
		std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
	}
	catch (MyException& e) {
		std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
	}
	/*catch (NoException&) {
		;
	}*/
}

void WrapDeduceDiffsByNatV2(const Database& db, const std::string& table_prefix, const std::string& table_prefix1, std::string& column, int t_max)
{
	try {
		DiffType d_type = GetDiffType(table_prefix);
		DiffType d_type1 = GetDiffType(table_prefix1);
		Timer timer;
		const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", t_max);
		const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", t_max);
		Staircases1d basis_ss = { db.load_basis_ss(table_prefix + "_ss", t_max), {} };

		Poly1d gen_to_E4t = db.load_gen_images(table_prefix + "_generators", column, column == "to_E4t" ? 189 : 230);

		const std::vector<Deg> gen_degs1 = db.load_gen_degs(table_prefix1 + "_generators");
		const grbn::GbWithCache gb1 = db.load_gb(table_prefix1 + "_relations", -1);
		const std::map<Deg, Mon1d> basis1 = db.load_basis(table_prefix1 + "_basis", -1);
		Staircases1d basis_ss1 = { db.load_basis_ss(table_prefix1 + "_ss", -1), {} };

		int count = 0;
		count = DeduceDiffsByNatV2(gb, basis, basis_ss, gen_to_E4t, d_type, gen_degs1, gb1, basis1, basis_ss1, d_type1);

		db.begin_transaction();
		db.update_basis_ss(table_prefix1 + "_ss", basis_ss1[1]);
		db.end_transaction();

		std::cout << "Changed differentials: " << count << '\n';
	}
	catch (SSException& e) {
		std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
	} 
	catch (MyException& e) {
		std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
	}
	/*catch (NoException&) {
		;
	}*/
}

void WrapDeduceDiffsForE4(const Database& db, int t_try, int t_test)
{
	std::string table_prefix = "E4";
	std::string table_prefix1 = "E4h0";
	std::string table_prefix2 = "E4t";

	std::string column1 = "loc";
	std::string column2 = "loc_sq0";
	std::string column3 = "to_E4t";
	try {
		DiffType d_type = GetDiffType(table_prefix);
		DiffType d_type1 = GetDiffType(table_prefix1);
		DiffType d_type2 = GetDiffType(table_prefix2);
		Timer timer;
		const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", -1);
		const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", -1);
		Staircases1d basis_ss = { db.load_basis_ss(table_prefix + "_ss", -1), {} };

		Poly1d gen_loc = db.load_gen_images(table_prefix + "_generators", column1, 230);
		Poly1d gen_loc_sq0 = db.load_gen_images(table_prefix + "_generators", column2, 230);
		Poly1d gen_to_E4t = db.load_gen_images(table_prefix + "_generators", column3, 189);

		const std::vector<Deg> gen_degs1 = db.load_gen_degs(table_prefix1 + "_generators");
		const grbn::GbWithCache gb1 = db.load_gb(table_prefix1 + "_relations", t_test + 1);
		const std::map<Deg, Mon1d> basis1 = db.load_basis(table_prefix1 + "_basis", t_test + 1);
		Staircases1d basis_ss1 = { db.load_basis_ss(table_prefix1 + "_ss", t_test + 1), {} };

		const std::vector<Deg> gen_degs2 = db.load_gen_degs(table_prefix2 + "_generators");
		const grbn::GbWithCache gb2 = db.load_gb(table_prefix2 + "_relations", t_test + 1);
		const std::map<Deg, Mon1d> basis2 = db.load_basis(table_prefix2 + "_basis", t_test + 1);
		Staircases1d basis_ss2 = { db.load_basis_ss(table_prefix2 + "_ss", t_test + 1), {} };

		int count = 0;
		count = DeduceDiffsForE4(gb, basis, basis_ss, gen_loc, gen_loc_sq0, gen_to_E4t, d_type,
			gen_degs1, gb1, basis1, basis_ss1, d_type1,
			gen_degs2, gb2, basis2, basis_ss2, d_type2,
			t_try, t_test);

		db.begin_transaction();
		db.update_basis_ss(table_prefix + "_ss", basis_ss[1]);
		db.end_transaction();

		std::cout << "Changed differentials: " << count << '\n';
	}
	/*catch (SSException& e) {
		std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
	}
	catch (MyException& e) {
		std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
	}*/
	catch (NoException&) {
		;
	}
}

/* generate the table of the spectral sequence */
void DeduceDiffs(const Database& db, const std::string& table_prefix, int t_max)
{
	try {
		DiffType d_type = GetDiffType(table_prefix);
		bool bForEt = IsEt(table_prefix);
		const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", t_max);
		const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", t_max);
		Staircases1d basis_ss = { db.load_basis_ss(table_prefix + "_ss", t_max) };

		while (true) {
			int count = DeduceDiffs(gb, basis, basis_ss, d_type, -1, -1, bForEt);
			std::cout << "DeduceDiffs: " << count << '\n';

			if (count == 0)
				break;
		}
		/*while (true) {
			int count1 = DeduceZeroDiffs(gb, basis, basis_ss, true);
			std::cout << "DeduceZeroDiffs: " << count1 << '\n';

			int count2 = DeduceDiffsForEt(gb, basis, basis_ss, -1);
			std::cout << "DeduceDiffsForEt: " << count2 << '\n';

			int count3 = DeduceDiffs(gb, basis, basis_ss, -1, -1, 1000, true);
			std::cout << "DeduceDiffs: " << count3 << '\n';

			if (count1 + count2 + count3 == 0)
				break;
		}
		while (true) {
			int count1 = DeduceZeroDiffs(gb, basis, basis_ss, true);
			std::cout << "DeduceZeroDiffs: " << count1 << '\n';

			int count2 = DeduceDiffsForEt(gb, basis, basis_ss, -1);
			std::cout << "DeduceDiffsForEt: " << count2 << '\n';

			int count3 = DeduceDiffs(gb, basis, basis_ss, -1, -1, 2000, true);
			std::cout << "DeduceDiffs: " << count3 << '\n';

			if (count1 + count2 + count3 == 0)
				break;
		}*/
		/* insert into the database */
		db.begin_transaction();
		//db.execute_cmd("DELETE FROM " + table_prefix + "_ss_tmp;");
		db.update_basis_ss(table_prefix + "_ss", basis_ss.front());
		db.end_transaction();
	}
	catch (SSException& e)
	{
		std::cerr << "Error code " << std::hex << e.id() << ": " << e.what() << '\n';
	}
	catch (MyException& e) {
		std::cerr << "MyError " << std::hex << e.id() << ": " << e.what() << '\n';
	}
}
