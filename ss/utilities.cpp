#include "main.h"
#include "benchmark.h"
#include "linalg.h"
#include "myio.h"
#include <iostream>
#include <sstream>


std::ostream& operator<<(std::ostream& sout, const Staircase& sc)
{
	sout << "Staircase:\n";
	for (size_t i = 0; i < sc.levels.size(); ++i) {
		sout << sc.levels[i] << "  " << sc.basis_ind[i] << "  " << sc.diffs_ind[i] << '\n';
	}
	return sout;
}

/* Return the newest version of the staircase in history */
const Staircase& GetRecentStaircase(const Staircases1d& basis_ss, const Deg& deg)
{
	for (auto p = basis_ss.rbegin(); p != basis_ss.rend(); ++p)
		if (p->find(deg) != p->end())
			return p->at(deg);
	throw MyException(0x553989e0U, "BUG: deg not found");
}

void ApplyAllChanges(Staircases1d& basis_ss, size_t nHistory /*= 0*/)
{
	for (auto p = basis_ss[nHistory].begin(); p != basis_ss[nHistory].end(); ++p) {
		auto& sc = GetRecentStaircase(basis_ss, p->first);
		if (&(p->second) != &sc)
			p->second = std::move(sc);
	}
	basis_ss.resize(nHistory + 1);
}

void ApplyRecentChanges(Staircases1d& basis_ss)
{
	size_t index_before_last = basis_ss.size() - 2;
	for (auto p = basis_ss.back().begin(); p != basis_ss.back().end(); ++p)
		basis_ss[index_before_last][p->first] = std::move(p->second);
	basis_ss.pop_back();
}

size_t GetFirstIndexOfLevel(const Staircase& sc, int level)
{
	return std::lower_bound(sc.levels.begin(), sc.levels.end(), level) - sc.levels.begin();
}

/* Return the first index of a level such that all levels above (inclusive) have known differentials */
size_t GetFirstIndexOfKnownLevels(const Staircase& sc, int level)
{
	size_t result = sc.levels.size();
	for (size_t i = sc.levels.size(); i-- > 0; ) {
		if (sc.diffs_ind[i] == array{ -1 } || sc.levels[i] < level)
			break;
		if (i == 0 || sc.levels[i - 1] != sc.levels[i])
			result = i;
	}
	return result;
}

/* Return if there is no null differential in degree t */
bool NoNullDiff(const Staircases1d& basis_ss, int t)
{
	auto p1 = basis_ss.front().lower_bound(Deg{ 0, t, 0 });
	auto p2 = basis_ss.front().lower_bound(Deg{ 0, t + 1, 0 });
	for (auto p = p1; p != p2; ++p) {
		const Staircase& sc = GetRecentStaircase(basis_ss, p->first);
		for (size_t i = 0; i < sc.diffs_ind.size(); ++i)
			if (sc.diffs_ind[i] == array{ -1 })
				return false;
	}
	return true;
}

/* Count the number of all possible d_r targets.
 * Return (count, index). */
std::pair<int, int> CountDrTgt(const Staircases1d& basis_ss, const Deg& deg_tgt, int r)
{
	std::pair<int, int> result;
	if (basis_ss.front().find(deg_tgt) != basis_ss.front().end()) {
		const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
		result.second = (int)GetFirstIndexOfLevel(sc_tgt, r);
		result.first = (int)GetFirstIndexOfKnownLevels(sc_tgt, kLevelMax - r + 2) - result.second;
	}
	else
		result = { 0, -1 };
	return result;
}

/* Count the number of all possible targets.
 * Return (count, deg, index).
 * count>=2 means there are >= 2 elements in the range and 
 * deg corresponds to the one with the least r. */
std::tuple<int, Deg, int> CountTgt(const Staircases1d& basis_ss, const Deg& deg, int r)
{
	Deg deg_tgt; int count = 0, index = -1;
	const Staircase& sc = GetRecentStaircase(basis_ss, deg);
	for (int r1 = r; r1 <= deg.v; r1 += 2) {
		Deg d_tgt = deg + Deg{ 1, 0, -r1 };
		auto [c, first_r1] = CountDrTgt(basis_ss, d_tgt, r1);
		if (c > 0) {
			if (count == 0) {
				deg_tgt = d_tgt;
				index = first_r1;
			}
			count += c;
			if (count >= 2)
				break;
		}
	}
	return std::make_tuple(count, deg_tgt, index);
}

/* Count the number of all possible sources
 * Return (Deg, count, index)
 * count>=2 means there are >= 2 elements in the range and 
 * deg corresponds to the one with the least r. */
std::tuple<int, Deg, int> CountSrc(const Staircases1d& basis_ss, const Deg& deg, int level)
{
	Deg deg_src; int count = 0, index = -1;
	int r_max = std::min(level, 2 * (deg.t - deg.s + 1) - deg.v);
	const Staircase& sc = GetRecentStaircase(basis_ss, deg);
	for (int r1 = kLevelMin; r1 <= r_max; r1 += 2) {
		Deg d_src = deg - Deg{ 1, 0, -r1 };
		if (basis_ss.front().find(d_src) != basis_ss.front().end()) {
			const Staircase& sc_src = GetRecentStaircase(basis_ss, d_src);
			int first_Nmr1 = (int)GetFirstIndexOfLevel(sc_src, kLevelMax - r1);
			int c = (int)GetFirstIndexOfKnownLevels(sc_src, kLevelMax - r1 + 2) - first_Nmr1;
			if (c > 0) {
				if (count == 0) {
					deg_src = d_src;
					index = first_Nmr1;
				}
				count += c;
				if (count >= 2)
					break;
			}
		}
	}
	return std::make_tuple(count, deg_src, index);
}

NullDiff GetNullDiff(const Staircases1d& basis_ss, const Deg& deg)
{
	NullDiff result{ -1, -1, -1 };
	const Staircase& sc = GetRecentStaircase(basis_ss, deg);
	for (size_t i = sc.diffs_ind.size(); i-- > 0;) {
		if (sc.diffs_ind[i] == array{ -1 } && sc.levels[i] > kLevelMax / 4 * 3) {
			int r = kLevelMax - sc.levels[i];
			Deg deg_tgt = deg + Deg{ 1, 0, -r };
			auto [num_tgts, first_r] = CountDrTgt(basis_ss, deg_tgt, r);
			result = { (int)i, num_tgts, first_r };
			break;
		}
	}
	return result;
}

void NullDiffs::InitNullDiffs(const Staircases1d& basis_ss, int t_max, bool bNew)
{
	const std::map<Deg, Staircase>& basis_ss0 = bNew ? basis_ss.front() : basis_ss.back();
	if (t_max == -1)
		t_max = INT_MAX;
	for (auto p = basis_ss0.begin(); p != basis_ss0.end(); ++p) {
		if (p->first.t > t_max)
			break;
		null_diffs_[p->first] = GetNullDiff(basis_ss, p->first);
	}
}

const NullDiff& CacheDeduction::GetRecentNullDiff(const Deg& deg) const
{
	for (auto p = null_diffs_.rbegin(); p != null_diffs_.rend(); ++p)
		if (p->find(deg) != p->end())
			return p->at(deg);
	throw MyException(0xe7610d58U, "BUG: deg not found");
}

/* Cache the data about null diffs. Cache the deg of the least targets. */
void CacheDeduction::InitNullDiffs(const Staircases1d& basis_ss, int t_max)
{
	bool bNew = null_diffs_.size() == 1;
	const std::map<Deg, Staircase>& basis_ss0 = bNew ? basis_ss.front() : basis_ss.back();
	if (t_max == -1)
		t_max = INT_MAX;

	/* Initialize null_diffs_.back() */
	for (auto p = basis_ss0.begin(); p != basis_ss0.end(); ++p) {
		if (p->first.t > t_max)
			break;
		null_diffs_.back()[p->first] = GetNullDiff(basis_ss, p->first);
	}

	/* Initialize degs_.back() */
	int min_num_tgts = INT_MAX;
	for (auto p = null_diffs_.front().begin(); p != null_diffs_.front().end(); ++p) {
		if (p->first.t > t_max)
			break;
		int num_tgts = GetRecentNullDiff(p->first).num_tgts;
		if (num_tgts >= 0 && (degs_.back().t == -1 || degs_.back().t > p->first.t || (degs_.back().t == p->first.t && min_num_tgts > num_tgts))) {
			degs_.back() = p->first;
			min_num_tgts = num_tgts;
		}
	}
}

/* Add x, dx, level and triangularize.
 * Output the image of a differential that should be moved to the next level */
Staircase triangularize(const Staircase& sc, size_t i_insert, array x, array dx, int level, array& image, int& level_image)
{
	level_image = -1;
	Staircase sc_new;

	/* Copy the first part */
	sc_new.basis_ind.insert(sc_new.basis_ind.end(), sc.basis_ind.begin(), sc.basis_ind.begin() + i_insert);
	sc_new.diffs_ind.insert(sc_new.diffs_ind.end(), sc.diffs_ind.begin(), sc.diffs_ind.begin() + i_insert);
	sc_new.levels.insert(sc_new.levels.end(), sc.levels.begin(), sc.levels.begin() + i_insert);

	/* Add x, dx, level at the insertion point */
	sc_new.basis_ind.push_back(std::move(x));
	sc_new.diffs_ind.push_back(std::move(dx));
	sc_new.levels.push_back(level);

	/* Triangularize the rest */
	for (size_t i = i_insert; i < sc.basis_ind.size(); ++i) {
		array base_ind = sc.basis_ind[i];
		array diff_ind = sc.diffs_ind[i];
		for (size_t j = i_insert; j < sc_new.basis_ind.size(); ++j) {
			if (std::binary_search(base_ind.begin(), base_ind.end(), sc_new.basis_ind[j][0])) {
				base_ind = lina::AddVectors(base_ind, sc_new.basis_ind[j]);
				if (sc.levels[i] == sc_new.levels[j] && diff_ind != array{ -1 })
					diff_ind = lina::AddVectors(diff_ind, sc_new.diffs_ind[j]);
			}
		}
		if (!base_ind.empty()) {
			sc_new.basis_ind.push_back(std::move(base_ind));
			sc_new.diffs_ind.push_back(std::move(diff_ind));
			sc_new.levels.push_back(sc.levels[i]);
		}
		else if (diff_ind != array{ -1 } && !diff_ind.empty()) {
			image = std::move(diff_ind);
			level_image = kLevelMax - sc.levels[i];
		}
	}
#ifdef _DEBUG
	if (sc_new.basis_ind.size() != sc.basis_ind.size())
		throw MyException(0xdd696043U, "BUG: triangularize()");
#endif
	return sc_new;
}

/* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
void AddDiff(Staircases1d& basis_ss, const Deg& deg_x, array x, array dx, int r)
{
	Deg deg_dx = deg_x + Deg{ 1, 0, -r };

	/* If x is zero then dx is in Im(d_{r-2}) */
	if (x.empty()) {
		if (dx != array{ -1 } && !dx.empty())
			AddImage(basis_ss, deg_dx, std::move(dx), { -1 }, r - 2);
		return;
	}

	/* If x is in Ker(d_r) then dx is in Im(d_{r-2}) */
	const Staircase& sc = GetRecentStaircase(basis_ss, deg_x);
	size_t first_Nmr = GetFirstIndexOfLevel(sc, kLevelMax - r);
	x = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_Nmr, x);
	if (x.empty()) {
		if (dx != array{ -1 } && !dx.empty())
			AddImage(basis_ss, deg_dx, std::move(dx), { -1 }, r - 2);
		return;
	}

	array image_new;
	int level_image_new = -1;
	if (dx == array{ -1 }) { /* If the target is unknown, insert it to the end of level N-r. */
		size_t first_Nmrp2 = GetFirstIndexOfLevel(sc, kLevelMax - r + 2);
		x = lina::Residue(sc.basis_ind.begin() + first_Nmr, sc.basis_ind.begin() + first_Nmrp2, x);
		if (!x.empty()) {
			basis_ss.back()[deg_x] = triangularize(sc, first_Nmrp2, x, { -1 }, kLevelMax - r, image_new, level_image_new);
		}
	}
	else if (dx.empty()) { /* If the target is zero, insert it to the end of level N-r-2 */
		basis_ss.back()[deg_x] = triangularize(sc, first_Nmr, x, { -1 }, kLevelMax - r - 2, image_new, level_image_new);
	}
	else { /* Otherwise insert it to the beginning of level N-r */
		basis_ss.back()[deg_x] = triangularize(sc, first_Nmr, x, dx, kLevelMax - r, image_new, level_image_new); //
	}

	if (level_image_new != -1) {
		if (level_image_new < kLevelMax / 2) { /* Add a d_{r1-2} image */
			Deg deg_image_new = deg_x + Deg{ 1, 0, -level_image_new };
			AddImage(basis_ss, deg_image_new, std::move(image_new), { -1 }, level_image_new - 2);
		}
		else { /* Add a d_r1 cycle */
			int r_image = kLevelMax - level_image_new;
			Deg deg_image_new = deg_x - Deg{ 1, 0, -r_image };
			AddDiff(basis_ss, deg_image_new, std::move(image_new), { -1 }, r_image + 2);
		}
	}

	/* Add image */
	if (dx != array{ -1 } && !dx.empty())
		AddImage(basis_ss, deg_dx, std::move(dx), std::move(x), r);
}

/* Add an image of d_r.
 * Assume dx is nonempty. */
void AddImage(Staircases1d& basis_ss, const Deg& deg_dx, array dx, array x, int r)
{
	Deg deg_x = deg_dx - Deg{ 1, 0, -r };

	/* If dx is in Im(d_{r-2}) then x is in Ker(d_r) */
	const Staircase& sc = GetRecentStaircase(basis_ss, deg_dx);
	size_t first_r = GetFirstIndexOfLevel(sc, r);
	dx = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_r, dx);
	if (dx.empty()) {
		if (x != array{ -1 } && !x.empty())
			AddDiff(basis_ss, deg_x, std::move(x), { -1 }, r + 2);
		return;
	}

	array image_new;
	int level_image_new = -1;
	if (x == array{ -1 }) { /* If the source is unknown, check if it can be hit and then insert it to the end of level r. */
		size_t first_rp2 = GetFirstIndexOfLevel(sc, r + 2);
		dx = lina::Residue(sc.basis_ind.begin() + first_r, sc.basis_ind.begin() + first_rp2, dx);
		if (!dx.empty())
			basis_ss.back()[deg_dx] = triangularize(sc, first_rp2, dx, x, r, image_new, level_image_new);
	}
	else { /* Otherwise insert it to the beginning of level r */
		basis_ss.back()[deg_dx] = triangularize(sc, first_r, dx, x, r, image_new, level_image_new);
	}

	if (level_image_new != -1) {
		if (level_image_new < kLevelMax / 2) { /* Add a d_{r1-2} image */
			Deg deg_image_new = deg_dx + Deg{ 1, 0, -level_image_new };
			AddImage(basis_ss, deg_image_new, std::move(image_new), { -1 }, level_image_new - 2);
		}
		else { /* Add a d_r1 cycle */
			int r_image = kLevelMax - level_image_new;
			Deg deg_image_new = deg_dx - Deg{ 1, 0, -r_image };
			AddDiff(basis_ss, deg_image_new, std::move(image_new), { -1 }, r_image + 2);
		}
	}
}

/* Add d_r(x)=dx and all its implications.
 * deg_x must be in basis_ss. */
void SetDiff(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, Deg deg_x, array x, array dx, int r, int t_max /*= -1*/)
{
	if (t_max == -1)
		t_max = basis.rbegin()->first.t;

	Deg deg_dx = deg_x + Deg{ 1, 0, -r };
	for (auto& [deg_y, basis_ss_d_original] : basis_ss.front()) {
		const Staircase& basis_ss_d = GetRecentStaircase(basis_ss, deg_y);
		Deg deg_dy = deg_y + Deg{ 1, 0, -r };
		if (deg_x.t + deg_y.t > t_max)
			break;
		for (size_t i = 0; i < basis_ss_d.levels.size(); ++i) {
			if (basis_ss_d.levels[i] == kLevelMax - r) { /* y is not a d_r cycle */
				if (basis_ss_d.diffs_ind[i] != array{ -1 }) {
					Poly poly_x = indices_to_Poly(x, basis.at(deg_x));
					Poly poly_y = indices_to_Poly(basis_ss_d.basis_ind[i], basis.at(deg_y));
					Poly poly_xy = grbn::Reduce(poly_x * poly_y, gb);
					array xy = poly_xy.empty() ? array{} : Poly_to_indices(poly_xy, basis.at(deg_x + deg_y));
					Poly poly_dx = dx.empty() ? Poly{} : indices_to_Poly(dx, basis.at(deg_dx));
					Poly poly_dy = basis_ss_d.diffs_ind[i].empty() ? Poly{} : indices_to_Poly(basis_ss_d.diffs_ind[i], basis.at(deg_dy));
					Poly poly_dxy = grbn::Reduce(poly_x * poly_dy + poly_dx * poly_y, gb);
					array dxy = poly_dxy.empty() ? array{} : Poly_to_indices(poly_dxy, basis.at(deg_x + deg_dy));

					AddDiff(basis_ss, deg_x + deg_y, std::move(xy), std::move(dxy), r);
				}
			}
			else if (r <= basis_ss_d.levels[i] && basis_ss_d.levels[i] < kLevelMax - r) { /* y is a d_r cycle */
				Poly poly_x = indices_to_Poly(x, basis.at(deg_x));
				Poly poly_y = indices_to_Poly(basis_ss_d.basis_ind[i], basis.at(deg_y));
				Poly poly_xy = grbn::Reduce(poly_x * poly_y, gb);
				array xy = poly_xy.empty() ? array{} : Poly_to_indices(poly_xy, basis.at(deg_x + deg_y));
				Poly poly_dx = dx.empty() ? Poly{} : indices_to_Poly(dx, basis.at(deg_dx));
				Poly poly_dxy = grbn::Reduce(poly_dx * poly_y, gb);
				array dxy = poly_dxy.empty() ? array{} : Poly_to_indices(poly_dxy, basis.at(deg_x + deg_dy));

				AddDiff(basis_ss, deg_x + deg_y, std::move(xy), std::move(dxy), r);
			}
		}
	}
}
