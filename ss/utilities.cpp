#include "main.h"
#include "benchmark.h"
#include "linalg.h"
#include "myio.h"
#include <iostream>

DiffType GetDiffType(const std::string& table_prefix)
{
	if (table_prefix == "E4")
		return DiffType::E4;
	else if (table_prefix == "E4t" || table_prefix == "E4t1")
		return DiffType::E4t;
	else if (table_prefix == "E4h0")
		return DiffType::E4h0;
	else
		throw MyException(0x10599a49U, "cannot determine the diff type for the table " + table_prefix);
};

bool IsEt(const std::string& table_prefix) //
{
	if (table_prefix == "E4")
		return false;
	else if (table_prefix == "E4t" || table_prefix == "E4h0" || table_prefix == "E4t1")
		return true;
	else
		throw MyException(0x3410c3c2U, "cannot determine the converging type for the table " + table_prefix);
}

bool IsEt(DiffType d_type)
{
	switch (d_type)
	{
	case DiffType::E4:
		return false;
	case DiffType::E4t:
	case DiffType::E4h0:
		return true;
	default:
		throw MyException(0x3567d341U, "unhandle difftype");
	}
}

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

/* Find the maximum level of a null (inv)differential
** Return -1 if not found */
int GetMaxLevelOfNull(const Staircase& sc)
{
	size_t i = sc.levels.size();
	while(i-- > 0) {
		if (sc.diffs_ind[i] == array{ -1 })
			return sc.levels[i];
	}
	return -1;
}

size_t GetFirstIndexOnLevel(const Staircase& sc, int level)
{
	return std::lower_bound(sc.levels.begin(), sc.levels.end(), level) - sc.levels.begin();
}

/* Find the position of the first base with null diff in level */
size_t GetFirstIndexOfNullOnLevel(const Staircase& sc, int level)
{
	array::const_iterator first = std::lower_bound(sc.levels.begin(), sc.levels.end(), level);
	array::const_iterator last = std::lower_bound(first, sc.levels.end(), level + 2);

	array::const_iterator it;
	ptrdiff_t count, step;
	count = std::distance(first, last);
	while (count > 0)
	{
		it = first; step = count / 2; advance(it, step);
		if (*(sc.diffs_ind.begin() + (it - sc.levels.begin())) != array{ -1 }) {
			first = ++it;
			count -= step + 1;
		}
		else count = step;
	}
	return first - sc.levels.begin();
}

/* Return if a possible source of an image can be found */
bool ExistsSrc(const Staircases1d& basis_ss, const Deg& deg, DiffType d_type, int r)
{
	auto [s_diff, t_diff] = GetST(d_type);
	int r_max = std::min(r, 2 * (deg.t - deg.s + 1) - deg.v);
	for (int r1 = kLevelMin; r1 <= r_max; r1 += 2) {
		Deg d_src = deg - Deg{ s_diff, t_diff, -r1 };
		if (basis_ss.front().find(d_src) != basis_ss.front().end()) {
			if (GetMaxLevelOfNull(GetRecentStaircase(basis_ss, d_src)) >= kLevelMax - r1)
				return true;
		}
	}
	return false;
}

/* Return if a null differential can be resolved by a null target */
bool ExistsTgt(const Staircases1d& basis_ss, const Deg& deg, DiffType d_type, int r)
{
	auto [s_diff, t_diff] = GetST(d_type);
	for (int r1 = r; r1 <= deg.v; r1 += 2) {
		Deg d_tgt = deg + Deg{ s_diff, t_diff, -r1 };
		if (basis_ss.front().find(d_tgt) != basis_ss.front().end()) {
			if (GetMaxLevelOfNull(GetRecentStaircase(basis_ss, d_tgt)) >= r1)
				return true;
		}
	}
	return false;
}

/* Return the first index of a source level such that all levels above are already fixed */
size_t GetFirstIndexOfFixedLevels(const Staircases1d& basis_ss, const Deg& deg, DiffType d_type, int level)
{
	auto [s_diff, t_diff] = GetST(d_type);
	const Staircase& sc = GetRecentStaircase(basis_ss, deg);
	size_t result = sc.levels.size();
	for (size_t i = sc.levels.size(); i-- > 0; ) {
		if (sc.diffs_ind[i] == array{ -1 } || sc.levels[i] < level)
			break;
		if (i == 0 || sc.levels[i - 1] != sc.levels[i]) {
			int r = kLevelMax - sc.levels[i];
			if (ExistsSrc(basis_ss, deg + Deg{ s_diff, t_diff, -r }, d_type, r - 2))
				break;
			else
				result = i;
		}
	}
	return result;
}

/* Return if there is no null differential in degree t */
bool NoNullDiff(const Staircases1d& basis_ss, int t) //E4h0
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
std::pair<int, int> CountDrTgt(const Staircases1d& basis_ss, const Deg& deg_tgt, DiffType d_type, int r) //
{
	std::pair<int, int> result;
	if (basis_ss.front().find(deg_tgt) != basis_ss.front().end()) {
		const Staircase& sc_tgt = GetRecentStaircase(basis_ss, deg_tgt);
		result.second = (int)GetFirstIndexOnLevel(sc_tgt, r);
		result.first = (int)GetFirstIndexOfFixedLevels(basis_ss, deg_tgt, d_type, kLevelMax - r + 2) - result.second;
	}
	else
		result = { 0, -1 };
	return result;
}

/* Count the number of all possible d_{r1} targets for r1>=r.
 * Return (count, deg, index).
 * count>=2 means there are >= 2 elements in the range and 
 * deg corresponds to the one with the least r1. */
std::tuple<int, Deg, int> CountTgt(const Staircases1d& basis_ss, const Deg& deg, DiffType d_type, int r)
{
	auto [s_diff, t_diff] = GetST(d_type);
	Deg deg_tgt; int count = 0, index = -1;
	const Staircase& sc = GetRecentStaircase(basis_ss, deg);
	for (int r1 = r; r1 <= deg.v; r1 += 2) {
		Deg d_tgt = deg + Deg{ s_diff, t_diff, -r1 };
		auto [c, first_r1] = CountDrTgt(basis_ss, d_tgt, d_type, r1);
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
std::tuple<int, Deg, int> CountSrc(const Staircases1d& basis_ss, const Deg& deg, DiffType d_type, int level)
{
	auto [s_diff, t_diff] = GetST(d_type);
	Deg deg_src; int count = 0, index = -1;
	int r_max = std::min(level, 2 * (deg.t - deg.s + 1) - deg.v); //
	const Staircase& sc = GetRecentStaircase(basis_ss, deg);
	for (int r1 = kLevelMin; r1 <= r_max; r1 += 2) {
		Deg d_src = deg - Deg{ s_diff, t_diff, -r1 };
		if (basis_ss.front().find(d_src) != basis_ss.front().end()) {
			const Staircase& sc_src = GetRecentStaircase(basis_ss, d_src);
			int first_Nmr1 = (int)GetFirstIndexOnLevel(sc_src, kLevelMax - r1);
			int c = (int)GetFirstIndexOfFixedLevels(basis_ss, d_src, d_type, kLevelMax - r1 + 2) - first_Nmr1;
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

NullDiff GetNullDiff(const Staircases1d& basis_ss, const Deg& deg, DiffType d_type)
{
	auto [s_diff, t_diff] = GetST(d_type);
	NullDiff result{ -1, -1, -1 };
	const Staircase& sc = GetRecentStaircase(basis_ss, deg);
	for (size_t i = sc.diffs_ind.size(); i-- > 0;) {
		if (sc.diffs_ind[i] == array{ -1 } && sc.levels[i] > kLevelPC) {
			int r = kLevelMax - sc.levels[i];
			Deg deg_tgt = deg + Deg{ s_diff, t_diff, -r };
			auto [num_tgts, first_r] = CountDrTgt(basis_ss, deg_tgt, d_type, r);
			result = { (int)i, num_tgts, first_r };
			break;
		}
	}
	return result;
}

void NullDiffs::InitNullDiffs(const Staircases1d& basis_ss, DiffType d_type, int t_max, bool bNew)
{
	const std::map<Deg, Staircase>& basis_ss0 = bNew ? basis_ss.front() : basis_ss.back();
	if (t_max == -1)
		t_max = INT_MAX;
	for (auto p = basis_ss0.begin(); p != basis_ss0.end(); ++p) {
		if (p->first.t > t_max)
			break;
		null_diffs_[p->first] = GetNullDiff(basis_ss, p->first, d_type);
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
void CacheDeduction::InitNullDiffs(const Staircases1d& basis_ss, DiffType d_type, int t_max)
{
	bool bNew = null_diffs_.size() == 1;
	const std::map<Deg, Staircase>& basis_ss0 = bNew ? basis_ss.front() : basis_ss.back();
	if (t_max == -1)
		t_max = INT_MAX;

	/* Initialize null_diffs_.back() */
	for (auto p = basis_ss0.begin(); p != basis_ss0.end(); ++p) {
		if (p->first.t > t_max)
			break;
		null_diffs_.back()[p->first] = GetNullDiff(basis_ss, p->first, d_type);
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
void triangularize(Staircase& sc, size_t i_insert, array x, array dx, int level, array& image, int& level_image)
{
	level_image = -1;

	size_t i = i_insert;
	while (!x.empty()) {
		std::swap(x, sc.basis_ind[i]);
		std::swap(dx, sc.diffs_ind[i]);
		std::swap(level, sc.levels[i]);

		++i;
		for (size_t j = i_insert; j < i; ++j) {
			if (std::binary_search(x.begin(), x.end(), sc.basis_ind[j][0])) {
				x = lina::AddVectors(x, sc.basis_ind[j]);
				if (level == sc.levels[j] && dx != array{ -1 })
					dx = lina::AddVectors(dx, sc.diffs_ind[j]);
			}
		}
	}
	if (dx != array{ -1 } && !dx.empty()) {
		image = std::move(dx);
		level_image = kLevelMax - level;
	}

	/* Triangularize the rest */
	for (; i < sc.basis_ind.size(); ++i) {
		for (size_t j = i_insert; j < i; ++j) {
			if (std::binary_search(sc.basis_ind[i].begin(), sc.basis_ind[i].end(), sc.basis_ind[j][0])) {
				sc.basis_ind[i] = lina::AddVectors(sc.basis_ind[i], sc.basis_ind[j]);
				if (sc.levels[i] == sc.levels[j] && sc.diffs_ind[i] != array{ -1 })
					sc.diffs_ind[i] = lina::AddVectors(sc.diffs_ind[i], sc.diffs_ind[j]);
			}
		}
#ifdef _DEBUG
		if (sc.basis_ind[i].empty())
			throw MyException(0xfe35902dU, "BUG: triangularize()");
#endif
	}
}

void UpdateStaircase(Staircases1d& basis_ss, const Deg& deg, const Staircase& sc_i, size_t i_insert, array x, array dx, int level, array& image, int& level_image)
{
	if (basis_ss.back().find(deg) == basis_ss.back().end())
		basis_ss.back()[deg] = sc_i;
	triangularize(basis_ss.back()[deg], i_insert, std::move(x), std::move(dx), level, image, level_image);
}

/* Return d_r(x) */
array GetSsDiff(Staircases1d& basis_ss, const Deg& deg_x, array x, int r)
{
	if (x.empty())
		return array{};
	const Staircase& sc = GetRecentStaircase(basis_ss, deg_x);
	size_t first_Nmr = GetFirstIndexOnLevel(sc, kLevelMax - r);
	size_t first_Nmrp2 = GetFirstIndexOnLevel(sc, kLevelMax - r + 2);
	x = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_Nmr, x); /* Compute x mod Ker(d_r) */
	if (lina::Residue(sc.basis_ind.begin() + first_Nmr, sc.basis_ind.begin() + first_Nmrp2, x).empty()) {
		return lina::GetImage(sc.basis_ind.begin() + first_Nmr, sc.basis_ind.begin() + first_Nmrp2,
			sc.diffs_ind.begin() + first_Nmr, sc.diffs_ind.begin() + first_Nmrp2, x);
	}
	else
		return array{ -1 };
}

bool IsNewDiff(Staircases1d& basis_ss, const Deg& deg_x, DiffType d_type, array x, array dx, int r)
{
	auto [s_diff, t_diff] = GetST(d_type);
	Deg deg_dx = deg_x + Deg{ s_diff, t_diff, -r };
	array dx1 = GetSsDiff(basis_ss, deg_x, x, r);
	if (dx1 == array{ -1 })
		return true;
	if (basis_ss.front().find(deg_dx) != basis_ss.front().end()) {
		const Staircase& sc = GetRecentStaircase(basis_ss, deg_dx);
		size_t first_r = GetFirstIndexOnLevel(sc, r);
		return !lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_r, lina::AddVectors(dx, dx1)).empty();
	}
	else
		return !dx1.empty();
}

/* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
void AddDiff(Staircases1d& basis_ss, const Deg& deg_x, DiffType d_type, array x, array dx, int r)
{
	auto [s_diff, t_diff] = GetST(d_type);
	Deg deg_dx = deg_x + Deg{ s_diff, t_diff, -r };
	int t_max = basis_ss.front().rbegin()->first.t;

	/* If x is zero then dx is in Im(d_{r-2}) */
	if (x.empty()) {
		if (dx != array{ -1 } && !dx.empty())
			AddImage(basis_ss, deg_dx, d_type, std::move(dx), { -1 }, r - 2);
		return;
	}

	const Staircase& sc = GetRecentStaircase(basis_ss, deg_x);
	size_t first_Nmr = GetFirstIndexOnLevel(sc, kLevelMax - r);
	x = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_Nmr, x);
	if (x.empty()) { /* If x is in Ker(d_r) then dx is in Im(d_{r-2}) */
		if (dx != array{ -1 } && !dx.empty())
			AddImage(basis_ss, deg_dx, d_type, std::move(dx), { -1 }, r - 2);
		return;
	}

	array image_new;
	int level_image_new = -1;
	/* If the target is unknown, insert it to the end of level N-r.
	** Check if null diff can be resolved when bForEt=True */
	if (dx == array{ -1 }) {
		size_t first_Nmrp2 = GetFirstIndexOnLevel(sc, kLevelMax - r + 2);
		x = lina::Residue(sc.basis_ind.begin() + first_Nmr, sc.basis_ind.begin() + first_Nmrp2, x);
		if (!x.empty()) {
			if (IsEt(d_type) && (deg_x.t - t_diff) <= t_max && !ExistsSrc(basis_ss, deg_x, d_type, kLevelMax) && !ExistsTgt(basis_ss, deg_x, d_type, r))
				throw SSException(0x23e587dU, "No target to resolve the null differential");
			UpdateStaircase(basis_ss, deg_x, sc, first_Nmrp2, x, { -1 }, kLevelMax - r, image_new, level_image_new);
		}
	}
	else if (dx.empty()) { /* If the target is zero, insert it to the end of level N-r-2 */
		if (IsEt(d_type) && (deg_x.t - t_diff) <= t_max && !ExistsSrc(basis_ss, deg_x, d_type, kLevelMax) && !ExistsTgt(basis_ss, deg_x, d_type, r + 2))
			throw SSException(0x32235858U, "No target to resolve the null differential");
		UpdateStaircase(basis_ss, deg_x, sc, first_Nmr, x, { -1 }, kLevelMax - r - 2, image_new, level_image_new);
	}
	else { /* Otherwise insert it to the beginning of level N-r */
		//size_t first_Nmrp2_null = GetFirstIndexOfNullOnLevel(sc, kLevelMax - r + 2);
		UpdateStaircase(basis_ss, deg_x, sc, first_Nmr, x, dx, kLevelMax - r, image_new, level_image_new); //
	}

	if (level_image_new != -1) {
		if (level_image_new < kLevelMax / 2) { /* Add a d_{r1-2} image */
			Deg deg_image_new = deg_x + Deg{ s_diff, t_diff, -level_image_new };
			AddImage(basis_ss, deg_image_new, d_type, std::move(image_new), { -1 }, level_image_new - 2);
		}
		else { /* Add a d_r1 cycle */
			int r_image = kLevelMax - level_image_new;
			Deg deg_image_new = deg_x - Deg{ s_diff, t_diff, -r_image };
			AddDiff(basis_ss, deg_image_new, d_type, std::move(image_new), { -1 }, r_image + 2);
		}
	}

	/* Add image */
	if (dx != array{ -1 } && !dx.empty())
		AddImage(basis_ss, deg_dx, d_type, std::move(dx), std::move(x), r);
}

/* Add an image of d_r.
 * Assume dx is nonempty. */
void AddImage(Staircases1d& basis_ss, const Deg& deg_dx, DiffType d_type, array dx, array x, int r) //
{
	auto [s_diff, t_diff] = GetST(d_type);
	Deg deg_x = deg_dx - Deg{ s_diff, t_diff, -r };

	/* If dx is in Im(d_{r-2}) then x is in Ker(d_r) */
	const Staircase& sc = GetRecentStaircase(basis_ss, deg_dx);
	size_t first_r = GetFirstIndexOnLevel(sc, r);
	dx = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_r, dx);
	if (dx.empty()) {
		if (x != array{ -1 } && !x.empty())
			AddDiff(basis_ss, deg_x, d_type, std::move(x), { -1 }, r + 2);
		return;
	}

	array image_new;
	int level_image_new = -1;
	if (x == array{ -1 }) { /* If the source is unknown, check if it can be hit and then insert it to the end of level r. */
		size_t first_rp2 = GetFirstIndexOnLevel(sc, r + 2);
		dx = lina::Residue(sc.basis_ind.begin() + first_r, sc.basis_ind.begin() + first_rp2, dx);
		if (!dx.empty()) {
			if (!ExistsSrc(basis_ss, deg_dx, d_type, r))
				throw SSException(0x75989376U, "No source for the image");
			UpdateStaircase(basis_ss, deg_dx, sc, first_rp2, dx, x, r, image_new, level_image_new);
		}
	}
	else { /* Otherwise insert it to the beginning of level r */
		UpdateStaircase(basis_ss, deg_dx, sc, first_r, dx, x, r, image_new, level_image_new); //
	}

	if (level_image_new != -1) {
		if (level_image_new < kLevelMax / 2) { /* Add a d_{r1-2} image */
			Deg deg_image_new = deg_dx + Deg{ s_diff, t_diff, -level_image_new };
			AddImage(basis_ss, deg_image_new, d_type, std::move(image_new), { -1 }, level_image_new - 2);
		}
		else { /* Add a d_r1 cycle */
			int r_image = kLevelMax - level_image_new;
			Deg deg_image_new = deg_dx - Deg{ s_diff, t_diff, -r_image };
			AddDiff(basis_ss, deg_image_new, d_type, std::move(image_new), { -1 }, r_image + 2);
		}
	}
}

/* Add d_r(x)=dx and all its implications.
 * deg_x must be in basis_ss. */
void SetDiff(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, Deg deg_x, DiffType d_type, array x, array dx, int r, int t_max)
{
	auto [s_diff, t_diff] = GetST(d_type);
	if (t_max == -1)
		t_max = basis.rbegin()->first.t;

	Deg deg_dx = deg_x + Deg{ s_diff, t_diff, -r };
	for (auto& [deg_y, basis_ss_d_original] : basis_ss.front()) {
		const Staircase& basis_ss_d = GetRecentStaircase(basis_ss, deg_y);
		Deg deg_dy = deg_y + Deg{ s_diff, t_diff, -r };
		if (deg_x.t + deg_y.t > t_max || deg_x.t + deg_y.t + t_diff > t_max)
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

					AddDiff(basis_ss, deg_x + deg_y, d_type, std::move(xy), std::move(dxy), r);
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

				AddDiff(basis_ss, deg_x + deg_y, d_type, std::move(xy), std::move(dxy), r);
			}
		}
	}
}

/* Check first if it is a new differential before adding differentials */
int SetDiffV2(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, const Deg& deg_x, const Deg& deg_dx, DiffType d_type, array x, array dx, int r, int t_max)
{
	if (t_max == -1)
		t_max = basis.rbegin()->first.t;
	if (deg_x.t > t_max) //TODO: Consider t_diff
		return 0;

	if (x.empty()) {
		if (!dx.empty()) {
			SetDiff(gb, basis, basis_ss, deg_dx, d_type, dx, {}, kLevelMax / 4, t_max); //
			AddImage(basis_ss, deg_dx, d_type, std::move(dx), array{ -1 }, r - 2); //
			return 1;
		}
	}
	else if (IsNewDiff(basis_ss, deg_x, d_type, x, dx, r)) {
		SetDiff(gb, basis, basis_ss, deg_x, d_type, std::move(x), std::move(dx), r, t_max);
		return 1;
	}
	return 0;
}
