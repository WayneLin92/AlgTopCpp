#include "main.h"
#include "groebner.h"
#include "benchmark.h"
#include "linalg.h"
#include "myexception.h"
#include "myio.h"
#include <iostream>
#include <sstream>

constexpr auto kLevelMin = 4;

/* A custom Exception class */
class SSException : public MyException
{
public:
	SSException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
	SSException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

class NoException : public MyException
{
public:
	NoException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
	NoException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

std::ostream& operator<<(std::ostream& sout, const Staircase& sc)
{
	sout << "Staircase:\n";
	for (size_t i = 0; i < sc.levels.size(); ++i) {
		sout << sc.levels[i] << "  " << sc.basis_ind[i] << "  " << sc.diffs_ind[i] << '\n';
	}
	return sout;
}

/* Return the newest version of the staircase in history */
Staircase& GetStaircase(std::vector<std::map<Deg, Staircase>>& basis_ss, const Deg& deg)
{
	for (auto p = basis_ss.rbegin(); p != basis_ss.rend(); ++p)
		if (p->find(deg) != p->end())
			return p->at(deg);
	throw MyException(0x553989e0U, "BUG: deg not found");
}

inline size_t GetFirstIndex(const Staircase& sc, int level)
{
	return std::lower_bound(sc.levels.begin(), sc.levels.end(), level) - sc.levels.begin();
}

/* Count undetermined differentials
 * Return (count, index) */
std::pair<int, int> CountUndertermined(const Staircase& sc, int level_min)
{
	int count_undermined = 0;
	int index = -1;
	for (int i = 0; i < (int)sc.diffs_ind.size(); ++i) {
		if (sc.diffs_ind[i] == array{ -1 } && sc.levels[i] >= level_min) {
			++count_undermined;
			index = i;
		}
	}
	return std::make_pair(count_undermined, index);
}

/* Return if an element in deg and level can be hit */
bool CanBeHit(std::vector<std::map<Deg, Staircase>>& basis_ss, const Deg& deg, int level)
{
	int r_max = std::min(level, 2 * (deg.t - deg.s + 1) - deg.v);
	for (int r = kLevelMin; r <= r_max; r += 2) {
		Deg d_src = deg - Deg{ 1, 0, -r };
		if (basis_ss.front().find(d_src) != basis_ss.front().end()) {
			auto pair = CountUndertermined(GetStaircase(basis_ss, d_src), kLevelMax - r); // Can be optimized
			if (pair.first > 0)
				return true;
		}
	}
	return false;
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
	sc_new.levels   .insert(sc_new.levels   .end(), sc.levels   .begin(), sc.levels   .begin() + i_insert);

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
		throw MyException(0xdd696043U, "BUG");
#endif
	return sc_new;
}

void AddImage(std::vector<std::map<Deg, Staircase>>& basis_ss, const Deg& deg_dx, array dx, array x, int r, bool bForSSToF2 = false);

/* Add d_r(x)=dx and d_r^{-1}(dx)=x. */
void AddDiff(std::vector<std::map<Deg, Staircase>>& basis_ss, const Deg& deg_x, array x, array dx, int r, bool bForSSToF2 = false)
{
	Deg deg_dx = deg_x + Deg{ 1, 0, -r };

	/* If x is zero then dx is in Im(d_{r-2}) */
	if (x.empty()) {
		if (dx != array{ -1 } && !dx.empty())
			AddImage(basis_ss, deg_dx, dx, { -1 }, r - 2, bForSSToF2);
		return;
	}

	/* If x is in Ker(d_r) then dx is in Im(d_{r-2}) */
	const Staircase& sc = GetStaircase(basis_ss, deg_x);
	size_t first_Nmr = GetFirstIndex(sc, kLevelMax - r);
	x = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_Nmr, x);
	if (x.empty()) {
		if (dx != array{ -1 } && !dx.empty())
			AddImage(basis_ss, deg_dx, dx, { -1 }, r - 2, bForSSToF2);
		return;
	}

	array image_new;
	int level_image_new = -1;
	if (dx == array{ -1 }) { /* If the target is unknown, insert it to the end of level N-r. */
		size_t first_Nmrp2 = GetFirstIndex(sc, kLevelMax - r + 2);
		x = lina::Residue(sc.basis_ind.begin() + first_Nmr, sc.basis_ind.begin() + first_Nmrp2, x);
		if (!x.empty()) {
			basis_ss.back()[deg_x] = triangularize(sc, first_Nmrp2, x, { -1 }, kLevelMax - r, image_new, level_image_new);
		}
	}
	else if (dx.empty()) { /* If the target is zero, insert it to the end of level N-r-2 */
		basis_ss.back()[deg_x] = triangularize(sc, first_Nmr, x, { -1 }, kLevelMax - r - 2, image_new, level_image_new);
	}
	else { /* Otherwise insert it to the beginning of level N-r */
		basis_ss.back()[deg_x] = triangularize(sc, first_Nmr, x, dx, kLevelMax - r, image_new, level_image_new);
	}
	
	if (level_image_new != -1) {
		if (level_image_new < kLevelMax / 2) { /* Add a d_{r1-2} image */
			Deg deg_image_new = deg_x + Deg{ 1, 0, -level_image_new };
			AddImage(basis_ss, deg_image_new, std::move(image_new), { -1 }, level_image_new - 2, bForSSToF2);
		}
		else { /* Add a d_r1 cycle */
			int r_image = kLevelMax - level_image_new;
			Deg deg_image_new = deg_x - Deg{ 1, 0, -r_image };
			AddDiff(basis_ss, deg_image_new, std::move(image_new), { -1 }, r_image + 2, bForSSToF2);
		}
	}

	/* Add image */
	if (dx != array{ -1 } && !dx.empty())
		AddImage(basis_ss, deg_dx, std::move(dx), std::move(x), r, bForSSToF2);
}

/* Add an image of d_r.
 * Assume dx is nonempty. */
void AddImage(std::vector<std::map<Deg, Staircase>>& basis_ss, const Deg& deg_dx, array dx, array x, int r, bool bForSSToF2)
{
	Deg deg_x = deg_dx - Deg{ 1, 0, -r };

	/* If dx is in Im(d_{r-2}) then x is in Ker(d_r) */
	const Staircase& sc = GetStaircase(basis_ss, deg_dx);
	size_t first_r = GetFirstIndex(sc, r);
	dx = lina::Residue(sc.basis_ind.begin(), sc.basis_ind.begin() + first_r, dx);
	if (dx.empty()) {
		if (x != array{ -1 } && !x.empty())
			AddDiff(basis_ss, deg_x, std::move(x), { -1 }, r + 2, bForSSToF2);
		return;
	}

	array image_new;
	int level_image_new = -1;
	if (x == array{ -1 }) { /* If the source is unknown, check if it can be hit and then insert it to the end of level r. */
		size_t first_rp2 = GetFirstIndex(sc, r + 2);
		dx = lina::Residue(sc.basis_ind.begin() + first_r, sc.basis_ind.begin() + first_rp2, dx);
		if (!dx.empty()) {
			if (bForSSToF2 && !CanBeHit(basis_ss, deg_dx, r)) {
				std::stringstream ss;
				ss << "impossible d_{r-2} image at " << "deg=" << deg_dx << " level=" << r << '\n';
				throw SSException(0xd485ef5eU, ss.str());
			}
			basis_ss.back()[deg_dx] = triangularize(sc, first_rp2, dx, x, r, image_new, level_image_new);
		}
	}
	else { /* Otherwise insert it to the beginning of level r */
		basis_ss.back()[deg_dx] = triangularize(sc, first_r, dx, x, r, image_new, level_image_new);
	}

	if (level_image_new != -1) {
		if (level_image_new < kLevelMax / 2) { /* Add a d_{r1-2} image */
			Deg deg_image_new = deg_dx + Deg{ 1, 0, -level_image_new };
			AddImage(basis_ss, deg_image_new, std::move(image_new), { -1 }, level_image_new - 2, bForSSToF2);
		}
		else { /* Add a d_r1 cycle */
			int r_image = kLevelMax - level_image_new;
			Deg deg_image_new = deg_dx - Deg{ 1, 0, -r_image };
			AddDiff(basis_ss, deg_image_new, std::move(image_new), { -1 }, r_image + 2, bForSSToF2);
		}
	}
}

/* Add d_r(x)=dx and all its implications.
 * deg_x must be in basis_ss. */
void SetDiff(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, std::vector<std::map<Deg, Staircase>>& basis_ss, Deg deg_x, array x, array dx, int r, bool bForSSToF2 = false)
{
	Deg deg_dx = deg_x + Deg{ 1, 0, -r };
	int t_max = basis.rbegin()->first.t;

	for (auto& [deg_y, basis_ss_d_original] : basis_ss.front()) {
		Staircase& basis_ss_d = GetStaircase(basis_ss, deg_y);
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

					AddDiff(basis_ss, deg_x + deg_y, xy, dxy, r, bForSSToF2);
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

				AddDiff(basis_ss, deg_x + deg_y, xy, dxy, r, bForSSToF2);
			}
		}
	}
}

/* Check if the element is the onlyone in the level */
bool isSingle(const Staircase& sc, int index)
{
	return (index == 0 || sc.levels[(size_t)index - 1] != sc.levels[index]) && (index == (int)sc.levels.size() - 1 || sc.levels[index] != sc.levels[(size_t)index + 1]);
}

int DeduceDiffsForSSToF2(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, std::vector<std::map<Deg, Staircase>>& basis_ss)
{
	int count_new_diffs = 0;
	/* Look for pairs of basis that can be only killed in one way */
	for (auto p = basis_ss.front().begin(); p != basis_ss.front().end(); ++p)
	{
		const Deg deg = p->first;
		const Staircase& staircase = GetStaircase(basis_ss, deg);

		/* Check if there is only one base with undermined differential */
		auto pair = CountUndertermined(staircase, 0);
		int count = pair.first, index = pair.second;
		int count_tgt = 0, index_tgt = -1;
		int count_src = 0, index_src = -1;
		Deg deg_tgt{ -1, -1, -1 }, deg_src{ -1, -1, -1 };
		if (count == 1) { /* If there is only one base, check if there is only one way to kill it */
			int level = staircase.levels[index];
			if (level > kLevelMax / 2) {
				/* Count how many chances it can be a source of a differential */
				for (int r = kLevelMax - level; r <= deg.v; r += 2) {
					Deg d_tgt = deg + Deg{ 1, 0, -r };
					if (basis_ss.front().find(d_tgt) != basis_ss.front().end()) {
						auto pair_tgt = CountUndertermined(GetStaircase(basis_ss, d_tgt), 0);
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
					auto pair_src = CountUndertermined(GetStaircase(basis_ss, d_src), kLevelMax - r);
					count_src += pair_src.first;
					if (pair_src.first > 0) {
						index_src = pair_src.second;
						deg_src = d_src;
					}
				}
			}

			if (count_tgt == 0 && count_src == 1) { /* It is a target */
				std::cout << "type1: (" << deg_src.t - deg_src.s << ", " << deg_src.v / 2 << ") --> (" << deg.t - deg.s << ", " << deg.v / 2 << ")\n";
				if (isSingle(staircase, index)) {
					++count_new_diffs;
					SetDiff(gb, basis, basis_ss, deg_src, GetStaircase(basis_ss, deg_src).basis_ind[index_src], staircase.basis_ind[index], deg_src.v - deg.v, true);
				}
			}
			else if (count_tgt == 1 && count_src == 0) { /* It is a source */
				std::cout << "type2: (" << deg.t - deg.s << ", " << deg.v / 2 << ") --> (" << deg_tgt.t - deg_tgt.s << ", " << deg_tgt.v / 2 << ")\n";
				Staircase sc_tgt = GetStaircase(basis_ss, deg_tgt);
				if (isSingle(sc_tgt, index_tgt)) {
					++count_new_diffs;
					SetDiff(gb, basis, basis_ss, deg, staircase.basis_ind[index], sc_tgt.basis_ind[index_tgt], deg.v - deg_tgt.v, true);
				}
			}
			else if ((level < kLevelMax / 2 || count_tgt == 0) && count_src == 0)
				throw SSException(0x9eacf10fU, "Nontrivial permanant cycle");
		}
	}
	return count_new_diffs;
}

/* generate the table of the spectral sequence */
void deduce_diffs(const Database& db, const std::string& table_prefix, int t_max)
{
	const grbn::GbWithCache gb = db.load_gb(table_prefix + "_relations", t_max);
	const std::map<Deg, Mon1d> basis = db.load_basis(table_prefix + "_basis", t_max);
	std::vector<std::map<Deg, Staircase>> basis_ss = { db.load_basis_ss(table_prefix + "_ss", t_max) };

	try {
		while (DeduceDiffsForSSToF2(gb, basis, basis_ss)) ;
	}
#ifdef _DEBUG
	catch (SSException& e)
	{
		std::cerr << "Error-" << std::hex << e.id() << ": " << e.what() << '\n';
	}
#else
	catch (NoException& e) {} /* This does not catch any exception */
#endif
	
	/* insert into the database */
	db.begin_transaction();
	db.execute_cmd("DELETE FROM " + table_prefix + "_ss_tmp;");
	db.save_basis_ss(table_prefix + "_ss_tmp", basis_ss.front());
	db.end_transaction();
}
