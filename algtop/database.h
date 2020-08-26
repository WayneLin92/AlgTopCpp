#ifndef DATABSE_H
#define DATABSE_H

#include "mymath.h"
#include "sqlite3/sqlite3.h"
#include <iostream>
#include <vector>
#include <array>
#include <map>

constexpr auto T_MAX = 10000;
inline const char* sqlite3_column_str(sqlite3_stmt* stmt, int iCol) { return reinterpret_cast<const char*>(sqlite3_column_text(stmt, iCol)); }
inline int sqlite3_bind_str(sqlite3_stmt* stmt, int iCol, const std::string& str) { return sqlite3_bind_text(stmt, iCol, str.c_str(), -1, SQLITE_TRANSIENT); }
inline int sqlite3_prepare_v100(sqlite3* db, const char* zSql, sqlite3_stmt** ppStmt) { return sqlite3_prepare_v2(db, zSql, int(strlen(zSql)) + 1, ppStmt, NULL); }
inline int sqlite3_prepare_v100(sqlite3* db, std::string sql, sqlite3_stmt** ppStmt) { return sqlite3_prepare_v2(db, sql.c_str(), int(sql.size()) + 1, ppStmt, NULL); }

/*--------- database.cpp ---------*/
std::string array_to_str(array::const_iterator pbegin, array::const_iterator pend);
array str_to_array(const char* str_mon);
std::string array2d_to_str(array2d::const_iterator pbegin, array2d::const_iterator pend);
array2d str_to_array2d(const char* str_poly);

inline std::string array_to_str(const array& a) { return array_to_str(a.begin(), a.end()); };
inline std::string array2d_to_str(const array2d& a) { return array2d_to_str(a.begin(), a.end()); };

template<typename DegOrderedList>
inline bool deg_in(const Deg& deg, const DegOrderedList& basis) {
	auto pEle = std::lower_bound(basis.begin(), basis.end(), deg, [](const typename DegOrderedList::value_type& ele, const Deg& deg) {
		return ele.deg < deg; });
	return pEle < basis.end() && pEle->deg == deg;
}

template<typename DegOrderedList>
inline auto& get_item_by_deg(DegOrderedList& basis, const Deg& deg)
{
	auto pEle = std::lower_bound(basis.begin(), basis.end(), deg, [](const typename DegOrderedList::value_type& ele, const Deg& deg) { 
		return ele.deg < deg; });
	return *pEle;
}

/*--------- generate_basis.cpp ---------*/

/* Execute simple commands.*/
void execute_cmd(sqlite3* conn, const std::string& cmd);
int main_generate_basis(int argc, char** argv);
void load_gen_degs(sqlite3* conn, const std::string& table_name, std::vector<Deg>& gen_degs);
void load_basis(sqlite3* conn, const std::string& table_name, std::map<Deg, array2d>& basis);
void load_leading_terms(sqlite3* conn, const std::string& table_name, array3d& leadings);

/*--------- generate_mon_diffs.cpp ---------*/

int main_generate_diff(int argc, char** argv);
void load_gen_diffs(sqlite3* conn, const std::string& table_name, array3d& diffs);
void load_gb(sqlite3* conn, const std::string& table_name, array3d& gb);
inline std::pair<int, int> get_ab(const std::vector<std::array<int, 5>>& indices_tsv, int s, int t, int v)
{
	auto row = std::lower_bound(indices_tsv.begin(), indices_tsv.end(), std::array<int, 5>({ t, s, v, 0, 0 }));
	return std::pair<int, int>({ (*row)[3], (*row)[4] });
}

/*--------- generate_ss.cpp ---------*/

int main_generate_ss(int argc, char** argv);

/* generate_next_page.cpp */

int main_generate_next_page(int argc, char** argv);
array poly_to_indices(const array2d& poly, const array2d& basis);
array2d indices_to_poly(const array& indices, const array2d& basis);

inline int get_index(const array2d& basis, const array& mon)
{
	auto index = std::lower_bound(basis.begin(), basis.end(), mon, cmp_mons);
#ifdef _DEBUG
	if (index == basis.end()) {
		std::cout << "index not found\n";
		throw "178905cf-781e-4b4b-b981-2e11a3e98ed3";
	}
#endif
	return int(index - basis.begin());
}

/*--------- generate_E4t.cpp ---------*/

int main_generate_E4t(int argc, char** argv);
int main_test(int argc, char** argv);

#endif /* DATABSE_H */