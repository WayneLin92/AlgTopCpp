#ifndef DATABSE_H
#define DATABSE_H

#include "mymath.h"
#include "sqlite3/sqlite3.h"
#include <iostream>
#include <vector>
#include <array>

constexpr auto N = 10000;
inline const char* sqlite3_column_str(sqlite3_stmt* stmt, int iCol) { return reinterpret_cast<const char*>(sqlite3_column_text(stmt, iCol)); }
inline int sqlite3_bind_str(sqlite3_stmt* stmt, int iCol, const std::string& str) { return sqlite3_bind_text(stmt, iCol, str.c_str(), -1, SQLITE_TRANSIENT); }
inline int sqlite3_prepare_v100(sqlite3* db, const char* zSql, sqlite3_stmt** ppStmt) { return sqlite3_prepare_v2(db, zSql, int(strlen(zSql)) + 1, ppStmt, NULL); }
inline int sqlite3_prepare_v100(sqlite3* db, std::string sql, sqlite3_stmt** ppStmt) { return sqlite3_prepare_v2(db, sql.c_str(), int(sql.size()) + 1, ppStmt, NULL); }

/* database.cpp */
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

/* generate_basis.cpp */
void execute_cmd(sqlite3* conn, const char* cmd);
/* num_mons[t] is the number of monomials up to degree t */
void load_mon_indices(sqlite3* conn, const char* table_name, array& num_mons);
int main_generate_basis(int argc, char** argv);

/* generate_diff.cpp */
void load_relations(sqlite3* conn, const std::string& table_name, array3d& relations);
void load_mon_indices_tsv(sqlite3* conn, const std::string& table_name, std::vector<std::array<int, 5>>& num_mons_tsv);
inline std::pair<int, int> get_ab(const std::vector<std::array<int, 5>>& indices_tsv, int s, int t, int v)
{
	auto row = std::lower_bound(indices_tsv.begin(), indices_tsv.end(), std::array<int, 5>({ t, s, v, 0, 0 }));
	return std::pair<int, int>({ (*row)[3], (*row)[4] });
}
int main_generate_diff(int argc, char** argv);

/* generate_ss.cpp */
int main_generate_ss(int argc, char** argv);

/* generate_next_page.cpp */
int main_generate_next_page(int argc, char** argv);

#endif /* DATABSE_H */