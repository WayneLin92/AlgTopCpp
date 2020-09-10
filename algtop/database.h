#ifndef DATABSE_H
#define DATABSE_H

#include "mymath.h"
#include "sqlite3/sqlite3.h"
#include <map>

constexpr auto T_MAX = 10000;

struct BasisComplex
{
	array2d boundary;
	array2d kernel;
};

struct BasisSS
{
	array2d basis_ind;
	array2d diffs_ind;
	array levels;
};

class Database
{
public:
	Database() : m_conn(nullptr) {}
	~Database() { sqlite3_close(m_conn); }
	void init(const char* filename) { if (sqlite3_open(filename, &m_conn) != SQLITE_OK) throw "8de81e80"; }
public:
	void execute_cmd(const std::string& cmd) const;
	void sqlite3_prepare_v100(const char* zSql, sqlite3_stmt** ppStmt) const { if (sqlite3_prepare_v2(m_conn, zSql, int(strlen(zSql)) + 1, ppStmt, NULL) != SQLITE_OK) throw "bce2dcfe"; }
	void sqlite3_prepare_v100(const std::string& sql, sqlite3_stmt** ppStmt) const { if (sqlite3_prepare_v2(m_conn, sql.c_str(), int(sql.size()) + 1, ppStmt, NULL) != SQLITE_OK) throw "da6ab7f6"; }
public:
	void begin_transaction() const { execute_cmd("BEGIN TRANSACTION"); }
	void end_transaction() const { execute_cmd("END TRANSACTION"); }
public:
	std::vector<Deg> load_gen_degs(const std::string& table_name) const;
	array3d load_gen_diffs(const std::string& table_name) const;
	array3d load_gen_reprs(const std::string& table_name) const;
	array3d load_leading_terms(const std::string& table_name) const;
	array3d load_gb(const std::string& table_name) const;
	std::map<Deg, array2d> load_basis(const std::string& table_name) const;
	std::map<Deg, array2d> load_mon_diffs_ind(const std::string& table_name) const;
	std::map<Deg, array3d> load_mon_diffs(const std::string& table_name, const std::map<Deg, array2d>& basis, int r) const;
	std::map<Deg, BasisComplex> load_basis_ss(const std::string& table_name_ss, int r) const;
public:
	void save_generators(const std::string& table_name, const std::vector<Deg>& gen_degs, const array3d& gen_reprs) const;
	void save_gb(const std::string& table_name, const array3d& gb, const std::vector<Deg>& gen_degs) const;
	void save_basis(const std::string& table_name, const std::map<Deg, array2d>& basis) const;
	void save_basis(const std::string& table_name, const std::map<Deg, array2d>& basis, const std::map<Deg, array2d>& mon_reprs) const;
	void save_basis_ss(const std::string& table_name, const std::map<Deg, BasisSS>& basis_ss) const;
private:
	sqlite3* m_conn;
};

class Statement
{
public:
	Statement() : m_stmt(nullptr) {}
	~Statement() { sqlite3_finalize(m_stmt); }
	void init(const Database& db, const std::string& cmd) { db.sqlite3_prepare_v100(cmd, &m_stmt); }
public:
	void bind_str(int iCol, const std::string& str) const { if (sqlite3_bind_text(m_stmt, iCol, str.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) throw "29cc3b21"; }
	void bind_int(int iCol, int i) const { if (sqlite3_bind_int(m_stmt, iCol, i) != SQLITE_OK) throw "a61e05b2"; }
	const char* column_str(int iCol) const { return reinterpret_cast<const char*>(sqlite3_column_text(m_stmt, iCol)); }
	int column_int(int iCol) const { return sqlite3_column_int(m_stmt, iCol); }
	int step() const { return sqlite3_step(m_stmt); }
	int reset() const { return sqlite3_reset(m_stmt); }
	void step_and_reset() const { step(); reset(); }
private:
	sqlite3_stmt* m_stmt;
};

/*--------- my_utilities.cpp ---------*/

std::string array_to_str(array::const_iterator pbegin, array::const_iterator pend);
array str_to_array(const char* str_mon);
std::string array2d_to_str(array2d::const_iterator pbegin, array2d::const_iterator pend);
array2d str_to_array2d(const char* str_poly);

inline std::string array_to_str(const array& a) { return array_to_str(a.begin(), a.end()); };
inline std::string array2d_to_str(const array2d& a) { return array2d_to_str(a.begin(), a.end()); };

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

array poly_to_indices(const array2d& poly, const array2d& basis);
array2d indices_to_poly(const array& indices, const array2d& basis);

inline const char* sqlite3_column_str(sqlite3_stmt* stmt, int iCol) { return reinterpret_cast<const char*>(sqlite3_column_text(stmt, iCol)); }//
inline int sqlite3_bind_str(sqlite3_stmt* stmt, int iCol, const std::string& str) { return sqlite3_bind_text(stmt, iCol, str.c_str(), -1, SQLITE_TRANSIENT); }//
inline int sqlite3_prepare_v100(sqlite3* db, const char* zSql, sqlite3_stmt** ppStmt) { return sqlite3_prepare_v2(db, zSql, int(strlen(zSql)) + 1, ppStmt, NULL); }//
inline int sqlite3_prepare_v100(sqlite3* db, std::string sql, sqlite3_stmt** ppStmt) { return sqlite3_prepare_v2(db, sql.c_str(), int(sql.size()) + 1, ppStmt, NULL); }//


#endif /* DATABSE_H */