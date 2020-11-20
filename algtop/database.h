#ifndef DATABSE_H
#define DATABSE_H

#include "mymath.h"
#include "sqlite3/sqlite3.h"
#include <map>

constexpr auto T_MAX = 10000;

struct BasisComplex
{
	array2d boundaries;
	array2d cycles;
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
	Database() : conn_(nullptr) {}
	~Database() { sqlite3_close(conn_); }
	void init(const char* filename) { if (sqlite3_open(filename, &conn_) != SQLITE_OK) throw "8de81e80"; }
public:
	void execute_cmd(const std::string& sql) const;
	int get_int(const std::string& sql) const;
	array get_ints(const std::string& table_name, const std::string& column_name, const std::string& conditions="") const;
	void sqlite3_prepare_v100(const char* zSql, sqlite3_stmt** ppStmt) const { if (sqlite3_prepare_v2(conn_, zSql, int(strlen(zSql)) + 1, ppStmt, NULL) != SQLITE_OK) throw "bce2dcfe"; }
	void sqlite3_prepare_v100(const std::string& sql, sqlite3_stmt** ppStmt) const { if (sqlite3_prepare_v2(conn_, sql.c_str(), int(sql.size()) + 1, ppStmt, NULL) != SQLITE_OK) throw "da6ab7f6"; }
public:
	void begin_transaction() const { execute_cmd("BEGIN TRANSACTION"); }
	void end_transaction() const { execute_cmd("END TRANSACTION"); }
public:
	std::vector<Deg> load_gen_degs(const std::string& table_name) const;
	Poly1d load_gen_diffs(const std::string& table_name) const;
	Poly1d load_gen_reprs(const std::string& table_name) const;
	Mon2d load_leading_terms(const std::string& table_name) const;
	Poly1d load_gb(const std::string& table_name) const;
	std::map<Deg, Mon1d> load_basis(const std::string& table_name) const;
	std::map<Deg, array2d> load_mon_diffs_ind(const std::string& table_name) const;
	std::map<Deg, Poly1d> load_mon_diffs(const std::string& table_name, const std::map<Deg, Mon1d>& basis, int r) const;
	std::map<Deg, BasisComplex> load_basis_ss(const std::string& table_name_ss, int r) const;
public:
	void save_generators(const std::string& table_name, const std::vector<Deg>& gen_degs, const Poly1d& gen_reprs, size_t i_start) const;
	void save_gb(const std::string& table_name, const Poly1d& gb, const std::vector<Deg>& gen_degs, size_t i_start) const;

	void save_generators(const std::string& table_name, const std::vector<Deg>& gen_degs, const Poly1d& gen_reprs) const { save_generators(table_name, gen_degs, gen_reprs, 0); }
	void save_gb(const std::string& table_name, const Poly1d& gb, const std::vector<Deg>& gen_degs) const { save_gb(table_name, gb, gen_degs, 0); }
	void save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis) const;
	void save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis, const std::map<Deg, array2d>& mon_reprs) const;
	void save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis, const std::map<Deg, Poly1d>& mon_reprs) const;
	void save_basis_ss(const std::string& table_name, const std::map<Deg, BasisSS>& basis_ss) const;
private:
	sqlite3* conn_;
};

class Statement
{
public:
	Statement() : stmt_(nullptr) {}
	~Statement() { sqlite3_finalize(stmt_); }
	void init(const Database& db, const std::string& sql) { db.sqlite3_prepare_v100(sql, &stmt_); }
public:
	void bind_str(int iCol, const std::string& str) const { if (sqlite3_bind_text(stmt_, iCol, str.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) throw "29cc3b21"; }
	void bind_int(int iCol, int i) const { if (sqlite3_bind_int(stmt_, iCol, i) != SQLITE_OK) throw "a61e05b2"; }
	const char* column_str(int iCol) const { return reinterpret_cast<const char*>(sqlite3_column_text(stmt_, iCol)); }
	int column_int(int iCol) const { return sqlite3_column_int(stmt_, iCol); }
	int column_type(int iCol) const { return sqlite3_column_type(stmt_, iCol); }
	int step() const { return sqlite3_step(stmt_); }
	int reset() const { return sqlite3_reset(stmt_); }
	void step_and_reset() const { step(); reset(); }
private:
	sqlite3_stmt* stmt_;
};

/*--------- my_utilities.cpp ---------*/

std::string array_to_str(array::const_iterator pbegin, array::const_iterator pend);
array str_to_array(const char* str_mon);
std::string Mon_to_str(MonInd pbegin, MonInd pend);
Mon str_to_Mon(const char* str_mon);
std::string Poly_to_str(Poly::const_iterator pbegin, Poly::const_iterator pend);
Poly str_to_Poly(const char* str_poly);

inline std::string array_to_str(const array& a) { return array_to_str(a.begin(), a.end()); };
inline std::string Mon_to_str(const Mon& mon) { return Mon_to_str(mon.begin(), mon.end()); };
inline std::string Poly_to_str(const Poly& poly) { return Poly_to_str(poly.begin(), poly.end()); };

inline int get_index(const Poly& basis, const Mon& mon)
{
	auto index = std::lower_bound(basis.begin(), basis.end(), mon);
#ifdef _DEBUG
	if (index == basis.end()) {
		std::cout << "index not found\n";
		throw "178905cf";
	}
#endif
	return int(index - basis.begin());
}

array Poly_to_indices(const Poly& poly, const Poly& basis);
Poly indices_to_Poly(const array& indices, const Poly& basis);


#endif /* DATABSE_H */