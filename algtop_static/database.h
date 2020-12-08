#ifndef DATABSE_H
#define DATABSE_H

//#define DATABASE_SAVE_LOGGING /* This is a switch to print what are saved to the database */

#include "algebras.h"
#include <string>
#include <vector>
#include <map>

/********** STRUCTS AND CLASSES **********/
/* Declarations */
constexpr auto T_MAX = 10000;
struct sqlite3;
struct sqlite3_stmt;

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

/* The wrapper for sqlite* */
class Database
{
public:
	Database() : conn_(nullptr) {}
	~Database();
	void init(const char* filename);
public:
	void execute_cmd(const std::string& sql) const;
	int get_int(const std::string& sql) const;
	array get_ints(const std::string& table_name, const std::string& column_name, const std::string& conditions="") const;
	void sqlite3_prepare_v100(const char* zSql, sqlite3_stmt** ppStmt) const;
	void sqlite3_prepare_v100(const std::string& sql, sqlite3_stmt** ppStmt) const;
public:
	void begin_transaction() const { execute_cmd("BEGIN TRANSACTION"); }
	void end_transaction() const { execute_cmd("END TRANSACTION"); }
public:
	std::vector<Deg> load_gen_degs(const std::string& table_name) const;
	Poly1d load_gen_diffs(const std::string& table_name) const;
	Poly1d load_gen_reprs(const std::string& table_name) const;
	Mon2d load_leading_terms(const std::string& table_name, int t_max) const;
	Poly1d load_gb(const std::string& table_name, int t_max) const;
	std::map<Deg, Mon1d> load_basis(const std::string& table_name, int t_max) const;
	std::map<Deg, array2d> load_mon_diffs_ind(const std::string& table_name, int t_max) const;
	std::map<Deg, Poly1d> load_mon_diffs(const std::string& table_name, const std::map<Deg, Mon1d>& basis, int r, int t_max) const;
	std::map<Deg, BasisComplex> load_basis_ss(const std::string& table_name_ss, int r, int t_max) const;
public:
	void save_generators(const std::string& table_name, const std::vector<Deg>& gen_degs, const Poly1d& gen_reprs, size_t i_start) const;
	void save_generators(const std::string& table_name, const std::vector<Deg>& gen_degs, const Poly1d& gen_reprs) const { save_generators(table_name, gen_degs, gen_reprs, 0); }
	void save_gb(const std::string& table_name, const Poly1d& gb, const std::vector<Deg>& gen_degs, size_t i_start) const;
	void save_gb(const std::string& table_name, const Poly1d& gb, const std::vector<Deg>& gen_degs) const { save_gb(table_name, gb, gen_degs, 0); }
	void save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis) const;
	void save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis, const std::map<Deg, array2d>& mon_reprs) const;
	void save_basis(const std::string& table_name, const std::map<Deg, Mon1d>& basis, const std::map<Deg, Poly1d>& mon_reprs) const;
	void save_basis_ss(const std::string& table_name, const std::map<Deg, BasisSS>& basis_ss) const;
private:
	sqlite3* conn_;
};

/* The wrapper for sqlite3_stmt* */
class Statement
{
public:
	Statement() : stmt_(nullptr) {}
	~Statement();
	void init(const Database& db, const std::string& sql);
public:
	void bind_str(int iCol, const std::string& str) const;
	void bind_int(int iCol, int i) const;
	const char* column_str(int iCol) const;
	int column_int(int iCol) const;
	int column_type(int iCol) const;
	int step() const;
	int reset() const;
	void step_and_reset() const;
private:
	sqlite3_stmt* stmt_;
};

array str_to_array(const char* str_mon);
Mon str_to_Mon(const char* str_mon);
Poly str_to_Poly(const char* str_poly);
std::string array_to_str(array::const_iterator pbegin, array::const_iterator pend);
inline std::string array_to_str(const array& a) { return array_to_str(a.begin(), a.end()); }
std::string Mon_to_str(MonInd pbegin, MonInd pend);
inline std::string Mon_to_str(const Mon& mon) { return Mon_to_str(mon.begin(), mon.end()); }
std::string Poly_to_str(Poly::const_iterator pbegin, Poly::const_iterator pend);
inline std::string Poly_to_str(const Poly& poly) { return Poly_to_str(poly.begin(), poly.end()); }

Poly indices_to_Poly(const array& indices, const Mon1d& basis);
array Poly_to_indices(const Poly& poly, const Mon1d& basis);

#endif /* DATABSE_H */
