#ifndef MAIN_H
#define MAIN_H

#include "groebner.h"
#include "database.h"
#include "myexception.h"

constexpr auto kLevelMin = 4;

/* A custom Exception class */
class SSException : public MyException
{
public:
	SSException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
	SSException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

/* Never be thrown. A dummy exception to catch */
class NoException : public MyException
{
public:
	NoException(unsigned int error_id, const char* message) : MyException(error_id, message) {}
	NoException(unsigned int error_id, const std::string& message) : MyException(error_id, message) {}
};

struct NullDiff { int index, num_tgts, first_r; };

class CacheDeduction
{
public:
	std::vector<std::map<Deg, NullDiff>> null_diffs_;
	std::vector<Deg> degs_;
	std::vector<int> indices_;
public:
	const NullDiff& GetRecentNullDiff(const Deg& deg) const;
	void push_back() { null_diffs_.push_back({}); degs_.push_back(Deg{ -1, -1, -1 }); indices_.push_back(-1); }
	void pop_back() { null_diffs_.pop_back(); degs_.pop_back(); indices_.pop_back(); }
	void InitNullDiffs(const std::vector<std::map<Deg, Staircase>>& basis_ss);
};

/* Staircase utilities */
const Staircase& GetRecentStaircase(const std::vector<std::map<Deg, Staircase>>& basis_ss, const Deg& deg);
void ApplyChanges(std::vector<std::map<Deg, Staircase>>& basis_ss);
std::ostream& operator<<(std::ostream& sout, const Staircase& sc);
size_t GetFirstIndex(const Staircase& sc, int level);
std::pair<int, int> CountUndetermined(const Staircase& sc, int level_min);
bool IsSingle(const Staircase& sc, int index);

/* Management of the data structure of spectral sequences */
void AddImage(std::vector<std::map<Deg, Staircase>>& basis_ss, const Deg& deg_dx, array dx, array x, int r, bool bForSSToF2 = false);
void AddDiff(std::vector<std::map<Deg, Staircase>>& basis_ss, const Deg& deg_x, array x, array dx, int r, bool bForSSToF2 = false);
void SetDiff(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, std::vector<std::map<Deg, Staircase>>& basis_ss, Deg deg_x, array x, array dx, int r, bool bForSSToF2 = false);

void generate_basis(const Database& db, const std::string& table_prefix, int t_max, bool drop_existing = false);
void generate_mon_diffs(const Database& db, const std::string& table_prefix, int r);
void generate_ss(const Database& db, const std::string& table_basis, const std::string& table_ss, int r);
void generate_next_page(const Database& db, const std::string& table_prefix, const std::string& table_H_prefix, int r);


#endif
