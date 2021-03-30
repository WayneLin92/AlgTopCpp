#ifndef MAIN_H
#define MAIN_H

#include "groebner.h"
#include "database.h"
#include "myexception.h"

constexpr auto kLevelMin = 4;
constexpr auto kLevelPC = kLevelMax / 4 * 3; /* Level of Permanant cycles */
enum struct DiffType
{
	E4 = 0,
	E4h0 = 1,
	E4t = 2,
};
using Staircases = std::map<Deg, Staircase>;
using Staircases1d = std::vector<Staircases>;

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

struct NullDiff { int index, num_tgts, first; };

class NullDiffs
{
public:
	std::map<Deg, NullDiff> null_diffs_; //TODO: This only captures a single null differential in the same degree
public:
	void InitNullDiffs(const Staircases1d& basis_ss, DiffType d_type, int t_max, bool bNew);
};

class CacheDeduction
{
public:
	std::vector<std::map<Deg, NullDiff>> null_diffs_; //TODO: This only captures a single null differential in the same degree
	std::vector<Deg> degs_;
	std::vector<int> indices_;
public:
	const NullDiff& GetRecentNullDiff(const Deg& deg) const;
	void push_back() { null_diffs_.push_back({}); degs_.push_back(Deg{ -1, -1, -1 }); indices_.push_back(-1); } //TODO: loop index backwards vs forwards
	void pop_back() { null_diffs_.pop_back(); degs_.pop_back(); indices_.pop_back(); }
	void InitNullDiffs(const Staircases1d& basis_ss, DiffType d_type, int t_max);
};

template<typename GbType>
Poly get_image(const Poly& poly, const Poly1d& gen_images, const GbType& gb)
{
	return grbn::evaluate(poly, [&gen_images](int i) {return gen_images[i]; }, gb);
}

inline std::pair<int, int> GetST(DiffType diff_type)
{
	switch (diff_type)
	{
	case DiffType::E4:
	case DiffType::E4t:
		return { 1, 0 };
	case DiffType::E4h0:
		return { 0, -1 };
	default:
		throw MyException(0x13d57cb1U, "unhandle difftype");
	}
};

DiffType GetDiffType(const std::string& table_prefix);
bool IsEt(const std::string& table_prefix);
bool IsEt(DiffType d_type);

/* Staircase utilities */
const Staircase& GetRecentStaircase(const Staircases1d& basis_ss, const Deg& deg);
void ApplyAllChanges(Staircases1d& basis_ss, size_t nHistory = 0);
void ApplyRecentChanges(Staircases1d& basis_ss);
std::ostream& operator<<(std::ostream& sout, const Staircase& sc);
bool NoNullDiff(const Staircases1d& basis_ss, int t);
NullDiff GetNextNullDiff(const Staircases1d& basis_ss, const Deg& deg, DiffType d_type, int index);
std::tuple<int, Deg, int> CountTgt(const Staircases1d& basis_ss, const Deg& deg, DiffType d_type, int r);
std::tuple<int, Deg, int> CountSrc(const Staircases1d& basis_ss, const Deg& deg, DiffType d_type, int level);

/* Management of the data structure of spectral sequences */
void AddImage(Staircases1d& basis_ss, const Deg& deg_dx, DiffType d_type, array dx, array x, int r);
void AddDiff(Staircases1d& basis_ss, const Deg& deg_x, DiffType d_type, array x, array dx, int r);
void SetDiff(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, Deg deg_x, DiffType d_type, array x, array dx, int r, int t_max = -1);
int SetDiffV2(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, const Deg& deg_x, const Deg& deg_dx, DiffType d_type, array x, array dx, int r, int t_max = -1);
void SetImage(const grbn::GbWithCache& gb, const std::map<Deg, Mon1d>& basis, Staircases1d& basis_ss, Deg deg_x, DiffType d_type, array x, int r, int t_max = -1);

void DeduceDiffs(const Database& db, const std::string& table_prefix, int t_max);

void generate_basis(const Database& db, const std::string& table_prefix, int t_max, bool drop_existing);
void generate_mon_diffs(const Database& db, const std::string& table_prefix, int r);
void generate_ss(const Database& db, const std::string& table_basis, const std::string& table_ss, int r);
void generate_next_page(const Database& db, const std::string& table_prefix, const std::string& table_H_prefix, int r);

void WrapDeduceZeroDiffs(const Database& db, const std::string& table_prefix, int t_max);
void WrapDeduceDiffsForEt(const Database& db, const std::string& table_prefix, int t_max);
void WrapDeduceDiffs(const Database& db, const std::string& table_prefix, int t_try, int t_test, int max_poss, int nCheck, double seconds);
void WrapDeduceDiffsByNat(const Database& db, const std::string& table_prefix, const std::string& table_prefix1, std::string& column,
	int t_try, int t_test);
void WrapDeduceDiffsByNatV2(const Database& db, const std::string& table_prefix, const std::string& table_prefix1, std::string& column, int t_max);
void WrapDeduceDiffsForE4(const Database& db, int t_try, int t_test);

#endif
