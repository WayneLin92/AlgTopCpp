#pragma once
#include "../algtop/main.h"
#include <fstream>
#include <future>


/********** STRUCTS AND CLASSES **********/

struct DgaBasis1
{
	Mon1d basis;
	Poly1d diffs;
};

/********** FUNCTIONS **********/

inline Poly get_repr(const Poly& poly, const Poly1d& gen_reprs, const Poly1d& gb)
{
	return evaluate(poly, [&gen_reprs](int i) {return gen_reprs[i]; }, gb);
}

std::map<Deg, DgaBasis1> load_dga_basis(const Database& db, const std::string& table_name, int r);
std::map<Deg, DgaBasis1> load_basis_X(const Database& db, const std::string& table_name, int t_max, int r);
void load_y(const Database& db, const std::string& table_name, Poly1d& y, array& t_y);
void save_y(const Database& db, const std::string& table_name, const Poly2d& y_t, int t);
void save_map_gen_id(const Database& db, const std::string& table_name, const array& map_gen_id, int i_start);
void SaveGb(const Database& db, const std::string& table_name, const Poly1d& gb, const array& gen_degs, const array& gen_degs1, int t);
void get_basis_B(const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X, Mon1d& basis_B, const Deg& deg);
void get_basis_with_diff_B(const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X, Mon1d& basis_B, Poly1d& mon_diffs_B, const Deg& deg);
/* Assume poly is a boundary. Return the chain with it as boundary */
Poly d_inv(const Poly& poly, const std::vector<Deg>& gen_degs, const Poly1d& diffs, const Poly1d& gb, const std::map<Deg, DgaBasis1>& basis_A, const std::map<Deg, DgaBasis1>& basis_X);
/* Assume poly is a cycle. Return the homology class */
Poly proj(const Poly& poly, const std::vector<Deg>& gen_degs, const Poly1d& gen_diffs, const Poly1d& gb, const std::map<Deg, DgaBasis1>& basis_A,
	const std::map<Deg, DgaBasis1>& basis_X, const Poly1d& gen_reprs, std::map<Deg, Mon1d>& basis_H);
template <typename Fn>
void add_rels_freemodule(Poly1d& gb, RelHeap& heap, Fn get_deg, int deg, int deg_max);

int main_generate_X_basis(int argc, char** argv);