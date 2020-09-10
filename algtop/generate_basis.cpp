#include "main.h"

void generate_basis(const Database& db, const std::string& table_prefix, int t_max, bool drop_existing=false)
{
	/* load gen_degs, leadings, basis */
	if (drop_existing)
		db.execute_cmd("DELETE FROM " + table_prefix + "_basis;");
	std::vector<Deg> gen_degs = db.load_gen_degs(table_prefix + "_generators");
	array3d leadings = db.load_leading_terms(table_prefix + "_relations");
	std::map<Deg, array2d> basis = db.load_basis(table_prefix + "_basis");

	/* starting t value */
	int t_min;
	if (!basis.empty())
		t_min = basis.rbegin()->first.t + 1;
	else {
		basis[Deg{ 0, 0, 0 }].push_back({}); /* if no monomial present insert the unit */
		t_min = 1;
		db.execute_cmd("INSERT INTO " + table_prefix + "_basis (mon_id, mon, s, t, v) VALUES (0, \"\", 0, 0, 0);");
	}

	/* Add new basis */
	for (int t = t_min; t <= t_max; t++) {
		std::map<Deg, array2d> basis_new;
		std::cout << "t=" << t << "          \r";
		for (int gen_id = (int)gen_degs.size() - 1; gen_id >= 0; --gen_id) {
			int t1 = t - gen_degs[gen_id].t;
			if (t1 >= 0) {
				auto p1 = basis.lower_bound(Deg{ 0, t1, 0 });
				auto p2 = basis.lower_bound(Deg{ 0, t1 + 1, 0 });
				for (auto p = p1; p != p2; ++p) {
					for (const auto& m : p->second) {
						if (m.empty() || gen_id <= m[0]) {
							array mon(mul(m, { gen_id, 1 }));
							if ((size_t)gen_id >= leadings.size() || std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(),
								[&mon](const array& _m) { return divides(_m, mon); }))
								basis_new[p->first + gen_degs[gen_id]].push_back(std::move(mon));
						}
					}
				}
			}
		}
		/* Insert the new basis in degree t into the database */
		db.save_basis(table_prefix + "_basis", basis_new);
		basis.merge(basis_new);
	}
}

int main_generate_basis(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	std::string table_prefix = "E4b6";
	int t_max = 74;
	generate_basis(db, table_prefix, t_max, true);
	return 0;
}