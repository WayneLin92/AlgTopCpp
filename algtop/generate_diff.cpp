#include "main.h"

void generate_mon_diffs(const Database& db, const std::string& table_prefix, int r)
{
	/* load gen_degs, gen_diffs, gb and basis */

	std::vector<Deg> gen_degs = db.load_gen_degs(table_prefix + "_generators");
	array3d gen_diffs = db.load_gen_diffs(table_prefix + "_generators");
	array3d gb = db.load_gb(table_prefix + "_relations");
	std::map<Deg, array2d> basis = db.load_basis(table_prefix + "_basis");
	std::map<Deg, array3d> mon_diffs = db.load_mon_diffs(table_prefix + "_basis", basis, r);

	/* compute diffs */

	Statement stmt;
	stmt.init(db, "UPDATE " + table_prefix + "_basis SET diff=?1 WHERE mon=?2;");

	db.begin_transaction();
	int prev_t = 0;
	for (auto& [d, basis_d] : basis) {
		if (mon_diffs.find(d) != mon_diffs.end())
			continue;
		if (d.t > prev_t) {
			prev_t = d.t;
			std::cout << "computing diffs, t=" << d.t << "          \r";
			db.end_transaction();
			db.begin_transaction();
		}
		for (size_t i = 0; i < basis_d.size(); ++i) {
			array diff_indices;
			if (basis_d[i].empty())
				mon_diffs[d].push_back({});
			else {
				int gen_id = basis_d[i].front();
				array mon1 = div(basis_d[i], { gen_id, 1 });
				Deg d1 = d - gen_degs[gen_id];
				size_t index_mon1 = std::lower_bound(basis.at(d1).begin(), basis.at(d1).end(), mon1, cmp_mons) - basis.at(d1).begin();
				mon_diffs[d].push_back(reduce(add(mul(gen_diffs[gen_id], mon1), mul({ { gen_id, 1 } }, mon_diffs[d1][index_mon1])), gb));
				Deg d_diff = d + Deg{ 1, 0, -r };
				diff_indices = poly_to_indices(mon_diffs[d][i], basis[d_diff]);
			}
			std::string str_diff(array_to_str(diff_indices));
			stmt.bind_str(1, str_diff);
			stmt.bind_str(2, array_to_str(basis_d[i]));
			stmt.step_and_reset();
		}
	}
	db.end_transaction();
}

int test1()
{
	return 0;
}

int main_generate_diff(int argc, char** argv)
{
	//return test1();

	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	std::string table_prefix = "E2";
	generate_mon_diffs(db, table_prefix, 2);
	return 0;
}