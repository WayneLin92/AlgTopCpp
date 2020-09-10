#include "main.h"

void generate_ss(const Database& db, const std::string& table_name_basis, const std::string& table_ss, int r)
{
	
	db.execute_cmd("DELETE FROM " + table_ss + ";");

	std::map<Deg, array2d> mon_diffs_ind = db.load_mon_diffs_ind(table_name_basis);
	std::map<Deg, BasisSS> basis_ss;

	/* fill basis_ss */
	int prev_t = 0;
	for (auto& [deg, mon_diffs_d] : mon_diffs_ind) {
		if (deg.t > prev_t) {
			prev_t = deg.t;
			std::cout << "Compute basis_ss, t=" << deg.t << "          \r";
		}
		BasisSS& basis_ss_d = basis_ss[deg];

		array range_ = range((int)mon_diffs_d.size());

		array lead_image;
		for (const auto& base : basis_ss_d.basis_ind)
			lead_image.push_back(base.front());
		std::sort(lead_image.begin(), lead_image.end());

		array x = add_vectors(range_, lead_image);

		array2d fx;
		for (int xi : x)
			fx.push_back(mon_diffs_d[xi]);

		array2d image, kernel, g;
		set_linear_map(x, fx, image, kernel, g);

		/* fill with other cycles after the boundaries */
		array lead_kernel;
		for (auto& cycle : kernel) {
			lead_kernel.push_back(cycle.front());
			basis_ss_d.basis_ind.push_back(cycle);
			basis_ss_d.diffs_ind.push_back({});
			basis_ss_d.levels.push_back(T_MAX - r);
		}
		std::sort(lead_kernel.begin(), lead_kernel.end());

		/* fill with boundaries */
		Deg deg_diff = deg + Deg{ 1, 0, -r };
		for (auto& boundary : image) {
			basis_ss[deg_diff].basis_ind.push_back(std::move(boundary));
			basis_ss[deg_diff].diffs_ind.push_back({});
			basis_ss[deg_diff].levels.push_back(r);
		}

		/* fill with the rest */

		array rest = add_vectors(x, lead_kernel);
		for (int i : rest) {
			basis_ss_d.basis_ind.push_back({ i });
			basis_ss_d.diffs_ind.push_back(std::move(mon_diffs_ind[deg][i]));
			basis_ss_d.levels.push_back(T_MAX);
		}
	}

	/* insert into the database */
	db.save_basis_ss(table_ss, basis_ss);
}

int main_generate_ss(int argc, char** argv)
{
	Database db;
	db.init(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	const char* table_basis = "E2_basis", * table_ss = "E2_ss";
	int t_max = 74;
	generate_ss(db, table_basis, table_ss, 2);
	return 0;
}