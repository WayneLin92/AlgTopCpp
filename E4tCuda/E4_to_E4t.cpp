#include "main_E4t.h"

void generate_map_E4_E4t()
{
	Database db_E4t(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t.db)");
	Database db_tmp(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	int t_max = 189;
	std::vector<Deg> gen_degs_E2t = db_E4t.load_gen_degs("B_generators");
	Poly1d gen_diffs_E2t = db_E4t.load_gen_diffs("B_generators");
	grbn::GbWithCache gb_E2t = db_E4t.load_gb("A_relations", t_max);
	const std::map<Deg, DgaBasis1> basis_E2 = load_dga_basis(db_E4t, "A_basis", 2, t_max);
	std::map<Deg, DgaBasis1> basis_X = load_basis_X(db_E4t, "X3_basis", t_max, 2);
	Poly1d gen_reprs_E4t = db_tmp.load_gen_reprs("E4t_generators");
	std::map<Deg, Mon1d> basis_E4t = db_tmp.load_basis("E4t_basis", t_max);
	
	Poly1d gen_reprs_E4 = db_tmp.load_gen_reprs("E4_generators");
	Poly1d image_E4_E4t;
	for (const auto& p : gen_reprs_E4) {
		if (get_deg_t(p, gen_degs_E2t) > t_max)
			break;
		Poly image_p = proj(p, gen_degs_E2t, gen_diffs_E2t, gb_E2t, basis_E2, basis_X, gen_reprs_E4t, basis_E4t);
		image_E4_E4t.push_back(std::move(image_p));
	}

	db_tmp.save_gen_images("E4_generators", "to_E4t", image_E4_E4t);
}

void generate_map_E4_E4t_loc()
{
	Database db_E4t(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t.db)");
	Database db_tmp(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	int t_max = 230;
	int t_max_E2t = 189;
	std::vector<Deg> gen_degs_E2t = db_E4t.load_gen_degs("B_generators");
	Poly1d gen_diffs_E2t = db_E4t.load_gen_diffs("B_generators");
	grbn::GbWithCache gb_E2t = db_E4t.load_gb("A_relations", t_max_E2t);
	const std::map<Deg, DgaBasis1> basis_E2 = load_dga_basis(db_E4t, "A_basis", 2, t_max_E2t);
	std::map<Deg, DgaBasis1> basis_X = load_basis_X(db_E4t, "X3_basis", t_max_E2t, 2);
	Poly1d gen_reprs_E4t = db_tmp.load_gen_reprs("E4t_generators");
	std::map<Deg, Mon1d> basis_E4t = db_tmp.load_basis("E4t_basis", t_max_E2t);

	Poly1d image_E2_loc = db_tmp.load_gen_images("E2_generators", "loc", -1);
	Poly1d gen_reprs_E4 = db_tmp.load_gen_reprs("E4_generators");
	Poly1d image_E4_E4t;
	for (const auto& p : gen_reprs_E4) {
		if (get_deg_t(p, gen_degs_E2t) > t_max)
			break;
		Poly p1 = get_image(p, image_E2_loc, gb_E2t);
		Poly image_p = proj(p1, gen_degs_E2t, gen_diffs_E2t, gb_E2t, basis_E2, basis_X, gen_reprs_E4t, basis_E4t);
		image_E4_E4t.push_back(std::move(image_p));
	}

	db_tmp.save_gen_images("E4_generators", "loc", image_E4_E4t);
}

void generate_map_E4_E4t_loc_sq0()
{
	Database db_E4t(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t.db)");
	Database db_tmp(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\tmp.db)");
	int t_max = 230;
	int t_max_E2t = 189;
	std::vector<Deg> gen_degs_E2t = db_E4t.load_gen_degs("B_generators");
	Poly1d gen_diffs_E2t = db_E4t.load_gen_diffs("B_generators");
	grbn::GbWithCache gb_E2t = db_E4t.load_gb("A_relations", t_max_E2t);
	const std::map<Deg, DgaBasis1> basis_E2 = load_dga_basis(db_E4t, "A_basis", 2, t_max_E2t);
	std::map<Deg, DgaBasis1> basis_X = load_basis_X(db_E4t, "X3_basis", t_max_E2t, 2);
	Poly1d gen_reprs_E4t = db_tmp.load_gen_reprs("E4t_generators");
	std::map<Deg, Mon1d> basis_E4t = db_tmp.load_basis("E4t_basis", t_max_E2t);

	Poly1d image_E2_loc = db_tmp.load_gen_images("E2_generators", "loc_sq0", -1);
	Poly1d gen_reprs_E4 = db_tmp.load_gen_reprs("E4_generators");
	Poly1d image_E4_E4t;
	for (const auto& p : gen_reprs_E4) {
		if (get_deg_t(p, gen_degs_E2t) > t_max)
			break;
		Poly p1 = get_image(p, image_E2_loc, gb_E2t);
		Poly image_p = proj(p1, gen_degs_E2t, gen_diffs_E2t, gb_E2t, basis_E2, basis_X, gen_reprs_E4t, basis_E4t);
		image_E4_E4t.push_back(std::move(image_p));
	}

	db_tmp.save_gen_images("E4_generators", "loc_sq0", image_E4_E4t);
}

int main_generate_E4_to_E4t(int argc, char** argv)
{
	generate_map_E4_E4t_loc();
	generate_map_E4_E4t_loc_sq0();
	return 0;
}