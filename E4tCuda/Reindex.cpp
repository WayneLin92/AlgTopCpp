#include "main_E4t.h"

Poly reindex_v2(const Poly& poly, array map_gen_id)
{
	Poly result;
	for (const Mon& m : poly) {
		Mon m1;
		for (GenPow ge : m)
			m1.push_back({ map_gen_id[ge.gen], ge.exp });
		std::sort(m1.begin(), m1.end(), [](const GenPow& lhs, const GenPow& rhs) {return lhs.gen < rhs.gen; });
		result.push_back(m1);
	}
	std::sort(result.begin(), result.end());
	return result;
}

/* This function reorders the generator by t and reproduce the Groebner basis */
std::pair<array, Poly1d> ReorderGens(const std::vector<Deg>& gen_degs, const Poly1d& gb, int t_max)
{
	array map_gen_id_inv = grbn::range((int)gen_degs.size()); /* the i`th new generator is the old map_gen_id_inv[i]`th generator */
	std::sort(map_gen_id_inv.begin(), map_gen_id_inv.end(), [&gen_degs](int i, int j) {return gen_degs[i].t < gen_degs[j].t; });
	array map_gen_id; map_gen_id.resize(gen_degs.size()); /* the i`th old generator becomes the map_gen_id[i]`th generator */
	array gen_degs_new;
	for (int i = 0; i < (int)gen_degs.size(); ++i) {
		map_gen_id[map_gen_id_inv[i]] = i;
		gen_degs_new.push_back(gen_degs[map_gen_id_inv[i]].t);
	}

	grbn::RelBuffer buffer;
	for (const Poly& g : gb)
		buffer[get_deg_t(g, gen_degs)].push_back(reindex_v2(g, map_gen_id));
	Poly1d gb_new;
	grbn::AddRelsB(gb_new, buffer, gen_degs_new, -1, t_max);
	return std::make_pair(std::move(map_gen_id_inv), std::move(gb_new));
}

void ReorderHA(const Database& db, int t_max)
{
	std::vector<Deg> gen_degs_HA = db.load_gen_degs("HA_generators");
	Poly1d gen_reprs_HA = db.load_gen_reprs("HA_generators");
	Poly1d gb_HA = db.load_gb("HA_relations", t_max);
	auto [map_gen_id_inv, gb_new] = ReorderGens(gen_degs_HA, gb_HA, t_max);

	std::vector<Deg> gen_degs_HA_new;
	Poly1d gen_reprs_HA_new;
	for (int i = 0; i < (int)gen_degs_HA.size(); ++i) {
		gen_degs_HA_new.push_back(gen_degs_HA[map_gen_id_inv[i]]);
		gen_reprs_HA_new.push_back(gen_reprs_HA[map_gen_id_inv[i]]);
	}

	db.execute_cmd("DELETE FROM HA_generators_ordered;");
	db.execute_cmd("DELETE FROM HA_relations_ordered;");

	db.save_generators("HA_generators_ordered", gen_degs_HA_new, gen_reprs_HA_new);
	db.save_gb("HA_relations_ordered", gb_new, gen_degs_HA_new);
}

int main_generate_Reindex(int argc, char** argv)
{
	Database db(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t.db)");

	ReorderHA(db, 200);
	return 0;
}