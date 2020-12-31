#include "groebner.h"
#include "benchmark.h"
#include "database.h"
#include "myio.h"
#include <iostream>

int leading_distribution()
{
	Timer timer;
	Database db(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t.db)");
	//std::vector<Deg> gen_degs = db.load_gen_degs("HA3_generators");
	Poly1d gb = db.load_gb("HA3_relations", 189);
	std::map<int, int> count_leads;
	for (const Poly& g : gb)
		++count_leads[g[0].back().gen];
	std::map<int, int> count_leads_sum;
	int sum = 0;
	for (auto [index, num] : count_leads) {
		sum += num;
		count_leads_sum[index] = sum;
	}
	for (auto [index, num] : count_leads_sum)
		std::cout << "index=" << index << ", num=" << num / 230688.0 << '\n';
	return 0;
}

int test_b8()
{
	Timer timer;

	array gen_degs;
	int n_max = 8;
	for (int d = 1; d <= n_max; d++) {
		for (int i = 0; i <= n_max - d; i++) {
			int j = i + d;
			gen_degs.push_back((1 << j) - (1 << i));
		}
}
	Poly1d rels;
	for (int d = 2; d <= n_max; d++) {
		for (int i = 0; i <= n_max - d; i++) {
			int j = i + d;
			Poly rel;
			for (int k = i + 1; k < j; k++) {
				int a = (1 << k) - (1 << i);
				int b = (1 << j) - (1 << k);
				auto p1 = std::find(gen_degs.begin(), gen_degs.end(), a);
				auto p2 = std::find(gen_degs.begin(), gen_degs.end(), b);
				int index1 = int(p1 - gen_degs.begin());
				int index2 = int(p2 - gen_degs.begin());
				rel = add(rel, index1 < index2 ? Poly{ {{index1, 1}, {index2, 1}} } : Poly{ {{index2, 1}, {index1, 1}} });
			}
			rels.push_back(std::move(rel));
		}
	}

	grbn::GbWithCache gb;
	std::sort(rels.begin(), rels.end(), [&gen_degs](const Poly& p1, const Poly& p2) {
		return get_deg(p1, gen_degs) < get_deg(p2, gen_degs); });
	grbn::AddRels(gb, std::move(rels), FnGetDeg{ gen_degs }, -1);
	size_t gb_size = gb.size();

	size_t answer = 163;
	std::cout << gb_size << '=' << answer << '\n';

	return 0;
}

int BenchmarkAddRels()
{
	Database db(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t.db)");
	std::string table_name = "E4t";
	int t_max = 130;

	std::vector<Deg> gen_degs = db.load_gen_degs(table_name + "_generators");
	array gen_degs_t;
	for (size_t i = 0; i < gen_degs.size(); ++i)
		gen_degs_t.push_back(gen_degs[i].t);
	Poly1d gb = db.load_gb(table_name + "_relations", t_max);

	Poly1d gb0, gb1, gb2;
	grbn::GbWithCache gb3, gb4;

	grbn::AddRels(gb0, gb, FnGetDeg{ gen_degs_t }, t_max);

	Timer timer;
	grbn::AddRels(gb1, gb, FnGetDeg{ gen_degs_t }, t_max);
	timer.print("Poly1d, GbBuffer: ");

	grbn::AddRels<Poly1d, FnGetDeg, grbn::GbBufferV2>(gb2, gb, FnGetDeg{ gen_degs_t }, t_max);
	timer.print("Poly1d, GbBufferV2: ");
	
	grbn::AddRels(gb3, gb, FnGetDeg{ gen_degs_t }, t_max);
	timer.print("GbWithCache, GbBuffer: ");

	grbn::AddRels<grbn::GbWithCache, FnGetDeg, grbn::GbBufferV2>(gb4, gb, FnGetDeg{ gen_degs_t }, t_max);
	timer.print("GbWithCache, GbBufferV2: ");

	std::cout << gb1.size() << ' ' << gb2.size() << ' ' << gb3.size() << ' ' << gb4.size() << '\n';

	return 0;
}

int main()
{
	return BenchmarkAddRels();
}
