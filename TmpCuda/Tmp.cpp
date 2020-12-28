#include "groebner.h"
#include "benchmark.h"
#include "database.h"
#include "myio.h"
#include <iostream>

int decomposables()
{
	Timer timer;
	Database db(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\E4t.db)");
	std::vector<Deg> gen_degs = db.load_gen_degs("HA3_generators");
	Poly1d gb = db.load_gb("HA3_relations", 189);
	int count = 0;
	for (int i = 0; i < (int)gen_degs.size(); ++i) {
		if (grbn::Reduce({ {{i, 1}} }, gb).empty())
			++count;
	}
	std::cout << "count=" << count << '\n';
	return 0;
}

int main()
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