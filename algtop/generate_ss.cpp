#include "database.h"

struct BasisSS
{
	array base_indices;
	array diff_indices;
	int level;
	int s;
	int t;
	int v;

	BasisSS(array&& base_indices_, array&& diff_indices_, int level_, int s_, int t_, int v_) :
		base_indices(base_indices_), diff_indices(diff_indices_), level(level_), s(s_), t(t_), v(v_) {};
	BasisSS() : level(-1), s(-1), t(-1), v(-1) {};
};

struct DgaBasisMon
{
	array diff_indices;
	int s;
	int t;
	int v;

	DgaBasisMon(array&& diff_indices_, int s_, int t_, int v_) :
		diff_indices(diff_indices_), s(s_), t(t_), v(v_) {};
};

void load_basis(sqlite3* conn, const std::string& table_name, std::vector<DgaBasisMon>& basis)
{
	sqlite3_stmt* stmt;
	std::string cmd = std::string("SELECT diff, s, t, v FROM ") + table_name + " ;";
	sqlite3_prepare_v100(conn, cmd, &stmt);
	while (sqlite3_step(stmt) == SQLITE_ROW) {
		basis.emplace_back(str_to_array(sqlite3_column_str(stmt, 0)),
			sqlite3_column_int(stmt, 1), sqlite3_column_int(stmt, 2), sqlite3_column_int(stmt, 3));
	}
	sqlite3_finalize(stmt);
}

bool cmp_ss(const BasisSS& b, const std::array<int, 3>& value)
{
	if (b.t < value[1])
		return true;
	else if (b.t == value[1])
		if (b.s < value[0])
			return true;
		else if (b.s == value[0])
			if (b.v < value[2])
				return true;
	return false;
}

void generate_ss(sqlite3* conn, const std::string& table_name_basis, const char* table_ss, int r)
{
	/* Compile SQL statements */
	sqlite3_stmt* stmt_update_basis;
	std::string cmd_update_basis = std::string("INSERT INTO ") + table_ss + " (base, diff, level, s, t, v) VALUES (?1, ?2, ?3, ?4, ?5, ?6);";
	sqlite3_prepare_v100(conn, cmd_update_basis, &stmt_update_basis);

	std::vector<std::array<int, 5>> indices_tsv;
	std::string table_indices_tsv = std::string(table_name_basis) + "_indices_stv";
	load_mon_indices_tsv(conn, table_indices_tsv, indices_tsv);
	std::cout << "num_mons_tsv loaded! Size=" << indices_tsv.size() << '\n';

	std::vector<DgaBasisMon> basis;
	load_basis(conn, table_name_basis, basis);
	std::cout << "basis loaded! Size=" << basis.size() << '\n';

	std::vector<BasisSS> basis_ss;
	basis_ss.resize(basis.size());

	/* fill basis_ss */
	int prev_t = 0;
	for (const auto& index : indices_tsv) {
		const int &t = index[0], &s = index[1], &v = index[2], &a = index[3], &b = index[4];
		if (t > prev_t) {
			prev_t = t;
			std::cout << "t=" << t << '\n';
		}

		array range;
		for (int i = 0; i < b - a; i++)
			range.push_back(i);

		int i;
		array lead_image;
		for (i = a; i < b; i++) {
			if (basis_ss[i].base_indices.empty())
				break;
			lead_image.push_back(basis_ss[i].base_indices[0]);
		}
		std::sort(lead_image.begin(), lead_image.end());

		array x;
		std::set_symmetric_difference(range.begin(), range.end(), lead_image.begin(), lead_image.end(), std::back_inserter(x));

		array2d fx;
		for (int xi : x)
			fx.push_back(basis[size_t(a) + xi].diff_indices);

		array2d image, kernel;
		get_image_kernel(x, fx, image, kernel);

		/* fill with the kernel after the image */
		array lead_kernel;
		for (size_t j = 0; j < kernel.size(); j++) {
			lead_kernel.push_back(kernel[j][0]);
			basis_ss[i + j] = BasisSS(std::move(kernel[j]), array(), N - r, s, t, v);
		}
		std::sort(lead_kernel.begin(), lead_kernel.end());

		/* fill with the image */
		auto ab_image = get_ab(indices_tsv, s + 1, t, v - r);
		for (size_t k = 0; k < image.size(); k++)
			basis_ss[ab_image.first + k] = BasisSS(std::move(image[k]), array(), r, s + 1, t, v - r);

		/* fill with the rest */
		array rest;
		std::set_symmetric_difference(x.begin(), x.end(), lead_kernel.begin(), lead_kernel.end(), std::back_inserter(rest));
		for (size_t j = 0; j < rest.size(); j++) {
			basis_ss[i + kernel.size() + j] = BasisSS({ rest[j] }, std::move(basis[size_t(a) + rest[j]].diff_indices), N, s, t, v);
		}
	}

	/* Insert into the database */
	execute_cmd(conn, "BEGIN TRANSACTION");
	prev_t = 0;
	for (const auto& b_ss : basis_ss) {
		if (b_ss.t > prev_t) {
			prev_t = b_ss.t;
			std::cout << "t=" << b_ss.t << '\n';
			execute_cmd(conn, "END TRANSACTION");
			execute_cmd(conn, "BEGIN TRANSACTION");
		}
		sqlite3_bind_str(stmt_update_basis, 1, array_to_str(b_ss.base_indices));
		sqlite3_bind_str(stmt_update_basis, 2, array_to_str(b_ss.diff_indices));
		sqlite3_bind_int(stmt_update_basis, 3, b_ss.level);
		sqlite3_bind_int(stmt_update_basis, 4, b_ss.s);
		sqlite3_bind_int(stmt_update_basis, 5, b_ss.t);
		sqlite3_bind_int(stmt_update_basis, 6, b_ss.v);
		sqlite3_step(stmt_update_basis);
		sqlite3_reset(stmt_update_basis);
		//std::cout << b_ss.base_indices << " " << b_ss.diff_indices << " " << b_ss.level << " " << b_ss.s << " " << b_ss.t << " " << b_ss.v << '\n';
	}
	execute_cmd(conn, "END TRANSACTION");

	sqlite3_finalize(stmt_update_basis);
}

int main_test()
{
	array2d spaceV = { {1, 2, 3, 4}, {2, 3, 4}, {3} };
	simplify_space(spaceV);
	std::cout << spaceV << '\n';
	return 0;
}

int main_generate_ss(int argc, char** argv)
{
	//return main_test();

	sqlite3* conn;
	sqlite3_open(R"(C:\Users\lwnpk\Documents\MyProgramData\Math_AlgTop\database\ss.db)", &conn);

	const char* table_basis, * table_ss;
	int t_max;
	if (argc == 1) {
		table_basis = "E2_basis";
		table_ss = "E2_ss";
		t_max = 50;
	}
	else {
		table_basis = argv[1];
		table_ss = argv[2];
		t_max = atoi(argv[3]);
	}

	generate_ss(conn, table_basis, table_ss, 2);

	sqlite3_close(conn);
	return 0;
}