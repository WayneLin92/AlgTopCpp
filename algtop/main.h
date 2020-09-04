#ifndef MAIN_H
#define MAIN_H

#include "database.h"
#include <iostream>
#include <vector>
#include <map>

/*--------- generate_basis.cpp ---------*/

/* Execute simple commands.*/
void execute_cmd(sqlite3* conn, const std::string& cmd);
int main_generate_basis(int argc, char** argv);
void load_gen_degs(sqlite3* conn, const std::string& table_name, std::vector<Deg>& gen_degs);
void load_basis(sqlite3* conn, const std::string& table_name, std::map<Deg, array2d>& basis);
void load_leading_terms(sqlite3* conn, const std::string& table_name, array3d& leadings);

/*--------- generate_mon_diffs.cpp ---------*/

int main_generate_diff(int argc, char** argv);
void load_gen_diffs(sqlite3* conn, const std::string& table_name, array3d& diffs);
void load_gb(sqlite3* conn, const std::string& table_name, array3d& gb);

/*--------- generate_ss.cpp ---------*/

int main_generate_ss(int argc, char** argv);

/* generate_next_page.cpp */

int main_generate_next_page(int argc, char** argv);

/*--------- generate_E4t.cpp ---------*/

int main_generate_E4t(int argc, char** argv);
int main_test(int argc, char** argv);


#endif
