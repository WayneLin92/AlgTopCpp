#ifndef MAIN_H
#define MAIN_H

#include "database.h"

void generate_basis(const Database& db, const std::string& table_prefix, int t_max, bool drop_existing = false);
void generate_mon_diffs(const Database& db, const std::string& table_prefix, int r);
void generate_ss(const Database& db, const std::string& table_name_basis, const std::string& table_ss, int r);
void generate_next_page(const Database& db, const std::string& table_prefix, const std::string& table_H_prefix, int r);
int main_generate_E4t(int argc, char** argv);
int main_generate_E4t_1(int argc, char** argv);
int main_test(int argc, char** argv);


#endif
