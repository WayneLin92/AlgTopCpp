/* m1 is a sparse vector while [p1, p2) is a dense vector starting with index `index` */
bool divides(const Mon& m, array::const_iterator p1, array::const_iterator p2, size_t index)
{
	for (MonInd p = m.begin(); p != m.end(); ++p)
		if (p->exp > * (p1 + (size_t(p->gen) - index)))
			return false;
	return true;
}

/* return a basis in degree `deg` */
Mon1d get_basis(const Mon2d& leadings, const std::vector<Deg>& gen_degs, Deg deg)
{
	Mon1d result;
	array mon_dense;
	mon_dense.resize(gen_degs.size());
	size_t index = mon_dense.size() - 1;

	while (true) {
		//std::cout << "index=" << index << '\n';
		mon_dense[index] = T_MAX;
		int e_max = T_MAX;
		for (const Mon& lead : leadings[index]) {
			if (divides(lead, mon_dense.begin() + index, mon_dense.end(), index)) {
				if (lead[0].exp - 1 < e_max)
					e_max = lead[0].exp - 1;
			}
		}
		int qs = gen_degs[index].s ? deg.s / gen_degs[index].s : T_MAX;
		int qt = gen_degs[index].t ? deg.t / gen_degs[index].t : T_MAX;
		int qv = gen_degs[index].v ? deg.v / gen_degs[index].v : T_MAX;
		e_max = std::min({ e_max, qs, qt, qv });
		//std::cout << "e_max=" << e_max << '\n';
		mon_dense[index] = e_max;
		deg -= gen_degs[index] * e_max;

		bool move_right = false;
		if (deg != Deg{ 0, 0, 0 }) {
			if (index > 0)
				index--;
			else
				move_right = true;
		}
		else {
			Mon mon;
			for (size_t i = index; i < mon_dense.size(); ++i)
				if (mon_dense[i])
					mon.emplace_back(int(i), mon_dense[i]);
			result.push_back(std::move(mon));
			if (index > 0) {
				deg += gen_degs[index];
				mon_dense[index--]--;
			}
			else
				move_right = true;
		}
		if (move_right) {
			for (deg += gen_degs[index] * mon_dense[index], ++index; index < mon_dense.size() && mon_dense[index] == 0; ++index);
			if (index == mon_dense.size()) {
				//std::cout << "mon_dense=" << mon_dense << '\n';
				break;
			}
			else {
				deg += gen_degs[index];
				mon_dense[index--]--;
			}
		}
	}
	return result;
}

/* build the basis for t<=t_max */
void get_basis(const Mon2d& leadings, const std::vector<Deg>& gen_degs, std::map<Deg, Mon1d>& basis, int t_max)
{
	int t_min;
	if (basis.empty()) {
		basis[Deg{ 0, 0, 0 }].push_back({}); /* if no monomial present insert the unit */
		t_min = 1;
	}
	else
		t_min = basis.rbegin()->first.t + 1;

	for (int t = t_min; t <= t_max; ++t) {
		std::map<Deg, Mon1d> basis_new;
		for (int gen_id = (int)gen_degs.size() - 1; gen_id >= 0; --gen_id) { /* 58 is the gen_id of b1 */
			int t1 = t - gen_degs[gen_id].t;
			if (t1 >= 0) {
				auto p1 = basis.lower_bound(Deg{ 0, t1, 0 });
				auto p2 = basis.lower_bound(Deg{ 0, t1 + 1, 0 });
				for (auto p = p1; p != p2; ++p) {
					for (auto p_m = p->second.begin(); p_m != p->second.end(); ++p_m) {
						if (p_m->empty() || gen_id <= p_m->front().gen) {
							Mon mon(mul(*p_m, { {gen_id, 1} }));
							if ((size_t)gen_id >= leadings.size() || std::none_of(leadings[gen_id].begin(), leadings[gen_id].end(),
								[&mon](const Mon& _m) { return divides(_m, mon); }))
								basis_new[p->first + gen_degs[gen_id]].push_back(std::move(mon));
						}
					}
				}
			}
		}
		basis.merge(basis_new);
	}
}

inline void hash_combine(std::size_t& seed, int v) { seed ^= std::hash<int>{}(v)+0x9e3779b9 + (seed << 6) + (seed >> 2); }
inline size_t hash_range(const array& a)
{
	size_t seed = 0;
	for (int i : a)
		hash_combine(seed, i);
	return seed;
}