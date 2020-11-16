#include "CppUnitTest.h"
#include "../algtop/database.h"
#include "../algtop/myparser.h"
#include "../algtop/utilities.cpp"
#include "../algtop/myparser.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

std::string ToString(array2d a) {
	std::stringstream ss;
	ss << a;
	return ss.str();
}

namespace UnitTestAlgTop
{
	TEST_CLASS(AlgebraMod2)
	{
	public:

		TEST_METHOD(test_pow)
		{
			Poly poly = { {{1, 1}}, {{0, 1}} };
			Poly1d gb = { {{{0, 1}, {1, 2}}, {{-1, 3}}} };
			Poly power = pow(poly, 3, gb);
			Poly anwer = { {{1, 3}}, {{0, 2}, {1, 1}}, {{0, 3}}, {{-1, 3}} };
			Assert::AreEqual(Poly_to_str(anwer), Poly_to_str(power));
		};

		TEST_METHOD(test_evaluate)
		{
			Poly p = { {{1, 3}}, {{0, 3}}, {{0, 1}, {2, 1}} };
			Poly1d images = { {{{1, 1}}, {{0, 1}}}, {}, {{{3, 1}}} };
			Poly fp = evaluate(p, [&images](int gen_id) {return images[gen_id]; }, {});
			Poly anwer = { {{1, 1}, {3, 1}}, {{1, 3}}, {{0, 1}, {3, 1}}, {{0, 1}, {1, 2}}, {{0, 2}, {1, 1}}, {{0, 3}} };
			Assert::AreEqual(Poly_to_str(anwer), Poly_to_str(fp));
		};
	};

	TEST_CLASS(LinearAlgebra)
	{
	public:
		
		TEST_METHOD(test_residue)
		{
			array2d spaceV = { {1, 3, 4}, {2, 4, 5} };
			array v = { 1, 2 };
			array r = residue(spaceV, v);
			array r1 = residue(spaceV, std::move(v));
			array expected = { 3, 5 };
			Assert::AreEqual(array_to_str(expected), array_to_str(r));
			Assert::AreEqual(array_to_str(expected), array_to_str(r1));
		};

		TEST_METHOD(test_set_linear_map)
		{
			array x = { 1, 2, 3 };
			array2d fx = { {1, 3, 4}, {2, 4, 5}, {1, 2, 3, 5} };
			array2d image, kernel, g;
			set_linear_map(x, fx, image, kernel, g);
			Assert::AreEqual(ToString(array2d({ {1, 3, 4}, {2, 4, 5} })), ToString(image));
			Assert::AreEqual(ToString(array2d({ {1, 2, 3} })), ToString(kernel));

			image.clear(); kernel.clear(); g.clear();
			set_linear_map(fx, image, kernel, g);
			Assert::AreEqual(ToString(array2d({ {0, 1, 2} })), ToString(kernel));
		};
		TEST_METHOD(test_get_image)
		{
			array2d spaceV = { {1, 2, 3, 4}, {2, 3, 4}, {3, 4}, {4} };
			array2d f = { {1}, {2, 5}, {3}, {4} };
			array v = { 2, 3 };
			array answer = { 2, 4, 5 };
			Assert::AreEqual(array_to_str(answer), array_to_str(get_image(spaceV, f, v)));
#ifdef _DEBUG
			v = { 2, 5 };
			bool bException = false;
			try {
				get_image(spaceV, f, v);
			}
			catch (const char* e) {
				Assert::AreEqual("6a4fe8a1", e);
				bException = true;
			}
			Assert::AreEqual(true, bException);
#endif
		};

		TEST_METHOD(test_quotient)
		{
			array2d spaceV = { {1, 2, 3, 4}, {2, 3, 4}, {3, 4}, {4} };
			array2d spaceW = { {2, 3}, {4} };
			array2d quotient = quotient_space(spaceV, spaceW);
			array2d answer = { {1}, {3} };
			Assert::AreEqual(ToString(answer), ToString(quotient));
#ifdef _DEBUG
			spaceW = { {2, 3}, {5} };
			bool bException = false;
			try {
				quotient = quotient_space(spaceV, spaceW);
			}
			catch (const char* e) {
				Assert::AreEqual("cec7f701", e);
				bException = true;
			}
			Assert::AreEqual(true, bException);
#endif
		};

		TEST_METHOD(test_simplify_space)
		{
			array2d spaceV = { {1, 2, 3, 4}, {2, 3, 4}, {3} };
			array2d answer = { {1}, {2, 4}, {3} };
			Assert::AreEqual(ToString(answer), ToString(simplify_space(spaceV)));
		}
	};


	TEST_CLASS(Groebner)
	{
	public:

		TEST_METHOD(test_groebner_B7)
		{
			array gen_degs;
			int n_max = 7;
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

			Poly1d gb;
			add_rels(gb, rels, gen_degs, -1);
			size_t gb_size = gb.size();
			size_t answer = 65;
			Assert::AreEqual(answer, gb_size);
		}

		TEST_METHOD(test_groebner)
		{
			array gen_degs = { 1, 1, 1, 1 };
			Poly1d rels = { {{{2, 3}}, {{1, 3}}, {{0, 1}, {1, 1}, {2, 1}}, {{0, 3}}} };
			Poly1d rels1 = { {{{2, 1}}, {{1, 1}}, {{0, 1}}} };

			Poly1d gb;
			add_rels(gb, rels, gen_degs, -1);
			add_rels(gb, rels1, gen_degs, -1);
			Assert::AreEqual(1, int(gb.size()));
		}

		TEST_METHOD(test_ann_seq)
		{
			array gen_degs = { 1, 1, 1, 1 }; /* x, y, z, w */
			Poly rel0 = { {{0, 4}} }; /* x^4 */
			Poly rel1 = { {{0, 2}, {1, 4}} }; /* x^2y^4 */
			Poly rel2 = { {{0, 1}, {2, 2}}, {{0, 2}, {1, 1}} }; /* x^2y + xz^2 */
			Poly rel3 = { {{3, 4}}, {{2, 4}} }; /* z^2 + w^4 */
			Poly1d gb;
			add_rels(gb, { rel0, rel1, rel2, rel3 }, gen_degs, -1);

			Poly2d ann = ann_seq(gb, { Poly{ {{0, 1}} } }, gen_degs, -1);
			Assert::AreEqual(3, int(ann.size()));
		}
	};
}
