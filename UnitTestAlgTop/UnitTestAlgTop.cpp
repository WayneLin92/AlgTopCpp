#include "CppUnitTest.h"
#include "../algtop/database.h"
#include "../algtop/database.cpp"
#include "../algtop/algmod2.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestAlgTop
{
	TEST_CLASS(AlgebraMod2)
	{
	public:

		TEST_METHOD(test_pow)
		{
			array2d poly = { {1, 1}, {0, 1} };
			array3d gb = { {{0, 1, 1, 2}, {-1, 3}} };
			array2d power = pow(poly, 3, gb);
			array2d anwer = { {1, 3}, {0, 2, 1, 1}, {0, 3}, {-1, 3} };
			Assert::AreEqual(array2d_to_str(anwer), array2d_to_str(power));
		};

		TEST_METHOD(test_evaluate)
		{
			array2d p = { {1, 3}, {0, 3}, {0, 1, 2, 1} };
			array3d images = { {{1, 1}, {0, 1}}, {}, {{3, 1}} };
			array2d fp = evaluate(p, [&images](int gen_id) {return images[gen_id]; }, {});
			array2d anwer = { {1, 1, 3, 1}, {1, 3}, {0, 1, 3, 1}, {0, 1, 1, 2}, {0, 2, 1, 1}, {0, 3} };
			Assert::AreEqual(array2d_to_str(anwer), array2d_to_str(fp));
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

		TEST_METHOD(test_get_image_kernel)
		{
			array x = { 1, 2, 3 };
			array2d fx = { {1, 3, 4}, {2, 4, 5}, {1, 2, 3, 5} };
			array2d image, kernel;
			get_image_kernel(x, fx, image, kernel);
			Assert::AreEqual(array2d_to_str(array2d({ {1, 3, 4}, {2, 4, 5} })), array2d_to_str(image));
			Assert::AreEqual(array2d_to_str(array2d({ {1, 2, 3} })), array2d_to_str(kernel));
			
			image.clear(); kernel.clear();
			get_image_kernel(fx, image, kernel);
			Assert::AreEqual(array2d_to_str(array2d({ {0, 1, 2} })), array2d_to_str(kernel));
		};

		TEST_METHOD(test_quotient)
		{
			array2d spaceV = { {1, 2, 3, 4}, {2, 3, 4}, {3, 4}, {4} };
			array2d spaceW = { {2, 3}, {4} };
			array2d quotient = quotient_space(spaceV, spaceW);
			array2d answer = { {1}, {3} };
			Assert::AreEqual(array2d_to_str(answer), array2d_to_str(quotient));
#ifdef _DEBUG
			spaceW = { {2, 3}, {5} };
			bool bException = false;
			try {
				quotient = quotient_space(spaceV, spaceW);
			}
			catch (const char* e) {
				Assert::AreEqual("cec7f701-0911-482a-a63f-1caaa646591b", e);
				bException = true;
			}
			Assert::AreEqual(true, bException);
#endif
		};

		TEST_METHOD(test_simplify_space)
		{
			array2d spaceV = { {1, 2, 3, 4}, {2, 3, 4}, {3} };
			array2d answer = { {1}, {2, 4}, {3} };
			Assert::AreEqual(array2d_to_str(answer), array2d_to_str(simplify_space(spaceV)));
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
			array3d rels;
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
						rel += index1 < index2 ? Mon({ index1, 1, index2, 1 }) : Mon({ index2, 1, index1, 1 });
					}
					rels.push_back(std::move(rel.data));
				}
			}

			array3d gb;
			add_rels(gb, rels, gen_degs, -1);
			int answer = 65;
			Assert::AreEqual(answer, int(gb.size()));
		}

		TEST_METHOD(test_groebner)
		{
			array gen_degs = { 1, 1, 1, 1 };
			array3d rels = { {{2, 3}, {1, 3}, {0, 1, 1, 1, 2, 1}, {0, 3}} };
			array3d rels1 = { {{2, 1}, {1, 1}, {0, 1}} };

			array3d gb;
			add_rels(gb, rels, gen_degs, -1);
			add_rels(gb, rels1, gen_degs, -1);
			Assert::AreEqual(1, int(gb.size()));
		}
	};
}
