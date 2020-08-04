#include "CppUnitTest.h"
#include "../algtop/database.h"
#include "../algtop/database.cpp"
#include "../algtop/algmod2.cpp"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestAlgTop
{
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
		};
	};
}
