module agm;
import std.math : sqrt;

struct AGM
{
	double[16] A;
	double[16] B;
	double[16] C;
	size_t length;

	this (double a, double b, double c, double ε = 1e-16) @safe @nogc pure nothrow
	{
		length = 0;
		do
		{
			assert(length < A.length);
			A[length] = a;
			B[length] = b;
			C[length] = c;
			++length;
			immutable double α = .5 * (a + b);
			immutable double β = sqrt(a * b);
			immutable double γ = .5 * (a - b);
			a = α;
			b = β;
			c = γ;
		}
		while (c > ε);
		assert(length < A.length);
		A[length] = a;
		B[length] = b;
		C[length] = c;
		++length;
	}

	@property double mean() const @safe @nogc pure nothrow
	{
		return A[length - 1];
	}
}