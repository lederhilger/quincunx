module ellipticity;

import std.math : sin, asin, cos, atan2, fabs, sqrt, PI;
import agm : AGM;

struct Ellipticity
{
	double k;
	double kPrime;
	double ε;

	private AGM agm;
	private double[16] ratio;
	private size_t N;
	private double scale;

	this (double modulus, double ε = 1e-16) @safe @nogc pure nothrow
	{
		k = modulus;
		assert(k >= 0 && k < 1, "Degenerate.");
		kPrime = sqrt((1 + k) * (1 - k));
		this.ε = ε;

		agm = AGM(1.0, kPrime, k, this.ε);
		N = agm.length - 1;

		foreach (i; 0 .. agm.length)
		{
			ratio[i] = agm.C[i] / agm.A[i];
		}
		scale = cast(double)(1UL << N) * agm.A[N];
	}

	@property double K() const @safe @nogc pure nothrow
	{
		return .5 * PI / agm.mean;
	}

	void ϕ(double x, out double ϕ0, out double ϕ1) const @safe @nogc pure nothrow
	{
		pragma(inline, true);
		
		double Φ = scale * x;
		double Ψ = Φ;

		foreach (n; 0 .. N)
		{
			immutable size_t index = N - n;
			immutable double next = .5 * (asin(ratio[index] * sin(Φ)) + Φ);
			Ψ = Φ;
			Φ = next;
		}
		ϕ0 = Φ;
		ϕ1 = Ψ;
	}

	double amplitude(double x) const @safe @nogc pure nothrow
	{
		double ϕ0, ϕ1;
		ϕ(x, ϕ0, ϕ1);
		return ϕ0;
	}

	double F(double φ) const @safe @nogc pure nothrow
	{
		enum double π = PI;
		assert(φ >= -π && φ <= π);
		if (k == 0 || φ == 0) {return φ;}
		immutable double sign = φ < 0 ? -1.0 : 1.0;
		double ϕ = fabs(φ);
		immutable double complete = K;
		if (ϕ == π) {return sign * (2 * complete);}
		if (ϕ == π / 2) {return sign * complete;}
		immutable bool reflect = ϕ > π / 2;
		if (reflect) {ϕ = π - ϕ;}
		size_t rotation = 0;
		foreach (n; 0 .. N)
		{
			ϕ += atan2(agm.B[n] * sin(ϕ), agm.A[n] * cos(ϕ));
			rotation *= 2;
			if (ϕ >= π)
			{
				ϕ -= π;
				++ rotation;
			}
		}
		immutable double integral = (ϕ + rotation * π) / scale;
		return sign * (reflect ? 2 * complete - integral : integral);
	}
}