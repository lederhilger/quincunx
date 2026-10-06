module peirce;

import ellipticity : Ellipticity;
import jacobi : sncndn;
import riemann : π, τ, Point, Angles, Pole, napkin, normalizeLongitude;
import std.complex : Complex, complex, complexAbs = abs, complexSqrt = sqrt;
import std.math : atan, atan2, sqrt, hypot, copysign, fabs;

private enum double modulus = 0.70710678118654752440;
private enum double sqrt2 = 1.41421356237309504880;
private enum double lemniscate = 0.64359425290558262474;
private enum double fudge = .02;
private enum double radius = lemniscate + fudge;

struct Peirce
{
	private Ellipticity ellipticity = Ellipticity(modulus);
	private Pole centered = Pole.south;
	private double meridian = 0;

	this(Pole center, double centralMeridian = 0) @safe pure
	{
		centered = center;
		meridian = normalizeLongitude(centralMeridian);
		ellipticity = Ellipticity(modulus);
	}

	@property Pole center() const @safe @nogc pure nothrow
	{
		return centered;
	}
	@property centralMeridian() const @safe @nogc pure nothrow
	{
		return meridian;
	}
	@property double K() const @safe pure
	{
		return ellipticity.K;
	}

	Point integral(double x, double z) const @safe pure
	{
		const exchange = fabs(z) > fabs(x);
		const ζ = exchange ? complex(fabs(z), fabs(x)) : complex(fabs(x), fabs(z));
		auto w = integrateOctant(ζ);
		if (exchange) {w = complex(w.im, w.re);}
		return Point(copysign(w.re, x), copysign(w.im, z));
	}

	private Complex!double integrateOctant(Complex!double ζ) const @safe pure
	{
		if (complexAbs(ζ) <= radius) {return series(ζ);}
		const ξ = complexSqrt(((1 - ζ) * (1 + ζ)) / (1 + ζ * ζ));
		if (complexAbs(ξ) <= radius)
		{
			return complex(K / sqrt2, 0.0) - series(ξ);
		}

		const A = hypot(ζ.re + 1, ζ.im);
		const B = hypot(ζ.re - 1, ζ.im);
		const α = .5 * A + .5 * B;
		const β = ζ.re / α;
		const numerator = (1 - β) * (1 + β);
		const Tϕ = numerator / (β * β);
		const q = ζ.im * ζ.im / (β * β * numerator);
		const d = 2 + Tϕ - q;
		const radical = sqrt(d * d + 8 * q);
		const Tμ = d >= 0 ? 4 * q / (radical + d) : .5 * (radical - d);
		const u = 1 + .5 * Tμ;
		const λ = atan(sqrt(Tϕ / u));
		const μ = atan(sqrt(Tμ));
		return complex((K - ellipticity.F(λ)) / sqrt2, ellipticity.F(μ) / sqrt2);
	}

	Point diamond(double x, double z) const @safe pure
	{
		const I = integral(x, z);
		return Point(-sqrt2 * I.x / K, -sqrt2 * I.z / K);
	}

	Point direct(double θ, double ϕ) const @safe pure
	{
		ϕ = centered == Pole.south ? ϕ : -ϕ;
		if (ϕ == -π / 2) {return Point(0,0);}
		if (ϕ == π / 2) {return Point(1, 1);}
		θ = normalizeLongitude(θ - meridian);
		const ζ = napkin(θ, ϕ);
		auto h = diamond(ζ.x, ζ.z);
		if (ϕ == 0)
		{
			if (θ == 0) {return Point(-1, 0);}
			if (θ == π / 2) {return Point(0, -1);}
			if (θ == π) {return Point(1, 0);}
			if (θ == 3 * π / 2) {return Point(0, 1);}
			const norm = fabs(h.x) + fabs(h.z);
			h.x /= norm;
			h.z /= norm;
		}
		if (ϕ <= 0) {return h;}

		const sx = θ < π / 2 || θ >= 3  *π / 2 ? -1.0 : 1.0;
		const sz = θ < π ? -1.0 : 1.0;
		return Point(sx * (1 - fabs(h.z)), sz * (1 - fabs(h.x)));
	}

	Angles inverse(double x, double z) const @safe pure
	{
		const outer = fabs(x) + fabs(z) > 1;
		double α = x;
		double β = z;
		if (outer)
		{
			α = copysign(1 - fabs(z), x);
			β = copysign(1 - fabs(x), z);
		}
		if (α == 0 && β == 0) {return sphere(meridian, outer ? π / 2 : -π / 2);}
		const u = sncndn(ellipticity, K * α);
		const v = sncndn(ellipticity, K * β);
		const snx = u.cn / u.dn;
		const cnx = -ellipticity.kPrime * u.sn / u.dn;
		const dnx = ellipticity.kPrime / u.dn;
		const denominator = 1 - dnx * dnx * v.sn * v.sn;
		const re = cnx * v.cn / denominator;
		const im = -snx * v.sn * dnx * v.dn / denominator;
		const r = hypot(re, im);
		double magnitude = 2 * atan((1 - r) / (1 + r));
		if (magnitude < 0) {magnitude = 0;}
		const ϕ = outer ? magnitude : -magnitude;
		const θ = normalizeLongitude(atan2(im, re));
		return sphere(normalizeLongitude(θ + meridian), ϕ);
	}

	private Angles sphere(double θ, double ϕ) const @safe @nogc pure nothrow
	{
		return Angles(θ, centered == Pole.south ? ϕ : -ϕ);
	}
}

private Complex!double series(Complex!double ζ) @safe pure
{
	const squared = ζ * ζ;
	const bisquared = squared * squared;
	const q = complexAbs(bisquared);
	auto product = ζ;
	auto sum = ζ;
	auto compensation = complex(0.0, 0.0);
	if (q == 0) {return sum;}
	foreach (n; 1 .. 128)
	{
		product *= bisquared * ((2.0 * n - 1) / (2.0 * n));
		const term = product / (4.0 * n + 1);
		const adjusted = term - compensation;
		const next = sum + adjusted;
		compensation = (next - sum) - adjusted;
		sum = next;
		if (complexAbs(term * q / (1 - q)) <= .125 * double.epsilon * complexAbs(sum)) {return sum;}
	}
	return sum;
}

Point peirce(double θ, double ϕ, Pole center = Pole.south, double centralMeridian = 0) @safe pure
{
	const projection = Peirce(center, centralMeridian);
	return projection.direct(θ, ϕ);
}