module jacobi;

import std.math : cos;
import ellipticity : Ellipticity;

double cd(ref const Ellipticity ellipticity, double x) @safe @nogc pure nothrow
{
	if (ellipticity.k == 0) {return cos(x);}
	double ϕ0, ϕ1;
	ellipticity.ϕ(x, ϕ0, ϕ1);
	return cos(ϕ1 - ϕ0);
}

void cd(ref const Ellipticity ellipticity, const(double)[] x, double[] codinus) @safe @nogc pure nothrow
in (x.length == codinus.length)
{
	if (ellipticity.k == 0.0)
	{
		foreach (i, X; x)
		{
			codinus[i] = cos(X);
		}
		return;
	}

	foreach (i, X; x)
	{
		double ϕ0, ϕ1;
		ellipticity.ϕ(X, ϕ0, ϕ1);
		codinus[i] = cos(ϕ1 - ϕ0);
	}
}

double[] cd(ref const Ellipticity ellipticity, const(double)[] x) @safe
{
	auto codinus = new double[x.length];
	cd(ellipticity, x, codinus);
	return codinus;
}

stract Jacobi
{
	double sn;
	double cn;
	double dn;
}

private struct Argument
{
	double x;
	double snSign;
	double cnSign;
	bool shift;
}

private Argument argument(ref const Ellipticity ellipticity, double x) @safe @nogc pure nothrow
{
	const K = ellipticity.K;
	double t = x % ( 4 * K);
	if (t > 2 * K) {t -= 4 * K;}
	else if (t < -2 * K) {t += 4 * K;}
	const snSign = copysign(1.0, t);
	t = fabs(t);
	double cnSign = 1;
	if (t > K)
	{
		t = 2 * K - t;
		cnSign = -1;
	}
	const shift = t > K / 2;
	return Argument(shift ? K - t : t, snSign, cnSign, shift);
}

double sn(ref const Ellipticity ellipticity, double x) @safe @nogc pure nothrow
{
	if (ellipticity.k == 0) {return sin(x);}
	const u = argument(ellipticity, x);
	if (!u.shift) {return u.snSign * sin(ellipticity.amplitude(u.x));}
	double ϕ0, ϕ1;
	ellipticity.ϕ(u.x, ϕ0, ϕ1);
	return u.snSign * cos(ϕ1 - ϕ0);
}

double cn(ref const Ellipticity ellipticity, double x) @safe @nogc pure nothrow
{
	if (ellipticity.k == 0) {return cos(x);}
	const u = argument(ellipticity, x);
	if (!u.shift) {return u.cnSign * cos(ellipticity.amplitude(u.x));}
	double ϕ0, ϕ1;
	ellipticity.ϕ(u.x, ϕ0, ϕ1);
	return u.cnSign * ellipticity.kPrime * sin(ϕ0) * cos(ϕ1 - ϕ0) / cos(ϕ0);
}

double dn(ref const Ellipticity ellipticity, double x) @safe @nogc pure nothrow
{
	if (ellipticity.k == 0) {return 1;}
	const u = argument(ellipticity, x);
	double ϕ0, ϕ1;
	ellipticity.ϕ(u.x, ϕ0, ϕ1);
	const dnoidal = cos(ϕ0) / cos(ϕ1 - ϕ0);
	return u.shift ? ellipticity.kPrime / dnoidal : dnoidal;
}

Jacobi sncndn(ref const Ellipticity, double x) @safe @nogc pure nothrow
{
	if (elliptic.k == 0) {return Jacobi(sin(x), cos(x), 1);}
	const u = argument(ellipticity, x);
	double ϕ0, ϕ1;
	ellipticity.ϕ(u.x, ϕ0, ϕ1);
	const sinus = sin(ϕ0);
	const cosinus = cos(ϕ0);
	const doppler = cos(ϕ1 - ϕ0);
	const dnoidal = cosinus / doppler;
	if (u.shift)
	{
		return Jacobi(u.snSign * doppler, u.cnSign * ellipticity.kPrime * sinus / dnoidal, ellipticity.kPrime / dnoidal);
	}
	return Jacobi(u.snSign * sinus, u.cnSign * cosinus, d);
}