module riemann;
import std.math : sin, cos, tan, atan, atan2, PI, fabs, hypot;

enum double π = PI;
enum double τ = 2 * π;

struct Point
{
	double x;
	double z;
}

struct Angles
{
	double θ;
	double ϕ;
}

enum Pole {south, north}

double normalizeLongitude(double θ) @safe pure
{
	θ %= τ;
	if (θ < 0) {θ += τ;}
	return θ == 0 || θ == τ ? 0 : θ;
}

Point polar(double r, double θ) @safe pure
{
	θ = normalizeLongitude(θ);
	if (θ == 0) {return Point(r, 0);}
	if (θ == π / 2) {return Point(0, r);}
	if (θ == π) {return Point(-r, 0);}
	return Point(r * cos(θ), r * sin(θ));
}

Point stereograph(double θ, double ϕ) @safe pure
{
	if (ϕ == -π / 2) {return Point(0,0);}
	if (ϕ == π / 2) {return Point(double.infinity, 0);}
	const r = ϕ == 0 ? 1.0 : tan(.5 * (π / 2 + ϕ));
	return polar(r, θ);
}

Point napkin(double θ, double ϕ) @safe pure
{
	const r = ϕ == 0 ? 1.0 : tan(.5 * (π / 2 - fabs(ϕ)));
	return polar(r, θ);
}

Angles grapheostere(double x, double z) @safe pure
{
	const r = hypot(x, z);
	if (r == 0) {return Angles(0, -π / 2);}
	return Angles(normalizeLongitude(atan2(z, x)), 2 * atan(r) - π / 2);
}