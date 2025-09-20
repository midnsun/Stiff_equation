#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <limits>
//#include <omp.h>

std::vector<double> operator+ (std::vector<double> l, std::vector<double> r) {
	std::vector<double> res(l.size());
	for (size_t i = 0; i < res.size(); ++i) res[i] = l[i] + r[i];
	return res;
}

std::vector<double> operator- (std::vector<double> l, std::vector<double> r) {
	std::vector<double> res(l.size());
	for (size_t i = 0; i < res.size(); ++i) res[i] = l[i] - r[i];
	return res;
}

std::vector<double> operator* (double a, std::vector<double> r) {
	std::vector<double> res(r.size());
	for (size_t i = 0; i < res.size(); ++i) res[i] = a * r[i];
	return res;
}

double norm2(const std::vector<double>& v) {
	double res = 0.0;
	for (size_t i = 0; i < v.size(); ++i) {
		res += v[i] * v[i];
	}
	return std::sqrt(res);
}

struct point {
	double x;
	std::vector<double> V;
	point() : x(0.0) {
		V.resize(2);
	}
	point(const point& p) : x(p.x), V(p.V) {}
	point(point&& p) : x(p.x), V(std::move(p.V)) {}
	point& operator= (const point& p) {
		x = p.x;
		V = p.V;
		return *this;
	}
	point& operator= (point&& p) {
		x = p.x;
		V = p.V;
		return *this;
	}
};

std::vector<double> rhs(double x, const std::vector<double>& V) {
	std::vector<double> Y(V.size());
	std::vector< std::vector<double> > M(V.size(), std::vector<double>(V.size()));

	M[0][0] = -500.005;		M[0][1] = 499.995;
	M[1][0] = 499.995;		M[1][1] = -500.005;

	for (size_t i = 0; i < M.size(); ++i) {
		for (size_t j = 0; j < M[i].size(); ++j) {
			Y[i] += M[i][j] * V[j];
		}
	}

	return Y;
}

std::vector<double> rhs(const point& P) {
	return rhs(P.x, P.V);
}

std::vector<double> answer(double x, const std::vector<double>& S) {
	std::vector<double> Y(S.size());
	double C2 = (S[0] + S[1]) / 2, C1 = (S[1] - S[0]) / 2;
	std::vector<double> v1(2), v2(2);

	v1[0] = -1.0;	v2[0] = 1.0;
	v1[1] = 1.0;	v2[1] = 1.0;

	Y = C1 * exp(-1000 * x) * v1 + C2 * exp(-0.01 * x) * v2;

	return Y;
}

bool isless_all(const point& r, const point& l) {
	if (r.x > l.x) return false;
	for (size_t i = 0; i < r.V.size(); ++i) if (r.V[i] > l.V[i]) return false;

	return true;
}

bool ismore_all(const point& r, const point& l) {
	if (r.x < l.x) return false;
	for (size_t i = 0; i < r.V.size(); ++i) if (r.V[i] < l.V[i]) return false;

	return true;
}

bool isless_any(const point& r, const point& l) {
	if (r.x < l.x) return true;
	for (size_t i = 0; i < r.V.size(); ++i) if (r.V[i] < l.V[i]) return true;

	return false;
}

bool ismore_any(const point& r, const point& l) {
	if (r.x > l.x) return true;
	for (size_t i = 0; i < r.V.size(); ++i) if (r.V[i] > l.V[i]) return true;

	return false;
}

point RK4(const point& P, double h) {
	point Res;
	Res.x = P.x + h;
	size_t n = Res.V.size();
	
	std::vector<double> k1(n);
	std::vector<double> k2(n);
	std::vector<double> k3(n);
	std::vector<double> k4(n);

	k1 = rhs(P.x, P.V);
	k2 = rhs(P.x + h / 4, P.V + h / 4 * k1);
	k3 = rhs(P.x + h / 2, P.V + h / 2 * k2);
	k4 = rhs(P.x + h, P.V + h * (k1 + (-2) * k2 + 2 * k3));

	Res.V = P.V + h / 6 * (1 * k1 + 0 * k2 + 4 * k3 + 1 * k4);

	return Res;
}

std::vector < std::vector< point > > solve_ivp(const point& S, double h, double tol, const point& minP, const point& maxP) {
	double eps = 1e-4;
	bool withOLP = (tol != 0.0);
	tol = std::abs(tol);

	point maxPeps, minPeps;
	std::vector<double> One(S.V.size());
	for (size_t i = 0; i < One.size(); ++i) One[i] = 1.0;
	maxPeps.x = maxP.x - eps;
	maxPeps.V = maxP.V + (-eps) * One;
	minPeps.x = minP.x + eps;
	minPeps.V = minP.V + eps * One;
	std::vector<point> V;
	std::vector<point> E;
	std::vector<point> e_appr;
	std::vector<point> U;
	std::vector<point> V_half;
	int p = 4; // order

	point curP, tmp1P, tmp2P;
	curP = S;
	point curE;
	point cure;
	point curU;
	bool abort = false;
	bool next;
	bool pomitsya = withOLP;
	double Spar;
	curE.x = S.x;
	cure.x = S.x;
	curE.V = One - One;
	cure.V = One - One;
	E.push_back(curE);
	e_appr.push_back(cure);
	V.push_back(S);
	U.push_back(S);
	V_half.push_back(S);

//#pragma omp parallel
	while (!abort) {
		abort = false;
		next = false;
		while (!next) {
			next = true;
			
			// next point calculation
			curP = RK4(curP, h);
			if (isless_all(curP, maxPeps) && ismore_all(curP, minPeps)) {
			}
			else if (ismore_any(curP, maxP) || isless_any(curP, minP)) {
				curP = V[V.size() - 1];
				if (pomitsya) {
					h /= 2;
					continue;
				}
				else {
					abort = true;
					break;
				}
			}
			else {
				abort = true;
			}

			// OLP calculation
			tmp1P = RK4(V[V.size() - 1], h / 2);
			tmp2P = RK4(tmp1P, h / 2);
			cure.V = tmp2P.V - curP.V;
			if (withOLP) {
				Spar = norm2(cure.V) / ((1ull << p) - 1);
				if (Spar > tol && std::abs(h) > 1e-16) {
					h /= 2;
					next = false; // recalculate point
					curP = V[V.size() - 1];
				}
				else if (Spar < tol / (1ull << (p + 1))) {
					h *= 2;
				}
			}
		}

		cure.V = ((1ull << p) / ((1ull << p) - 1)) * (tmp2P.V - curP.V);
		curE.x = curP.x;
		cure.x = curP.x;
		curU.x = curP.x;
		curU.V = answer(curP.x, S.V);
		curE.V = curU.V - curP.V;
		if (withOLP) curP.V = curP.V + cure.V; //

		U.push_back(curU);
		E.push_back(curE);
		e_appr.push_back(cure);
		V.push_back(curP);
		V_half.push_back(tmp2P);
	}

	std::vector <std::vector <point> > ret;
	ret.push_back(V);
	ret.push_back(E);
	ret.push_back(e_appr);
	ret.push_back(U);
	ret.push_back(V_half);
	return ret;
}

int main(int argc, char **argv) {
	try {
		if (argc != 12) {
			std::cerr << "Invalid count of agruments, expetced: 12" << std::endl;
			return 0;
		}
		std::vector< std::string > input(argc);
		for (int i = 0; i < argc; ++i) input[i] = argv[i];

		point S;
		double h;
		double tol;
		point minP, maxP;

		S.x = std::stod(input[1]);
		S.V[0] = std::stod(input[2]);
		S.V[1] = std::stod(input[3]);
		h = std::stod(input[4]);
		tol = std::stod(input[5]);
		minP.x = std::stod(input[6]);
		minP.V[0] = std::stod(input[7]);
		minP.V[1] = std::stod(input[8]);
		maxP.x = std::stod(input[9]);
		maxP.V[0] = std::stod(input[10]);
		maxP.V[1] = std::stod(input[11]);

		//S.x = 0.0;
		//S.V[0] = 7.0;
		//S.V[1] = 13.0;
		//h = 0.001;
		//tol = 0.0;
		//minP.x = -0.1;
		//minP.V[0] = -1000.0;
		//minP.V[1] = -1000.0;
		//maxP.x = 0.1;
		//maxP.V[0] = 1000.0;
		//maxP.V[1] = 1000.0;

		std::vector <std::vector <point> > points = solve_ivp(S, h, tol, minP, maxP);
		if (points.size() == 0 || points[0].size() == 0) return 0;

		std::ofstream pFile("C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/pFile.txt");
		std::ofstream gTolFile("C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/gTolFile.txt");
		std::ofstream rTolFile("C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/rTolFile.txt");
		std::ofstream ansFile("C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/ansFile.txt");
		std::ofstream spravka("C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/spravka.txt");
		std::ofstream svodnaya_tablica("C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/svodnaya_tablica.txt");

		pFile << std::scientific << std::setprecision(8);
		gTolFile << std::scientific << std::setprecision(8);
		rTolFile << std::scientific << std::setprecision(8);
		ansFile << std::scientific << std::setprecision(8);
		spravka << std::scientific << std::setprecision(16);
		svodnaya_tablica << std::scientific << std::setprecision(16);
		double maxGlobTol = 0.0;
		double maxGlobTolX = 0.0;
		double howCloseToZero = std::numeric_limits<double>::infinity();

		svodnaya_tablica << "n" << '\t' << std::setw(22) << "Xn" << '\t' << std::setw(22) << "Vn1" << '\t' << std::setw(22) << "Vn2" << '\t' << std::setw(22) << "V^n1" << '\t' << std::setw(22) << "V^n2" << '\t' << std::setw(22) << "En (norm2)" << '\t' << std::setw(22) << "S* (norm2)" << '\t' << std::setw(22) << "Hn-1" << std::endl;
		size_t n = 0;
		svodnaya_tablica << n << '\t' << points[0][0].x << '\t' << points[0][0].V[0] << '\t' << points[0][0].V[1] << '\t' << "-" << '\t' << std::setw(18) << "-" << '\t' << std::setw(18) << 0.0 << '\t' << 0.0 << '\t' << "-" << std::endl;
		for (n = 1; n < points[0].size(); ++n) {
			svodnaya_tablica << n << '\t' << points[0][n].x << '\t' << points[0][n].V[0] << '\t' << points[0][n].V[1] << '\t' << points[4][n].V[0] << '\t' << points[4][n].V[1] << '\t' << norm2(points[1][n].V) << '\t' << norm2(points[2][n].V) << '\t' << points[0][n].x - points[0][n - 1].x << std::endl;
		}

		for (size_t i = 0; i < points[0].size(); ++i) {
			pFile << points[0][i].x << "\t" << points[0][i].V[0] << "\t" << points[0][i].V[1] << std::endl;
			gTolFile << points[1][i].x << "\t" << points[1][i].V[0] << "\t" << points[1][i].V[1] << std::endl;
			rTolFile << points[2][i].x << "\t" << points[2][i].V[0] << "\t" << points[2][i].V[1] << std::endl;
			ansFile << points[3][i].x << "\t" << points[3][i].V[0] << "\t" << points[3][i].V[1] << std::endl;
			if (norm2(points[1][i].V) > maxGlobTol) {
				maxGlobTol = norm2(points[1][i].V);
				maxGlobTolX = points[1][i].x;
			}
			howCloseToZero = norm2(points[0][i].V);
		}
		spravka << "Count of steps: " << points[0].size() - 1 << '\n' << "Global tolerance: point and max value: " << maxGlobTolX << '\t' << maxGlobTol << '\n' << "How close to zero: " << howCloseToZero << std::endl;

		pFile.close();
		gTolFile.close();
		rTolFile.close();
		ansFile.close();
		spravka.close();
		svodnaya_tablica.close();
	}
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
	}

	return 0;
}

