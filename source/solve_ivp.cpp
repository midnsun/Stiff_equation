#include <vector>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>

std::vector<double> operator+ (std::vector<double> l, std::vector<double> r) {
	std::vector<double> res(l.size());
	for (size_t i = 0; i < res.size(); ++i) res[i] = l[i] + r[i];
	return res;
}

std::vector<double> operator* (double a, std::vector<double> r) {
	std::vector<double> res(r.size());
	for (size_t i = 0; i < res.size(); ++i) res[i] = a * r[i];
	return res;
}

struct point {
	double t;
	std::vector<double> X;
	point() : t(0.0) {
		X.resize(2);
	}
	point(const point& p) : t(p.t), X(p.X) {}
	point(point&& p) : t(p.t), X(std::move(p.X)) {}
	point& operator= (const point& p) {
		t = p.t;
		X = p.X;
		return *this;
	}
	point& operator= (point&& p) {
		t = p.t;
		X = p.X;
		return *this;
	}
};

std::vector<double> rhs(double t, const std::vector<double>& X) {
	std::vector<double> Y(X.size());
	std::vector< std::vector<double> > M(X.size(), std::vector<double>(X.size()));

	M[0][0] = -500.005;		M[0][1] = 499.995;
	M[1][0] = 499.995;		M[1][1] = -500.005;

	for (size_t i = 0; i < M.size(); ++i) {
		for (size_t j = 0; j < M[i].size(); ++j) {
			Y[i] += M[i][j] * X[j];
		}
	}

	return Y;
}

std::vector<double> rhs(const point& P) {
	return rhs(P.t, P.X);
}

std::vector<double> answer(double t, const std::vector<double>& S) {
	std::vector<double> Y(S.size());
	double C2 = (S[0] + S[1]) / 2, C1 = (S[1] - S[0]) / 2;
	std::vector<double> v1(2), v2(2);

	v1[0] = -1.0;	v2[0] = 1.0;
	v1[1] = 1.0;	v2[1] = 1.0;

	Y = C1 * exp(-1000 * t) * v1 + C2 * exp(-0.01 * t) * v2;

	return Y;
}

bool isless(const point& x, const point& l) {
	if (x.t > l.t) return false;
	for (size_t i = 0; i < x.X.size(); ++i) if (x.X[i] > l.X[i]) return false;

	return true;
}

bool ismore(const point& x, const point& r) {
	if (x.t < r.t) return false;
	for (size_t i = 0; i < x.X.size(); ++i) if (x.X[i] < r.X[i]) return false;

	return true;
}


point RK4(const point& P, double h) {
	point Res;
	Res.t = P.t + h;
	size_t n = Res.X.size();
	
	std::vector<double> k1(n);
	std::vector<double> k2(n);
	std::vector<double> k3(n);
	std::vector<double> k4(n);

	k1 = rhs(P.t, P.X);
	k2 = rhs(P.t + h / 4, P.X + h / 4 * k1);
	k3 = rhs(P.t + h / 2, P.X + h / 2 * k2);
	k4 = rhs(P.t + h, P.X + h * (k1 + (-2) * k2 + 2 * k3));

	Res.X = P.X + h / 6 * (1 * k1 + 0 * k2 + 4 * k3 + 1 * k4);

	return Res;
}

double norm2(const std::vector<double>& v) {
	double res = 0.0;
	for (size_t i = 0; i < v.size(); ++i) {
		res += v[i] * v[i];
	}
	return std::sqrt(res);
}

std::vector < std::vector< point > > solve_ivp(const point& S, double h, double tol, const point& minP, const point& maxP) {
	double eps = 1e-4;
	bool withOLP = (tol != 0.0);
	tol = std::abs(tol);
	if (tol < 1e-12 && withOLP) throw std::logic_error("Too small tolerance");

	point maxPeps, minPeps;
	std::vector<double> One(S.X.size());
	for (size_t i = 0; i < One.size(); ++i) One[i] = 1.0;
	maxPeps.t = maxP.t - eps;
	maxPeps.X = maxP.X + (-eps) * One;
	minPeps.t = minP.t + eps;
	minPeps.X = minP.X + eps * One;
	std::vector<point> res;
	std::vector<point> globalTol;
	std::vector<point> relatieTol;
	std::vector<point> ans;
	int p = 4; // order

	point curP, tmp1P, tmp2P;
	curP = S;
	point curGTol;
	point curRTol;
	point curAns;
	bool abort = false;
	bool next;
	double OLP;
	curGTol.t = S.t;
	curRTol.t = S.t;
	curGTol.X = One + (-1) * One;
	curRTol.X = One + (-1) * One;
	globalTol.push_back(curGTol);
	relatieTol.push_back(curRTol);
	res.push_back(S);
	ans.push_back(S);

	while (!abort) {
		abort = false;
		next = false;
		while (!next) {
			next = true;
			
			// next point calculation
			curP = RK4(curP, h);
			if (isless(curP, maxPeps) && ismore(curP, minPeps)) {
			}
			else if (ismore(curP, maxP) || isless(curP, minP)) {
				abort = true;
				break;
			}
			else {
				abort = true;
			}

			// OLP calculation
			tmp1P = RK4(res[res.size() - 1], h / 2);
			tmp2P = RK4(tmp1P, h / 2);
			curRTol.X = tmp2P.X + (-1) * curP.X;
			if (withOLP) {
				OLP = norm2(curRTol.X) / ((1ull << p) - 1);
				if (OLP > tol) {
					h /= 2;
					next = false; // recalculate point
					curP = res[res.size() - 1];
				}
				else if (OLP < tol / (1ull << (p + 1))) {
					h *= 2;
				}
			}
		}

		curGTol.t = curP.t;
		curRTol.t = curP.t;
		curAns.t = curP.t;
		curAns.X = answer(curP.t, S.X);
		curGTol.X = curAns.X + (-1) * curP.X;
		ans.push_back(curAns);
		globalTol.push_back(curGTol);
		relatieTol.push_back(curRTol);
		res.push_back(curP);
	}

	std::vector <std::vector <point> > ret;
	ret.push_back(res);
	ret.push_back(globalTol);
	ret.push_back(relatieTol);
	ret.push_back(ans);
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

		S.t = std::stod(input[1]);
		S.X[0] = std::stod(input[2]);
		S.X[1] = std::stod(input[3]);
		h = std::stod(input[4]);
		tol = std::stod(input[5]);
		minP.t = std::stod(input[6]);
		minP.X[0] = std::stod(input[7]);
		minP.X[1] = std::stod(input[8]);
		maxP.t = std::stod(input[9]);
		maxP.X[0] = std::stod(input[10]);
		maxP.X[1] = std::stod(input[11]);

		//S.t = 0.0;
		//S.X[0] = 7.0;
		//S.X[1] = 13.0;
		//h = 0.001;
		//tol = 0.0;
		//minP.t = -0.1;
		//minP.X[0] = -1000.0;
		//minP.X[1] = -1000.0;
		//maxP.t = 0.1;
		//maxP.X[0] = 1000.0;
		//maxP.X[1] = 1000.0;

		std::vector <std::vector <point> > points = solve_ivp(S, h, tol, minP, maxP);

		std::ofstream pFile("C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/pFile.txt");
		std::ofstream gTolFile("C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/gTolFile.txt");
		std::ofstream rTolFile("C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/rTolFile.txt");
		std::ofstream ansFile("C:/Users/chehp/OneDrive/Desktop/all/SE/KSR1/data/ansFile.txt");

		pFile << std::scientific << std::setprecision(8);
		gTolFile << std::scientific << std::setprecision(8);
		rTolFile << std::scientific << std::setprecision(8);
		ansFile << std::scientific << std::setprecision(8);

		for (size_t i = 0; i < points[0].size(); ++i) {
			pFile << points[0][i].t << "\t" << points[0][i].X[0] << "\t" << points[0][i].X[1] << std::endl;
			gTolFile << points[1][i].t << "\t" << points[1][i].X[0] << "\t" << points[1][i].X[1] << std::endl;
			rTolFile << points[2][i].t << "\t" << points[2][i].X[0] << "\t" << points[2][i].X[1] << std::endl;
			ansFile << points[3][i].t << "\t" << points[3][i].X[0] << "\t" << points[3][i].X[1] << std::endl;
		}

		pFile.close();
		gTolFile.close();
		rTolFile.close();
		ansFile.close();
	}
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
	}

	return 0;
}

