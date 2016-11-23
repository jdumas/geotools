////////////////////////////////////////////////////////////////////////////////
#include "poisson_disk.hpp"
#include "vec.h"
#include "chrono.h"
// -----------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

typedef PoissonSampling<double, 3> Sampler3D;

void write_samples(const std::vector<Vec3d> &pts, std::ostream &out) {
	out << pts.size() << std::endl;
	for (Vec3d x : pts) {
		out << x[0] << " " << x[1] << " " << x[2] << std::endl;
	}
	out << std::endl;
}

int main(int argc, char** argv) {
	if (argc < 3) {
		std::cout << "Usage: " << argv[0] << " min_dist out_file" << std::endl;
		return 0;
	}

	Chrono timer("Timer");

	Vec3d unit_box = Vec::constant<double, 3>(1);
	double min_dist = std::stod(argv[1]);

	std::cout << "- Initialization..." << std::endl;
	Sampler3D sampler(min_dist, unit_box);

	std::cout << "- Sampling (Naive)..." << std::endl;
	std::vector<Vec3d> result_naive;
	timer.tic();
	sampler.naive(result_naive);
	timer.toc();

	std::cout << "- Sampling (Bridson)..." << std::endl;
	std::vector<Vec3d> result;
	timer.tic();
	sampler.box(30, result);
	timer.toc();

	std::cout << "- Saving result..." << std::endl;
	std::ofstream fout(argv[2]);
	write_samples(result, fout);

	return 0;
}
