#include <cmath>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <limits>
#include <boost/scoped_array.hpp>
#include </home/kumpei/app/kdtree/install/include/kdtree.h>

#include "flowFormat.hpp"
#include "input/solverConfig.hpp"
#include "mesh/mesh.hpp"
#include "variables.hpp"

struct Point {
  double x, y, z;
  Point(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

double norm2(const Point &l, const Point &r) ;

static double useKDTree(const std::vector<Point> &left, const std::vector<Point> &right);

void calcWallDistance_kdtree(solverConfig& cfg , mesh& msh , variables& var );