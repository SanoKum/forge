#include "input/calcWallDistance_kdtree.hpp"

geom_float norm2_(const Point &l, const Point &r) {
  const geom_float dx = l.x - r.x;
  const geom_float dy = l.y - r.y;
  const geom_float dz = l.z - r.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static void useKDTree(const std::vector<Point> &left, const std::vector<Point> &right, std::vector<geom_float> &distance) {
    kdtree *tree = kd_create(3);

    boost::scoped_array<int> indexes(new int[left.size()]);
    for (int i = 0; i < (int)left.size(); ++i) {
        indexes[i] = i;
        kd_insert3(tree, left[i].x, left[i].y, left[i].z, &indexes[i]);
    }

    int l, r;
    geom_float minimum = std::numeric_limits<geom_float>::max();
    for (int j = 0; j < (int)right.size(); ++j) {
        kdres *set = kd_nearest3(tree, right[j].x, right[j].y, right[j].z);

        int i = *(int *)kd_res_item_data(set);
        kd_res_free(set);
        geom_float d = norm2_(left[i], right[j]);

        distance.push_back(d);
    }

    kd_free(tree);
    //return minimum;
}

void calcWallDistance_kdtree(solverConfig& cfg , mesh& msh , variables& var) {

    std::vector<Point> wall_points; 
    std::vector<Point> cell_points; 
    std::vector<geom_float> distance; 

    int wall_count = 0;
    for (auto& bc : msh.bconds) {
        if (bc.bcondKind == "wall_isothermal" or bc.bcondKind == "wall") {
            for (auto& ip : bc.iPlanes) {
                geom_float x = msh.planes[ip].centCoords[0];
                geom_float y = msh.planes[ip].centCoords[1];
                geom_float z = msh.planes[ip].centCoords[2];
                wall_points.emplace_back(x,y,z);
            }
            wall_count += 1;
        } 
    }

    if (wall_count > 0) {
        for (auto& icell : msh.cells) {
            geom_float x = icell.centCoords[0];
            geom_float y = icell.centCoords[1];
            geom_float z = icell.centCoords[2];
            cell_points.emplace_back(x,y,z);
        }

        useKDTree(wall_points, cell_points, distance);

        std::copy(distance.begin(), distance.end(), var.c["wall_dist"].begin());
    } else {
        distance.resize(msh.nCells_all);
        geom_float zero = 0.0;
        std::fill(distance.begin(), distance.end(), zero);
    }

    std::list<std::string> names = {"wall_dist"};
    var.copyVariables_cell_H2D(names);
}