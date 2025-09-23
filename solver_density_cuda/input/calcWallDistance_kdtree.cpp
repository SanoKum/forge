#include "input/calcWallDistance_kdtree.hpp"

namespace {
inline geom_float dist(const Point &a, const Point &b) {
    const geom_float dx = a.x - b.x;
    const geom_float dy = a.y - b.y;
    const geom_float dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

#ifdef HAVE_KDTREE
void compute_with_kdtree(const std::vector<Point> &walls,
                                                 const std::vector<Point> &cells,
                                                 std::vector<geom_float> &distance) {
    kdtree *tree = kd_create(3);
    boost::scoped_array<int> indices(new int[walls.size()]);
    for (int i = 0; i < static_cast<int>(walls.size()); ++i) {
        indices[i] = i;
        kd_insert3(tree, walls[i].x, walls[i].y, walls[i].z, &indices[i]);
    }
    distance.reserve(cells.size());
    for (int j = 0; j < static_cast<int>(cells.size()); ++j) {
        kdres *set = kd_nearest3(tree, cells[j].x, cells[j].y, cells[j].z);
        int i = *(int *)kd_res_item_data(set);
        kd_res_free(set);
        distance.push_back(dist(walls[i], cells[j]));
    }
    kd_free(tree);
}
#endif

void compute_bruteforce(const std::vector<Point> &walls,
                                                const std::vector<Point> &cells,
                                                std::vector<geom_float> &distance) {
    distance.reserve(cells.size());
    for (const auto &c : cells) {
        geom_float dmin = std::numeric_limits<geom_float>::max();
        for (const auto &w : walls) {
            dmin = std::min(dmin, dist(w, c));
        }
        distance.push_back(dmin == std::numeric_limits<geom_float>::max() ? 0.0 : dmin);
    }
}
} // namespace

void calcWallDistance_kdtree(solverConfig &cfg, mesh &msh, variables &var) {
    (void)cfg; // currently unused

    std::vector<Point> wall_points;
    std::vector<Point> cell_points;
    std::vector<geom_float> distance;

    int wall_count = 0;
    for (auto &bc : msh.bconds) {
        if (bc.bcondKind == "wall_isothermal" || bc.bcondKind == "wall") {
            for (auto &ip : bc.iPlanes) {
                geom_float x = msh.planes[ip].centCoords[0];
                geom_float y = msh.planes[ip].centCoords[1];
                geom_float z = msh.planes[ip].centCoords[2];
                wall_points.emplace_back(x, y, z);
            }
            ++wall_count;
        }
    }

    if (wall_count > 0) {
        cell_points.reserve(msh.cells.size());
        for (auto &icell : msh.cells) {
            geom_float x = icell.centCoords[0];
            geom_float y = icell.centCoords[1];
            geom_float z = icell.centCoords[2];
            cell_points.emplace_back(x, y, z);
        }

#ifdef HAVE_KDTREE
        compute_with_kdtree(wall_points, cell_points, distance);
#else
        compute_bruteforce(wall_points, cell_points, distance);
#endif

        // Copy into variable array (assumes var.c["wall_dist"] sized nCells_all)
        if (var.c.count("wall_dist") && var.c["wall_dist"].size() >= distance.size()) {
            std::copy(distance.begin(), distance.end(), var.c["wall_dist"].begin());
        }
    } else {
        distance.assign(msh.nCells_all, geom_float(0));
        if (var.c.count("wall_dist") && var.c["wall_dist"].size() >= distance.size()) {
            std::copy(distance.begin(), distance.end(), var.c["wall_dist"].begin());
        }
    }

    std::list<std::string> names = {"wall_dist"};
    var.copyVariables_cell_H2D(names);
}