// mesh_partition.cpp : 再構築版 + VALUE属性のXMF出力 + globalCellIDs出力 + /COMM 構造追加

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <string>
#include <fstream>
#include <filesystem>
#include <metis.h>

#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

using namespace std;
using namespace HighFive;

using geom_int = int;
using geom_float = float;

//-------------------------------------------------------------
// ユーティリティ関数
//-------------------------------------------------------------
set<geom_int> extract_nodes_from_conne(const vector<vector<geom_int>>& conne_cells, const vector<geom_int>& cell_ids);
unordered_map<geom_int, geom_int> build_local_node_map(const set<geom_int>& global_nodes);
void write_xmf(const string& xmf_name, int nParts,
               const vector<geom_int>& nNodes_per_rank,
               const vector<geom_int>& nCells_per_rank,
               const vector<geom_int>& conne_sizes,
               const vector<string>& value_keys);

//-------------------------------------------------------------
// main 関数
//-------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: ./mesh_partition input.h5 nParts" << endl;
        return 1;
    }

    string input_file = argv[1];
    int nParts = stoi(argv[2]);
    cout << "[INFO] Loading mesh from " << input_file << endl;

    File file(input_file, File::ReadOnly);
    Group mesh_group = file.getGroup("/MESH");

    geom_int nNodes = mesh_group.getAttribute("nNodes").read<geom_int>();
    geom_int nCells = mesh_group.getAttribute("nCells").read<geom_int>();

    vector<geom_float> coords;
    file.getDataSet("/MESH/COORD").read(coords);

    vector<geom_int> conne_raw;
    file.getDataSet("/MESH/CONNE").read(conne_raw);

    vector<vector<geom_int>> conne_cells(nCells);
    size_t pos = 0;
    for (geom_int i = 0; i < nCells; ++i) {
        geom_int type = conne_raw[pos++];
        int npts = (type == 9 ? 8 : type == 8 ? 6 : type == 7 ? 5 : 4);
        vector<geom_int> line;
        line.push_back(type);
        for (int j = 0; j < npts; ++j)
            line.push_back(conne_raw[pos++]);
        conne_cells[i] = line;
    }

    vector<geom_float> volume;
    file.getDataSet("/CELLS/volume").read(volume);
    vector<geom_float> centroids;
    file.getDataSet("/CELLS/centCoords").read(centroids);

    // `/CELLS/STRUCT` を直接読み込む
    vector<geom_int> cell_struct;
    file.getDataSet("/CELLS/STRUCT").read(cell_struct);

    // セルごとの情報を解析
    pos = 0;
    vector<vector<geom_int>> parsed_struct(nCells);
    for (geom_int i = 0; i < nCells; ++i) {
        vector<geom_int> cell_info;

        // ノード数とノードID
        geom_int nNode = cell_struct[pos++];
        cell_info.push_back(nNode);
        for (geom_int j = 0; j < nNode; ++j) {
            cell_info.push_back(cell_struct[pos++]);
        }

        // 面数と面ID
        geom_int nPlane = cell_struct[pos++];
        cell_info.push_back(nPlane);
        for (geom_int j = 0; j < nPlane; ++j) {
            cell_info.push_back(cell_struct[pos++]);
        }

        // 面符号
        geom_int nPlaneSign = cell_struct[pos++];
        cell_info.push_back(nPlaneSign);
        for (geom_int j = 0; j < nPlaneSign; ++j) {
            cell_info.push_back(cell_struct[pos++]);
        }

        // セルタイプ
        geom_int cellType = cell_struct[pos++];
        cell_info.push_back(cellType);

        parsed_struct[i] = cell_info;
    }

    map<string, vector<geom_float>> value_fields;
    vector<string> value_keys;
    if (file.exist("/VALUE")) {
        Group valgrp = file.getGroup("/VALUE");
        for (const auto& name : valgrp.listObjectNames()) {
            vector<geom_float> field;
            valgrp.getDataSet(name).read(field);
            value_fields[name] = field;
            value_keys.push_back(name);
        }
    }

    map<geom_int, vector<geom_int>> node_to_cells;
    for (geom_int ic = 0; ic < nCells; ++ic) {
        for (size_t j = 1; j < conne_cells[ic].size(); ++j)
            node_to_cells[conne_cells[ic][j]].push_back(ic);
    }

    vector<set<geom_int>> cell_adj(nCells);
    for (const auto& [nid, connected] : node_to_cells) {
        for (size_t i = 0; i < connected.size(); ++i) {
            for (size_t j = i + 1; j < connected.size(); ++j) {
                cell_adj[connected[i]].insert(connected[j]);
                cell_adj[connected[j]].insert(connected[i]);
            }
        }
    }

    vector<idx_t> xadj(nCells + 1);
    vector<idx_t> adjncy;
    for (geom_int i = 0; i < nCells; ++i) {
        xadj[i+1] = xadj[i] + cell_adj[i].size();
        for (auto nei : cell_adj[i])
            adjncy.push_back(nei);
    }

    // *** METIS での分割 ***
    vector<idx_t> part(nCells);
    idx_t objval;
    idx_t ncon = 1;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_UFACTOR] = 30;

    cout << "[INFO] Partitioning with METIS..." << endl;
    METIS_PartGraphKway(&nCells, &ncon, xadj.data(), adjncy.data(), NULL, NULL, NULL,
                        &nParts, NULL, NULL, options, &objval, part.data());
    cout << "[INFO] METIS partitioning done." << endl;

    vector<vector<geom_int>> part_cells(nParts);
    for (geom_int i = 0; i < nCells; ++i)
        part_cells[part[i]].push_back(i);


    vector<geom_int> nNodes_per_rank(nParts);
    vector<geom_int> nCells_per_rank(nParts);
    vector<geom_int> conne_sizes(nParts);

    filesystem::create_directory("partitioned");


    // 各プロセスのデータを保存
    for (int rank = 0; rank < nParts; ++rank) {
        const auto& owned = part_cells[rank];
        set<geom_int> owned_set(owned.begin(), owned.end());
        set<geom_int> halo_set;

        for (auto cid : owned) {
            for (size_t j = 1; j < conne_cells[cid].size(); ++j) {
                geom_int nid = conne_cells[cid][j];
                for (auto adj_cid : node_to_cells[nid]) {
                    if (!owned_set.count(adj_cid))
                        halo_set.insert(adj_cid);
                }
            }
        }

        vector<geom_int> all_cells(owned.begin(), owned.end());
        all_cells.insert(all_cells.end(), halo_set.begin(), halo_set.end());

        set<geom_int> used_nodes = extract_nodes_from_conne(conne_cells, all_cells);
        auto global_to_local = build_local_node_map(used_nodes);

        vector<geom_float> local_coords(global_to_local.size() * 3);
        for (const auto& [gid, lid] : global_to_local)
            for (int d = 0; d < 3; ++d)
                local_coords[3 * lid + d] = coords[3 * gid + d];

        vector<geom_int> conne_out;
        for (auto cid : all_cells) {
            conne_out.push_back(conne_cells[cid][0]);
            for (size_t j = 1; j < conne_cells[cid].size(); ++j)
                conne_out.push_back(global_to_local[conne_cells[cid][j]]);
        }

        vector<geom_float> local_volume;
        vector<geom_float> local_centroids;
        vector<geom_int> globalCellIDs;
        for (auto cid : all_cells) {
            local_volume.push_back(volume[cid]);
            for (int d = 0; d < 3; ++d)
                local_centroids.push_back(centroids[3 * cid + d]);
            globalCellIDs.push_back(cid);
        }

        // グローバルノードIDを保存
        vector<geom_int> globalNodeIDs(global_to_local.size());
        for (const auto& [gid, lid] : global_to_local) {
            globalNodeIDs[lid] = gid; // ローカルIDに対応するグローバルIDを保存
        }
        
        map<string, vector<geom_float>> local_values;
        for (const auto& [name, data] : value_fields) {
            vector<geom_float> part_data;
            for (auto cid : all_cells)
                part_data.push_back(data[cid]);
            local_values[name] = part_data;
        }

        // `/CELLS/STRUCT` の分割と保存
        vector<geom_int> local_struct;
        for (auto cid : all_cells) {
            size_t pos = 0;
        
            // ノード数とノードIDを追加
            geom_int nNode = parsed_struct[cid][pos++];
            local_struct.push_back(nNode);
            for (geom_int j = 0; j < nNode; ++j) {
                geom_int global_node_id = parsed_struct[cid][pos++];
                local_struct.push_back(global_to_local[global_node_id]); // ローカルノードIDに変換
            }
        
            // 面数と面IDを追加
            geom_int nPlane = parsed_struct[cid][pos++];
            local_struct.push_back(nPlane);
            for (geom_int j = 0; j < nPlane; ++j) {
                geom_int plane_id = parsed_struct[cid][pos++];
                local_struct.push_back(plane_id); // 面ID（グローバルIDのまま）
            }
        
            // 面符号を追加
            geom_int nPlaneSign = parsed_struct[cid][pos++];
            local_struct.push_back(nPlaneSign);
            for (geom_int j = 0; j < nPlaneSign; ++j) {
                geom_int plane_sign = parsed_struct[cid][pos++];
                local_struct.push_back(plane_sign); // 面符号（±1）
            }
        
            // セルタイプを追加
            geom_int cellType = parsed_struct[cid][pos++];
            local_struct.push_back(cellType);
        }
        
        string fname = "partitioned/rank_" + to_string(rank) + ".h5";
        File outf(fname, File::Overwrite);

        Group mgrp = outf.createGroup("/MESH");
        mgrp.createAttribute("nNodes", static_cast<geom_int>(used_nodes.size()));
        mgrp.createAttribute("nCells", static_cast<geom_int>(all_cells.size()));
        mgrp.createDataSet("globalNodeIDs", globalNodeIDs);

        outf.createDataSet("/MESH/COORD", local_coords);
        outf.createDataSet("/MESH/CONNE", conne_out);


        Group cgrp = outf.createGroup("/CELLS");
        cgrp.createDataSet("volume", local_volume);
        cgrp.createDataSet("centCoords", local_centroids);
        cgrp.createDataSet("globalCellIDs", globalCellIDs);
        cgrp.createDataSet("STRUCT", local_struct);                                                                                                                                        

        Group vgrp = outf.createGroup("/VALUE");
        for (const auto& [name, vals] : local_values) {
            vgrp.createDataSet(name, vals);
        }

        nNodes_per_rank[rank] = used_nodes.size();
        nCells_per_rank[rank] = all_cells.size();
        conne_sizes[rank] = conne_out.size();

        cout << "[rank " << rank << "] written to " << fname << endl;
    }

    write_xmf("partition.xmf", nParts, nNodes_per_rank, nCells_per_rank, conne_sizes, value_keys);
    return 0;
}

set<geom_int> extract_nodes_from_conne(const vector<vector<geom_int>>& conne_cells, const vector<geom_int>& cell_ids) {
    set<geom_int> used_nodes;
    for (auto cid : cell_ids)
        for (size_t i = 1; i < conne_cells[cid].size(); ++i)
            used_nodes.insert(conne_cells[cid][i]);
    return used_nodes;
}

unordered_map<geom_int, geom_int> build_local_node_map(const set<geom_int>& global_nodes) {
    unordered_map<geom_int, geom_int> global_to_local;
    geom_int lid = 0;
    for (auto gid : global_nodes)
        global_to_local[gid] = lid++;
    return global_to_local;
}

void write_xmf(const string& xmf_name, int nParts,
               const vector<geom_int>& nNodes_per_rank,
               const vector<geom_int>& nCells_per_rank,
               const vector<geom_int>& conne_sizes,
               const vector<string>& value_keys) {
    ofstream xmf(xmf_name);
    xmf << "<?xml version=\"1.0\" ?>\n";
    xmf << "<Xdmf Version=\"3.0\">\n";
    xmf << "  <Domain>\n";
    xmf << "    <Grid Name=\"PartitionedMesh\" GridType=\"Collection\" CollectionType=\"Spatial\">\n";

    for (int rank = 0; rank < nParts; ++rank) {
        string fname = "partitioned/rank_" + to_string(rank) + ".h5";
        xmf << "      <Grid Name=\"rank_" << rank << "\" GridType=\"Uniform\">\n";
        xmf << "        <Topology TopologyType=\"Mixed\" NumberOfElements=\""
            << nCells_per_rank[rank] << "\">\n";
        xmf << "          <DataItem Dimensions=\"" << conne_sizes[rank] << "\" Format=\"HDF\">\n";
        xmf << "            " << fname << ":/MESH/CONNE\n";
        xmf << "          </DataItem>\n";
        xmf << "        </Topology>\n";
        xmf << "        <Geometry GeometryType=\"XYZ\">\n";
        xmf << "          <DataItem Dimensions=\"" << nNodes_per_rank[rank] << " 3\" Format=\"HDF\">\n";
        xmf << "            " << fname << ":/MESH/COORD\n";
        xmf << "          </DataItem>\n";
        xmf << "        </Geometry>\n";

        for (const auto& name : value_keys) {
            xmf << "        <Attribute Name=\"" << name << "\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
            xmf << "          <DataItem Dimensions=\"" << nCells_per_rank[rank] << "\" Format=\"HDF\">\n";
            xmf << "            " << fname << ":/VALUE/" << name << "\n";
            xmf << "          </DataItem>\n";
            xmf << "        </Attribute>\n";
        }

        xmf << "      </Grid>\n";
    }

    xmf << "    </Grid>\n";
    xmf << "  </Domain>\n";
    xmf << "</Xdmf>\n";
    xmf.close();
    cout << "[DONE] XMF file written to " << xmf_name << endl;
}
