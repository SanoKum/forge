#include <iostream>
#include <vector>
#include <set>
#include <sstream>
#include <stdexcept>
#include "H5Cpp.h"
#include "metis.h"

using namespace H5;
using namespace std;

int main(int argc, char* argv[]) {
    // --- 引数による設定 ---
    string filename;
    idx_t nparts;

    // ファイル名の取得
    if (argc > 1) {
        filename = argv[1];
    } else {
        cout << "メッシュファイル名 (HDF5形式) を入力してください: ";
        cin >> filename;
        if (filename.empty()) {
            cerr << "ファイル名が入力されていません。" << endl;
            return -1;
        }
    }

    // パーティション数の取得
    if (argc > 2) {
        stringstream ss(argv[2]);
        ss >> nparts;
        if (!ss || nparts < 1) {
            cerr << "無効なパーティション数が指定されました: " << argv[2] << endl;
            return -1;
        }
    } else {
        cout << "パーティション数を入力してください: ";
        cin >> nparts;
        if (!cin || nparts < 1) {
            cerr << "無効なパーティション数が入力されました" << endl;
            return -1;
        }
    }

    try {
        // --- Step 1: HDF5 からメッシュデータの読み込み ---
        // 指定されたファイルを読み込み（データセット "elements" と "nodes" を仮定）
        H5File file(filename, H5F_ACC_RDONLY);

        // セル（要素）の接続情報を読み込み
        DataSet elementDataset = file.openDataSet("elements");
        DataSpace elementSpace = elementDataset.getSpace();
        hsize_t elementDims[2];
        elementSpace.getSimpleExtentDims(elementDims, NULL);
        int numElements = static_cast<int>(elementDims[0]);
        int nodesPerElement = static_cast<int>(elementDims[1]);
        vector<int> elements(numElements * nodesPerElement);
        elementDataset.read(elements.data(), PredType::NATIVE_INT);

        // ノード座標の読み込み（"nodes" データセット、各ノードは (x,y,z) と仮定）
        DataSet nodeDataset = file.openDataSet("nodes");
        DataSpace nodeSpace = nodeDataset.getSpace();
        hsize_t nodeDims[2];
        nodeSpace.getSimpleExtentDims(nodeDims, NULL);
        int numNodes = static_cast<int>(nodeDims[0]);
        vector<double> nodes(numNodes * 3);
        nodeDataset.read(nodes.data(), PredType::NATIVE_DOUBLE);

        cout << "読み込み完了: セル数 = " << numElements 
             << ", ノード数 = " << numNodes 
             << ", 各セルのノード数 = " << nodesPerElement << endl;

        // --- Step 2: メッシュからグラフ表現の作成 ---
        // 各ノードが属するセルをマッピング
        vector<vector<int>> nodeToCell(numNodes);
        for (int elem = 0; elem < numElements; elem++) {
            for (int j = 0; j < nodesPerElement; j++) {
                int node = elements[elem * nodesPerElement + j];
                if (node < 0 || node >= numNodes) {
                    throw runtime_error("ノード番号が範囲外です。");
                }
                nodeToCell[node].push_back(elem);
            }
        }
        // 各セルについて、共通ノードを持つセルを隣接セルと定義
        vector<set<int>> cellAdj(numElements);
        for (int node = 0; node < numNodes; node++) {
            for (size_t i = 0; i < nodeToCell[node].size(); i++) {
                for (size_t j = 0; j < nodeToCell[node].size(); j++) {
                    if (i == j) continue;
                    int cell_i = nodeToCell[node][i];
                    int cell_j = nodeToCell[node][j];
                    cellAdj[cell_i].insert(cell_j);
                }
            }
        }
        cout << "グラフ表現の作成完了" << endl;

        // --- Step 3: METIS 用の入力データ作成と領域分割 ---
        // METIS 用の配列 xadj と adjncy を作成
        vector<idx_t> xadj(numElements + 1, 0);
        vector<idx_t> adjncy;
        for (int i = 0; i < numElements; i++) {
            xadj[i + 1] = xadj[i] + cellAdj[i].size();
            for (int neighbor : cellAdj[i]) {
                adjncy.push_back(neighbor);
            }
        }
        
        // 領域分割結果を格納する配列
        vector<idx_t> part(numElements, 0);
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        idx_t ncon = 1;  // 制約数（通常1）
        idx_t objval;    // カットエッジ重み

        // METIS_PartGraphKway により領域分割を実行
        int ret = METIS_PartGraphKway(&numElements,
                                      &ncon,
                                      xadj.data(),
                                      adjncy.data(),
                                      NULL,   // 頂点重み（NULLなら均等）
                                      NULL,   // 頂点サイズ
                                      NULL,   // 辺重み
                                      &nparts,
                                      NULL,   // 目標パーティションサイズ（NULLなら均等）
                                      NULL,   // 許容不均衡率
                                      options,
                                      &objval,
                                      part.data());
        if(ret != METIS_OK) {
            cerr << "METIS による領域分割に失敗しました！" << endl;
            return -1;
        }
        cout << "METIS による領域分割完了 (カットエッジ重み: " << objval << ")" << endl;

        // 分割結果をもとに、各パーティションのゴーストセル（境界セル）の抽出
        vector<set<int>> ghostCells(nparts);
        for (int elem = 0; elem < numElements; elem++) {
            for (int neighbor : cellAdj[elem]) {
                if (part[elem] != part[neighbor]) {
                    ghostCells[part[elem]].insert(neighbor);
                }
            }
        }

        // 結果の表示
        for (int elem = 0; elem < numElements; elem++) {
            cout << "セル " << elem << " -> パーティション " << part[elem] << endl;
        }
        for (idx_t p = 0; p < nparts; p++) {
            cout << "パーティション " << p << " のゴーストセル: ";
            for (int ghost : ghostCells[p]) {
                cout << ghost << " ";
            }
            cout << endl;
        }
    }
    catch (FileIException &error) {
        error.printErrorStack();
        return -1;
    }
    catch (DataSetIException &error) {
        error.printErrorStack();
        return -1;
    }
    catch (DataSpaceIException &error) {
        error.printErrorStack();
        return -1;
    }
    catch (exception &ex) {
        cerr << "エラー: " << ex.what() << endl;
        return -1;
    }
    
    return 0;
}
