Project: solver_density_cuda — Developer Guidelines

Scope
- This document captures project-specific build, configuration, testing, and development practices to speed up onboarding and avoid common pitfalls. It assumes an advanced C++/CUDA developer using CLion with the pre-configured CMake profile provided in this repository.

Environment and External Dependencies
- CUDA: The project uses CUDA (CUDA 12.x on the reference machine). Several executables depend on device code.
- HDF5 + HighFive headers: HDF5 libraries are linked directly; HighFive headers are included from an absolute path configured in CMake.
- YAML-CPP: Used for yaml configs (solverConfig.yaml). Linked in CMake.
- kdtree: External KD-tree library used in pre-processing (wall distance).
- METIS: Required for the mesh partitioning utility (mesh_part target) only.
- Paths are currently hard-coded in CMakeLists.txt:
  - hdf5_inc, hdf5_libdir
  - kdtree_dir
  Adjust these if your local paths differ. The reference machine paths are:
    hdf5_inc: /home/kumpei/app/HighFive/include, /home/kumpei/app/hdf/HDF5-1.14.3-Linux/HDF_Group/HDF5/1.14.3/include
    hdf5_libdir: /home/kumpei/app/hdf/HDF5-1.14.3-Linux/HDF_Group/HDF5/1.14.3/lib
    kdtree_dir: /home/kumpei/app/kdtree/install/

CMake Targets
- Executables:
  - forge: Main CFD solver (CUDA-heavy).
  - convertGmshToForge: Mesh conversion utility (Gmsh → Forge HDF5) and pre-processing.
  - mesh_part: Mesh partitioning using METIS.
- Libraries:
  - cuda_forge: Device-side kernels and wrappers.
  - probe: Probe functionality.

Build and Configuration Instructions (CLion profile)
- Use the existing CLion-configured CMake profile (Debug) and its build directory:
  - Profile: Debug
  - Build dir: /home/kumpei/work/forge/solver_density_cuda/cmake-build-debug
- Do NOT create new build directories; use the existing profile.
- To build a specific target from terminal within CLion’s environment:
  cmake --build /home/kumpei/work/forge/solver_density_cuda/cmake-build-debug --target forge -j 6
- Notes:
  - CUDA architectures: Not explicitly set. If you need to target a specific GPU arch, set CMAKE_CUDA_ARCHITECTURES at configure-time, or add an explicit set in CMakeLists.txt (prefer local-only changes or presets). The reference build produced sm_52 code paths.
  - OpenMP is enabled for host code via -fopenmp on the forge target.
  - If you change dependency locations (HDF5, HighFive, kdtree, METIS), update the include_directories/link_directories in the top-level CMakeLists.

Runtime Configuration
- The solver (forge) expects a solverConfig.yaml in the working directory and mesh data in HDF5 format.
- Key config fields consumed early:
  - cfg.meshFormat: must be "hdf5" currently.
  - cfg.meshFileName: HDF5 mesh file path.
  - cfg.valueFileName: initial field values (HDF5).
  - cfg.nStep, cfg.nStage, cfg.gpu, etc.
- Boundary conditions are read via readBcondConfig, and probe points are handled by probe library.
- Output: HDF5 + XDMF files via outputH5_XDMF.

Testing: How to Add and Run a Simple Test (Ad-hoc pattern)
- There is no built-in unit test framework. Use the ad-hoc executable pattern for quick checks. Example procedure:
  1) Create a temporary test source that exercises a small component without external runtime deps (avoid GPU/IO for speed). Example content:

     // File: tests/vector_util_smoke.cpp (temporary)
     #include <iostream>
     #include <vector>
     #include <cassert>
     #include <algorithm>   // for std::find used by common/vectorUtil.hpp
     #include "common/vectorUtil.hpp"
     int main() {
         std::vector<int> a{1,2,3}, b{3,2,1}, c{1,2,4};
         assert(ifEqualComponent(a,b));
         assert(!ifEqualComponent(a,c));
         assert(findIndex(a,2) == 1);
         assert(findIndex(a,5) == -1);
         std::cout << "OK" << std::endl; return 0;
     }

  2) Add a temporary target in CMakeLists.txt (do not commit long-term):

     # temporary test target
     add_executable(vector_util_smoke tests/vector_util_smoke.cpp)
     set_target_properties(vector_util_smoke PROPERTIES LINKER_LANGUAGE CXX)

  3) Build and run in one command using the existing Debug profile:

     cmake --build /home/kumpei/work/forge/solver_density_cuda/cmake-build-debug --target vector_util_smoke -j 6 && \
     /home/kumpei/work/forge/solver_density_cuda/cmake-build-debug/vector_util_smoke

  4) Clean up: remove the test file and the temporary target changes in CMakeLists.txt before committing.

- Important: Keep such tests self-contained and CPU-only when possible to avoid GPU/HDF5 dependencies in quick checks.
- If you must test CUDA kernels, prefer wrapping device calls behind small host functions and test their invariants with small device allocations. Consider adding CHECK_LAST_CUDA_ERROR() after kernel launches during development.

Demonstrating the Test Flow (verified)
- The above ad-hoc vector utility smoke test was executed successfully on this environment by:
  - Adding the test source and target temporarily.
  - Building and running with the Debug profile in a single command.
  - Removing the temporary files/changes afterward.

Guidelines for Adding More Substantial Tests
- For multi-file tests, group them under tests/ and gate them with a CMake option (e.g., -DForge_DEV_TESTS=ON) so they don’t pollute production builds. Since we avoid adding new presets or directories here, keep such options local to your machine.
- Prefer testing pure host logic or small wrappers around device code to keep iteration fast.
- For integration runs of forge, prepare a small HDF5 mesh and a minimal solverConfig.yaml; run a few steps (nStep small, e.g., 1-2) and verify output files are produced. Do not check large data into the repo.

Debugging and Development Notes
- Error handling:
  - The macro CHECK_LAST_CUDA_ERROR (main.cpp) is useful to catch kernel launch errors early; consider adding it after kernel launches during debugging.
- Performance: 
  - Inner loops run on device; watch register pressure and local memory in CUDA kernels. Use nvcc --ptxas-options=-v or Nsight Compute for diagnostics.
  - Host-side uses -O3 and OpenMP; ensure thread-safe access to shared structures.
- Mesh and BCs:
  - Mesh is expected in HDF5 format; convert from Gmsh using convertGmshToForge if needed.
  - Boundary conditions are configured via YAML; ensure names match those handled in boundaryCond and device counterparts in cuda_forge/boundaryCond_d.*
- Probes and Output:
  - Probe initialization is done early; ensure probe configuration is consistent with the mesh.
  - Outputs go to HDF5/XDMF; verify filesystem performance when writing frequently.
- Code style / conventions:
  - Mixed C++/CUDA codebase; headers often assume transitive includes (e.g., std::find usage in vectorUtil.hpp). When writing tests or new TU’s, include missing standard headers explicitly to avoid ADL-dependent compilation.
  - Prefer explicit includes in new code to reduce ODR and build fragility.

Common Pitfalls
- Changing machines breaks hard-coded HDF5/kdtree include/lib paths → update CMakeLists accordingly.
- Running forge without solverConfig.yaml or with meshFormat != hdf5 will fail early.
- Building mesh_part requires METIS; if not installed, limit your builds to forge or convertGmshToForge targets.

Reproducibility
- Verified on the provided Debug profile:
  - forge target builds successfully.
  - An ad-hoc CPU-only smoke test target builds and runs successfully when added temporarily as described above.

Housekeeping
- Do not commit temporary tests; delete them and revert CMake changes after verifying behavior. This document (.junie/guidelines.md) is the only persistent artifact added by this process.

Troubleshooting: Junieのチャット欄でコピペ（コピー/貼り付け）ができないとき
- まず試すこと（環境依存せず有効）
  - テキストをドラッグで選択し、右クリックメニューから Copy を選択（または Edit メニュー → Copy）。
  - フォーカスを明示的に与える：コピーしたいメッセージ内を一度クリックしてから Ctrl+C（macOS は Cmd+C）。
  - 連続メッセージの一部だけを取りたい場合は、ダブル/トリプルクリックで行/段落選択を活用。
  - クリップボード監視系の常駐アプリを一時停止（セキュリティソフトやクリップボード履歴ツールによってブロックされる場合があります）。

- IDE（JetBrains系/CLion）での典型要因と対処
  - Keymap 競合：Settings/Preferences → Keymap → Copy/Paste のショートカットを確認。Vim系プラグイン（IdeaVim）が有効だとノーマルモードでは Ctrl/Cmd+C が効きません。対処：
    - Insert/Visual モードで選択してから Cmd/Ctrl+C、または
    - IdeaVim で "+y（システムクリップボード）を使う、または
    - 一時的に IdeaVim を無効化。
  - ターミナルとショートカットの干渉：内蔵ターミナルがフォーカスを奪っていると Copy がターミナル用のショートカットに奪われることがあります。チャット欄をクリックしてから操作、または右クリックメニューを使用。
  - 右クリックが「Copy」表示にならない：選択がされていない可能性。必ず範囲選択してから再度右クリック。
  - クリップボードへのアクセスが OS 側で制限：Linux/Snap/Flatpak 版だと制限がかかる場合があります。JetBrains Toolbox 版や公式 tar.gz 版の CLion を推奨。

- Linux 特有（Wayland/X11）
  - Wayland でのクリップボード制限が原因のことがあります。対処例：
    - IDE を XWayland で起動（環境によっては既定）、または Wayland 対応の最新 CLion を使用。
    - wl-clipboard（wl-copy/wl-paste）等のツールを導入し、システム側のクリップボード動作を確認。
    - デスクトップ環境のキーボードショートカットと衝突していないか確認。
  - X11 環境では中クリック（PRIMARY 選択の貼り付け）が使える場合があります（選択→中クリック）。

- Windows/macOS
  - Windows：Clipchamp/PowerToys/サードパーティ常駐がショートカットを横取りしていないか確認。必要に応じて CLion を再起動または管理者権限の有無を切り替え。
  - macOS：システム設定 → キーボード → キーボードショートカットでアプリ固有の上書きがないか確認。IME 切り替え系のショートカットと衝突しがちです。

- 代替手段（どうしてもコピーできない場合）
  - チャットの右上メニュー（…）に「コピー」「全体をコピー」「エクスポート」がある場合はそれを使用。
  - テキスト選択後にメインメニューの Edit → Copy を使う（ショートカット無効化回避）。
  - 必要箇所を選択し、スクリーンショットで共有（緊急時の回避策）。

- それでも解決しない場合
  - どの OS/ディストリ/ディスプレイサーバ（Wayland/X11）/CLion バージョン/プラグイン（IdeaVim など）/導入形態（Toolbox, Snap, Flatpak, tar.gz）かをメモし、再現手順（どの領域で選択し、どの操作をしたか）を共有してください。環境に合わせた追加のワークアラウンドを提案します。
