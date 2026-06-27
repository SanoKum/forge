# architecture-axisym-nozzle-geometry

## メタ

- **area**: `architecture`
- **status**: `done`
- **related_docs**:
  - [`methods/axisymmetric/theory.md`](../../methods/axisymmetric/theory.md)
  - [`methods/axisymmetric/implementation.md`](../../methods/axisymmetric/implementation.md)
- **related_plans**:
  - [`architecture-axisymmetric.md`](architecture-axisymmetric.md)
- **created**: `2026-05-31`
- **owner**: `Copilot`

## 1. 目的

`case/23.axi_nozzle/mesh_axisym_m2/` の検証ノズルは、収束部は滑らかでも、
発散部が sharp-throat minimum-length **planar** MOC を土台にしていたため、
軸対称ノズルとしては出口半径の相似補正に頼っていた。これを
source-corrected compatibility を使う **axisymmetric MOC** と、喉部直後の
5 次多項式ブレンドへ更新し、軸対称ノズル向けの geometry generator を持つ。

## 2. スコープ

- **やる**:
  - `mesh_axisym_m2/generate_moc_contour.py` に axisymmetric MOC mode を追加する。
  - `mesh_axisym_m2/generate_moc_contour.py` に喉部直後の 5 次ブレンドを追加する。
  - 既存の sharp-throat planar MOC は比較用として残し、axisymmetric mode を選べるようにする。
  - `mesh_axisym_m2/` に true 2D mesh (`gmsh -2`) の `.geo` / build script を追加する。
  - geometry diagnostics に接続点・曲率情報を出す。
  - `case/23.axi_nozzle` に public reference contour として、収束部も含めて公開式だけで閉じた bell nozzle 比較用 mesh path を追加する。
  - 関連 docs を更新する。
- **やらない**:
  - solver 本体の物理モデル変更。
  - 他ケースの既存ノズル生成スクリプトの全面統一。

## 3. 関連 docs と前提

- 軸対称 solver 側の理論・実装は [`methods/axisymmetric/`](../../methods/axisymmetric/) を参照。
- 本 plan は validation nozzle の形状生成改善であり、solver core の離散化変更は含まない。
- 現行 generator は「Witoszynski 収束部 + sharp-throat planar MOC + 軸対称出口半径への相似補正」。
- public comparison path では、Witoszynski を流用せず「直線 chamber + conical convergent + 1.5Rt throat entrant arc + Rao bell divergent + straight test section」で完全に公開式だけの contour を使う。
- 追加する axisymmetric MOC は、planar invariants をそのまま使うのではなく、
  characteristic 上の compatibility に一次の source correction を入れる近似実装とする。

## 4. 設計方針

- 喉部位置を `x=0`、喉部半径を `r_t` とする。
- axisymmetric MOC では characteristic invariant の代わりに、
  `K_- = theta + nu`, `K_+ = theta - nu` に対して
  `dK_-/ds = +S`, `dK_+/ds = -S`, `S = sin(mu) sin(theta) / r` の
  一次積分を使う。
- interior / axis / wall 点は、planar MOC と同じ幾何交点構築を使いながら、
  上記 source correction を fixed-point で 2-3 回更新して決める。
- 収束部終端 (`x=0`) では Witoszynski 曲線の `y, y', y''` を境界条件に使う。
- 発散部側の接続点は、sharp corner そのものではなく、少し下流の滑らかな MOC 点を使う。
- 区間 `[0, x_1]` を 5 次多項式で埋め、両端の `y, y', y''` を一致させる。
- `x_1` は MOC 壁点列の index で指定し、初期値は先頭数点を置換する保守的な設定とする。

## 5. 実装ステップ

1. `methods/axisymmetric/theory.md`, `methods/axisymmetric/implementation.md` に検証ノズル形状の位置づけを追記する。
2. `.github/plans/README.md` に本 plan を追加する。
3. `case/23.axi_nozzle/mesh_axisym_m2/generate_moc_contour.py` に planar / axisymmetric MOC 切替と source-corrected compatibility を追加する。
4. 同 generator に 5 次ブレンド関数と CLI 引数を追加する。
5. `case/23.axi_nozzle/mesh_axisym_m2/build_mesh.sh` の既定を axisymmetric MOC に切り替える。
6. `generate_moc_contour.py` 実行で contour/geo/preview が正常生成されることを確認する。
7. public Rao contour generator で contour/geo/preview と fresh run が正常生成されることを確認する。
8. public Rao contour が `--convergent-model public` 既定で、既存 Witoszynski も比較用オプションとして保持されることを確認する。

## 6. 検証

- **単体 / ビルド**: `python3 generate_moc_contour.py --design-model axisymmetric --moc-model axisymmetric` が完走する。
- **検証ケース**: `case/23.axi_nozzle/mesh_axisym_m2/`。
- **判定基準**:
  - 喉部から発散部への曲線が単調に接続される。
  - `contour_preview.png` で喉部直後の折れが消える。
  - axisymmetric MOC mode が planar mode と異なる発散壁を返す。
  - `characteristics_preview.png` で C+ / C- の特性網が確認できる。
  - true 2D mesh (`gmsh -2`) が生成でき、fresh run で solver が完走する。
  - `axi_nozzle_points.geo` と CSV が正常出力される。

## 7. 影響範囲

- 触るファイル:
  - [`methods/axisymmetric/theory.md`](../../methods/axisymmetric/theory.md)
  - [`methods/axisymmetric/implementation.md`](../../methods/axisymmetric/implementation.md)
  - [`case/23.axi_nozzle/mesh_axisym_m2/generate_moc_contour.py`](../../case/23.axi_nozzle/mesh_axisym_m2/generate_moc_contour.py)
  - [`case/23.axi_nozzle/mesh_axisym_m2/build_mesh.sh`](../../case/23.axi_nozzle/mesh_axisym_m2/build_mesh.sh)
  - [`.github/plans/README.md`](../README.md)
- 既存 solver 実行手順への影響はない。
- validation nozzle の壁形状と mesh preview は変わる。
- public reference contour の比較用 mesh/run が追加される。
- public reference path は上流側も既存 in-house contour から切り離される。

## 8. 完了条件

- [x] 関連 `methods/axisymmetric/theory.md` 更新済み
- [x] 関連 `methods/axisymmetric/implementation.md` 更新済み
- [x] 実装・検証完了 (本 plan の §6 を満たす)
- [x] `.github/plans/README.md` の状態を `done` に更新
- [x] 本 plan の `status` を `done` に変更し、§9 に変更ログを記載

## 9. 変更ログ

- `2026-05-31` — 初稿。
- `2026-05-31` — axisymmetric MOC + quintic blend 方針へ更新。
- `2026-05-31` — true 2D mesh path (`gmsh -2`) を追加。
- `2026-05-31` — characteristic net の可視化 PNG を追加。
- `2026-06-02` — public reference contour として Rao bell nozzle 比較 path を追加。
- `2026-06-02` — `run_0032_slau_axisym_m4_publicrao_10k` を fresh run で完走確認し、public Rao contour でも 10000 step / CFL 0.8 の安定計算を確認。
- `2026-06-02` — public reference contour を、Witoszynski 流用ではなく収束部も含めて公開式のみで閉じた bell nozzle path へ拡張。
- `2026-06-02` — 本 plan を完了扱いに更新。
  - fully public bell contour を fresh run `run_0033` で検証済み。
  - 同一メッシュの SU2 比較 `run_0034` でも出口面の radial nonuniformity が再現された。
  - `.github/plans/README.md` の状態を `done` に同期。
