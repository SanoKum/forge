"""ParaView 5.11 + Python 3.12 用の互換 shim (PYTHONPATH で読ませる)。

ParaView 5.11 の `paraview/detail/pythonalgorithm.py` は Python 3.11 で削除された
`inspect.getargspec` を import するため、Python 3.12 上では Python プラグイン
(Tools > Manage Plugins からの .py 読み込み・自動ロード) が全て失敗する。

このファイルを PYTHONPATH の先頭に置くと、インタプリタ起動時 (site の処理) に
`inspect.getargspec` が補われるので、ParaView 本体を書き換えずにプラグインが読める。
ParaView 側での用途は `getargspec(func).args` の参照のみなので `getfullargspec` で等価。

利用は tools/run_paraview_native.sh 経由が簡単:
    ./tools/run_paraview_native.sh case/.../res_1000.xmf
"""

import inspect

if not hasattr(inspect, "getargspec"):
    inspect.getargspec = inspect.getfullargspec

# 同名の sitecustomize が他にもある場合に備えて、後続のものも読み込んでおく
try:
    import os
    import sys

    _here = os.path.dirname(os.path.abspath(__file__))
    for _path in sys.path:
        try:
            _same = os.path.samefile(_path, _here)
        except OSError:
            _same = False
        if _same:
            continue
        _cand = os.path.join(_path, "sitecustomize.py")
        if os.path.isfile(_cand):
            with open(_cand) as _fp:
                exec(compile(_fp.read(), _cand, "exec"), {"__file__": _cand, "__name__": "sitecustomize"})
            break
except Exception:  # shim 本体の役目は果たしているので、ここでの失敗は無視する
    pass
