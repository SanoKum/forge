"""ユーザ指定 (2026-08-19) の Wyslouzil Fig.3 ノズル上壁プロファイル。単位 mm、throat x=0、y は半高。"""
X_STR = -60.0          # 入口直管始点
X_CONV = -38.0         # 収縮開始
X_CUB = -21.272727272727272727
X_TH = 0.0
X_BLEND = 8.093175942
X_EXIT = 95.0

def y_user_mm(x):
    if x <= X_CONV:
        return 12.7
    if x < X_CUB:
        return 0.16 - 0.33 * x
    if x < X_TH:
        return 2.5 + 0.015512821 * x * x + 0.000243078 * x ** 3
    if x < X_BLEND:
        return 2.5 + 0.00194105 * x * x - 7.99459e-5 * x ** 3
    return 2.457620744 + 0.015709255 * x
