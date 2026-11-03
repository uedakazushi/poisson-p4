# B18 odd一次行列の座標表記：raw7とadjusted6の区別

2026-11-03。独立reviewから、前稿の「Y28方向rank6」という略記の座標を明示する必要が指摘された。
前稿の凍結済みsource・証明書・noteは変更せず、本追記で明確にする。

CA/CBの共通raw29座標は(Y,row61)である。
前稿で定義したadjusted座標は(Y,eta)、eta=row61−L*Yであり、Lは族ごとに異なる。

* actual odd一次行列の全29方向rankは両族とも7。
* **raw row61を固定**したraw Y28列のrankも両族とも7。
* **adjusted etaを固定**したY28列、すなわち前稿の保存`tangent`行列のrankは6。
* 座標によらず表すと、actual一次行列をdiscriminant微分dnuのkernelという28次元接超平面へ制限したrankが6。

actual行列へraw←adjusted座標の微分を掛けると、保存tangent列＋normal列と全係数で一致する。
従って数値の矛盾ではなく、固定する29番目の座標の違いである。
前稿のrank6を「raw row61一定のY28」と解釈してはならない。

Nでもactual odd全28方向rankは7、adjusted27方向・intrinsic discriminant tangent制限のrankは6である。
Nのextra61 odd係数、dnu_Nとの独立性、共通R29へのprojected Nのordinary-pinch完成型は影響を受けない。
以前のmodulo discriminantでの六つの一次消去もintrinsic接超平面に関する結果なので変わらない。

* code: `code/ambient_B18_odd_firstjet_coordinate_scope_20261103_v63.py`
* receipt: `reports/AMBIENT_B18_ODD_FIRSTJET_COORDINATE_SCOPE_20261103_v63.json`
