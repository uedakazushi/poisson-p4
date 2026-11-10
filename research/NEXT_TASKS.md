# 次の課題 — 非被約接錐の高次初期方程式

2026-11-10 / v12。再開はRESEARCH_POLICY.mdから。GitHub main正本、PDF非追跡、Actions不使用、Libraryは実push失敗時のみ。

証明済み：30次元接空間、全7本の二次初期関係、h=a(s2+c)による真の接錐の非被約性、唯一の最大次元接錐V(a,c)、局所次元28と重複度1、混合ac障害の明示的な周囲二次補正。新原稿 papers/isotropic_log_tangent_cone_ja.tex と自己完結producer check_isotropic_quadratic_cone.py。

最優先：微分の95 normal変数を二次まで解き、三次残差を二次関係のsyzygyと合わせて全初期イデアルへ進む。二次障害模型の26/24次元成分が真の接錐で何になるか、hのannihilatorと台がどこまで縮むかを決める。その後、既知25次元長鎖成分との交叉、元の局所環Rの被約性、全次数の変形を研究する。接錐非被約とR非被約、二次模型と真の接錐、混合acの二次消去と全次数の持上げを混同しない。

手動検証：bash scripts/verify_github_source.sh tangent-cone。
