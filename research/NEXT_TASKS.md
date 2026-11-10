# 次の課題 — 非被約接錐の高次初期方程式

2026-11-10 / v16。再開はRESEARCH_POLICY.mdから。GitHub main正本、PDF非追跡、Actions不使用、Libraryは実push失敗時のみ。

証明済み：30次元接空間、全7本の二次初期関係、h=a(s2+c)による真の接錐の非被約性、唯一の最大次元接錐V(a,c)、局所次元28と重複度1、混合ac障害の明示的な周囲二次補正。新原稿 papers/isotropic_log_tangent_cone_ja.tex と自己完結producer check_isotropic_quadratic_cone.py。

新たに証明：唯一の28次元Log1112成分はpで滑らか、T_p Zlog=V(a,c)、完備局所環は28変数形式冪級数環。20次元PGL5軌道がLog成分と25次元long rank4成分の交叉に含まれ、交叉次元は23〜24。交叉の厳密値は未解決。gl5微分のrank21は新しいcheck_isotropic_orbit20.pyで独立検算。

共通23次元族を構成し、rank2側のLog留数極限とrank4側の非定数g長鎖全六括弧を確認。24次minor=10871635968、ゲージを含めたrank48。下界23自体は一般交叉次元不等式からも出るが、明示的源を得た。未解決は交叉が23か24か、及び高次scheme構造。新コードcheck_isotropic_intersection_dim23.pyと証明書。

新しいv15結果：Aを固定する5パラメータの滑らかなqではPλが普遍的にstrict Poissonで、Pfaffian=(3λ/2)u1²z dq。λ=0はLog1112閉包、λ≠0は長鎖非定数g成分。qの五係数はPGL5 shearでη'=η/(1−δη)一係数へ規格化される。したがって固定Aのqを複雑化するだけでは24次元方向は増えない。現行論文papers/isotropic_nilpotent_pencil.texに証明、自己完結producer・JSON・手動照合を追加。

v16：g_ρ=diag(1,1,ρ,1,ρ^-1)は(η,λ)↦(ρ²η,ρλ)と作用し、任意の固定A正規Poisson pencilを完全等方rank2点へ退化させる。二次形式はrank(q|kerψ)=0/1の二軌道のみで、λ≠0のη/λ²は一変数群の不変量。証明papers/isotropic_nilpotent_pencil.tex、コードcheck_nilpotent_torus_attraction.py、対応証明書を手動検証入口に統合。全交叉の次元23/24は未決定。

最優先：固定Aと[ψ]を保つ二軌道の外の共通変形、または三次以上の接錐障害で交叉23/24を判定し、微分の95 normal変数を二次まで解き、三次残差を二次関係のsyzygyと合わせて全初期イデアルへ進む。二次障害模型の26/24次元成分が真の接錐で何になるか、hのannihilatorと台がどこまで縮むかを決める。その後、既知25次元長鎖成分との交叉、元の局所環Rの被約性、全次数の変形を研究する。接錐非被約とR非被約、二次模型と真の接錐、混合acの二次消去と全次数の持上げを混同しない。

手動検証：bash scripts/verify_github_source.sh tangent-cone。
