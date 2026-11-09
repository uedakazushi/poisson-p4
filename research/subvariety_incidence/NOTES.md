# Poisson部分多様体incidence：現行研究ノート

2026-11-10 / v4。数学的証明の全文は`papers/poisson_subvariety_incidence_ja.tex`と`papers/quadric_log_extension.tex`。本書には現行の障害・具体例・次の着手方法を整理する。

## strict stabilizerとJacobiによる延長の原理

smooth Q=V(q)に対して全延長は標準のq-central lift s(σ)から`K_T=qT-(1/2)E∧v_T`を加えた10変数源。Poisson条件から`v_T|Q`はσのPoisson自己同型を与えるので、まず`aut(Q,σ)⊂so(q)`を線形計算し、残る小次元源のJacobi理想を解く。modular場はq因子の重複度を決める。Lie代数の次元が跳ぶことは非自明な延長の**必要条件に過ぎず**、十分条件ではない。

## 一般三成分と特殊Logの実例

- 一般Rat(1,2)：strict stabilizer 0、全延長はrank2の被約一点。非零minor 27980515917168640。
- 一般Log(1,1,1)：universal q_U・residue weightでrank2延長のA1。strict so(q) stabilizerはrotation Rのみ（rank9、minor 10240000/243）。全Jacobi/Pfaffianの多項式恒等式と一般開条件で証明。
- 一般Aff：被約二点。rank2とrank4 (D=2Q+H)。
- 特殊Log境界：σの1形式はx2(x1dx0−x0dx1)（定数倍無視）。strict stabilizer R,Sの2次元。全Jacobi理想は(v(t+2))、交叉二直線。rank4第二枝のPfaffianはQ*x2*(x2²+x3²+x4²)。

## 新しい直交性wallの意味

q_c=Σxi²+2c x0x2、c²≠1、共鳴重みxi=(1,1,-2)。
c≠0では重みの二次元固有空間とx2がquadricに関して直交せず、strict stabilizerはRだけ。q_cに適合する全so(q_c) generatorをG_c^−1 Zとして構成し、[Ax,P0,c]=0の73×10係数行列の9次minor
[
13824c/((c-1)^4(c+1)^4)
]
が全非零cで非零。したがって全延長はrank2のA1。c=0では直交性によってSが新対称性として現れ、全ファイバーがノード状二直線となる。**対数一形式の零因子x2はどちらにも存在する**ので、零因子そのものでは跳躍を検出できない。本文`thm:quadric-resonant-log-wall`と`check_log_wall.py`を参照。

## 残る障害と方向転換基準

今回のwallは一変数族の完全なファイバー計算であって、全特殊Logの一般必要十分条件ではない。一般Logの重み一致、直交配置、非退化quadricの幾何を不変な条件として整理することが次の課題。一般分類が冗長な場合は、核のLie代数分解や多項式延長理想に対する概念的な補題へ切り替える。その際、解決したこと（一般Rat/Log/Affとこのwall）を保持し、全特殊層は未解決と明示する。

旧del Pezzo一般剛性は20係数Jacobiのrank60／90 syzygyと二つの非零minorによる。特殊pencilや反標準切断での真の延長例外は未分類。P3の横断対称性の全6族への適用も未完。

## 正確な資料と保存

LPT: arXiv:1112.3871, DOI 10.1002/mana.201100354。Mayanskiy: arXiv:1303.5916, DOI 10.1016/j.jalgebra.2014.08.049。文献上の優先権は未確認。
GitHub mainの5つの`code/subvariety_incidence/*.py`と5つの`research/subvariety_incidence/certificates/*.json`、`scripts/reproduce_subvariety_incidence.sh`だけから現在の論文と有限証明書を再生成する。GitHub Actions禁止、Libraryは実push失敗時の緊急退避のみ。
