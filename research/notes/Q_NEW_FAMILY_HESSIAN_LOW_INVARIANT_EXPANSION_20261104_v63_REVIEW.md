# 四つの追加 source word map：直接Hessianによる独立検証

2026-11-04。実際の1032-node回路からbeta16を再構成し、元の表現とcocharacter
を使ってT²=Bのpreconditionを照合する。source点は評価しない。

thetaは元のmoment定義、H16は受理済み全117項のquarticから
`Omega^-1 Hess(F4/2)` を直接計算する。元の全136 quadratic monomial上で
2thetaと2H16が整数係数であることを確認する。この入口は著者のrank10
vector分解と異なる。全entryのCartan weightを検査し、有限torusでの共変性と
全trace wordの不変性を束縛する。純trace6/8のcharacter恒等式も独立に計算
するが、下のworkerはその恒等式を計算の短縮に使わない。

preconditionは `bbar=T^-1 torus(T)b(B=T²)` である。beta次数dの不変式は
`I(bbar)=B^-d/2 I(b)`、従って実際のbetaへの値は
`I(beta)=I(bbar)/[(1080H)^d B^(d/2)]`。cを1にするのは記号的で、最終c冪は
d/2。全16座標で非負T指数0/1となること、元のnative inputの全係数を検査する。

独立integer workerは2thetaと2Hessianの**16個の列を直接**反復する。
rank10のZ recursionも、純character項への分割も使用しない。CD二wordも
別の積順序で算出する。全4つの保存多項式の各係数を、整数とTの冪を払った
完全なQQ恒等式として比較する。T偶性、u=aX後のX偶性も確認する。
term数3449、3618、3948、10819、合計21834個の全係数を対象とする。

この結果は追加四wordの正確なsource mapのみ。新しいbeta点、全67列の関係
検査、NEWのdegree24方程式は含めない。一つのworkerが上限に達した場合は
入力・途中logを保存し、自動retryや予算拡大を行わない。

```sh
. /workspace/cloud-setup/activate.sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_hessian_source_review_20261104_v63.py --output work/hessian-source-independent --work work/hessian-source-independent-native
```

全体60秒・2GB上限。direct workerのwall上限50秒。
