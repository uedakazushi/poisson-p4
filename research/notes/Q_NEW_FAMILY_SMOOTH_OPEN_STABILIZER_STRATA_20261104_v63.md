# 滑らかなuvp領域全体のactual stabilizer

2026-11-04。p=32u+v、領域は `u*v*p!=0`。一般点の既受理packetは変更しない。
新しい点・prime・GB・同値solver探索は使わず、全parameterの係数恒等式を扱う。

## 一般点のconstant-pencil論証はそのまま全域へ延ばせない

元の全Pからconstant bivector六係数に対するSchouten行列20×6を再構築した。
lex順の係数をkxy,kxz,kxw,kyz,kyw,kzwとすると、v≠0でそのkernelは正確に

    kxy=kxw=kyw=0, kxz=−u*kzw,
    (u²−64)*kzw=0, kyzは自由。

従ってu=±8ではcompatible constant空間が二次元へ跳ぶ。そのbasisは
`∂y∧∂z` と `−u∂x∧∂z+∂z∧∂w=∂z∧(u∂x+∂w)`。
この跳びを隠してgenericなCの一意性を全域で使ってはいけない。

一方、HとQ|Hのvertexはuvp領域全体でintrinsic。
Hを保つaffine変換Mがprojective tensorをscalar κで保てば、divergence covarianceは
homogeneous modular行列について `N M=κ M N`。従ってImNを保つ。
Nのrows(y,z),columns(h,x)のminorは4u²/3なので、全u≠0で
`ImN=span(ey,ez)`。これでCの一意性に依存せず、同じ全affine旗を得る。

## 全二次case

既受理の一般affine旗からの必要条件はu,v,Dが非零であることだけを使うので、
この領域全体で有効。D=a²−64c²、κ=1/D、f=D、j=D²、
e=(D−a)/v、g=−64c/v。三つの二次分子Q1,Q2,Q3は既受理表示を使う。

* u²≠256なら `ac=0`。さらにv+16u≠0ならc=0,a=±1。
* u=16では `Q2=(a+8c)²(1−(a−8c)²)`、
  `Q3+vQ2=16c(v+256)(a+8c)`。D≠0だからv≠−256なら同じくc=0,a=±1。
  u=−16も符号を逆にした同じ証明が使える。
* v=−16u、u²≠256では、a=±1,c=0に加えてa=0,c=±1/8だけ。
* u=16s、v=−256s（s=±1）では、二次条件の全解は
  `a=(t+epsilon)/2, c=s*(t−epsilon)/16`、t∈Gm、epsilon=±1。
  このときD=t*epsilon。Q1,Q3はQ2の倍となり、他の二次caseはない。
  ただしこれはまだ全tensorを保つ条件ではない。

## 残る全係数と群構造

Tを自由なy平行移動として、残るH,K,Zの解を固定された三つのconstant係数
Pxz,Pzw,Pyzから求める。その3×3 determinantは **−256u/(243D)**。
全領域でunitなので、余分なrank例外はない。
J=0の追加二familyについて、全150 homogeneous係数がexact係数環で零となる。
その一意なH,K,Zと全行列を保存した。CB等のsourceではなく、元のfull Poisson tensor
そのものへの有限作用である。

二つの共鳴点の四つのGm候補では、同じunit系でH,K,Zを解いても、全tensorの
**一次係数が追加条件t²=1を強制する**。例えばu=16,v=−256,epsilon=1では

    residual Pxz=(t²−1)(−h*x/192+h*w/24).

残るPyz,Pzwのresidualもt²−1の定数倍であり、各符号の六個の全residual係数を保存した。
四caseのそれぞれに非零有理unitの係数があるので、t²=1は必要十分。
t=±1は上記の四つの離散blockへ戻る。従って**extra Gmは存在しない**。
二次条件だけからのGm仮説は棄却され、最初の実際の非零residualも別fileとして保持する。

J=0非共鳴部でa=0,c=1/8,T=0の行列をSとする。
S²=E²=I、SE=ES、SN=−NS、EN=NEを全係数で確認した。
全四componentはGa,GaE,GaS,GaESで尽くされる。

同じJ=0二familyの恒等式にはu²−256の分母が一切ないので、共鳴点にもそのまま適用する。
全一次係数を検査した後に、共鳴点に別のconnected componentが増えることはない。

従って標数零でのactual algebraic projective stabilizerは正確に

| 領域 | stabilizer |
|---|---|
| v+16u≠0 | Ga × μ2 |
| v=−16u（u²=256を含む） | (Ga ⋊ μ2) × μ2（最初のμ2はGaを反転） |

QQ上の各stratumでもsplitし、Cおよび幾何的base extensionでも同じ。
parameter間のKlein four同定と、同じtensorのstabilizerを混同しない。
u=0、v=0、p=0の境界は未分類。専用の独立監査前にledgerへは採用しない。

再現：`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_smooth_open_stabilizer_strata_20261104_v63.py`
（60秒/2GB上限）。
