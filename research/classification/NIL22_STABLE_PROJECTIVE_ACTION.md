# nil22 image の stable 二源：全磁気商と平面 projective Lie action

2026-11-07。主義務 E_multiple_hyperplane。
全十一源のうち未帰属の `stable_line_conic` と `stable_three_lines` を扱う。
以下は全六磁気係数を保持する正規化と幾何学的還元である。
既知成分への所属、新成分性、又は E の完了は結論していない。

## 1. 元の完全源と全磁気商

座標順を `(q1,q2,x,y)` とし、`p=(x,y)`、`D3=xy(x−y)` とする。
既存の必要十分表示は

```
Hess = S+Hessian(D3)/2,
M = −adj(Hess)+p*(V0−p/3)^T,
{qi,pj}=Mij, {x,y}=0,
{q1,q2}=R+V0_2*q1−V0_1*q2,
R=m0+m1*x+m2*y+m3*x²+m4*xy+m5*y².
```

二源はそれぞれ

| 源 | S | V0 |
| --- | --- | --- |
| line_conic | `[[3,−15],[−15,3]]` | `(3/2,−3/2)` |
| three_lines | `[[−6,3],[3,−6]]` | `(0,0)` |

である。`q'=q+L(x,y)`、L は任意 affine 二成分、とした全六変数の作用行列を
作った。rank は 4 と 3 で、以下が完全な cokernel 不変量となる。

```
line_conic:
 I4 = m4−2m1/3−8m3,
 I5 = m5−m1/12−m2/12−m3,
 Rnormal=I4*xy+I5*y².

three_lines:
 I3 = m3+m0/27+m1/9+2m2/9,
 I4 = m4−m0/27+2m1/9−2m2/9,
 I5 = m5+m0/27−2m1/9−m2/9,
 Rnormal=I3*x²+I4*xy+I5*y².
```

正規化する L は元の六 m の**多項式**である。パラメータを一つも反転せず、
全六括弧の pushforward 恒等式を確認した。従って不変量が零になる境界も
同じ二変数・三変数源に文字通り残る。単一の磁気代表を一般源と取り違えていない。
各不変量が独立な座標を与えることも、cokernel 行列の選んだ残余列部分が単位行列
となることで確認した。

## 2. 二 anchor は通常の projective linear field

M の各行を base 平面の vector field と読む。
二次部分はそれぞれ `−x(x∂x+y∂y)/3` と `−y(x∂x+y∂y)/3` である。
従って両 anchor は P2 の通常の projective linear vector field であり、
trace 零の 3×3 matrix A,B で表せる。
この表示の全係数を直接再構成し、`h=1` で元の二 anchor と一致することを確認した。

three_lines では `[A,B]=0` であり、次の三一次式で同時対角化される。

```
(L1,L2,L3)=(3h+x, −3h+y, −3h+x−y),
Rbase=[[3,1,0], [−3,0,1], [−3,1,−1]],
Rbase*A*Rbase⁻¹=diag(2,−1,−1),
Rbase*B*Rbase⁻¹=diag(1,−2,1).
```

両場は projective に独立なので、Cartan 部分代数そのものである。

line_conic では通常の matrix commutator が `K=[A,B]=3(A+B)/2` となる。
`K³=0`, `K²≠0`、`D=2A/3` は固有値 `−1,0,1` を持ち、`[D,K]=K` である。
従って三次元既約 sl2 表現を Borel 部分代数へ制限した形である。
Borel 表現自体を既約と呼んではいけない。その表現は不変 flag を持つ。

このため、残る E 認識の対象を全く任意の anchor ではなく、上の二つの
projective Lie action と二・三磁気変数へ縮められた。
Lie action の型だけから既知 Poisson 成分への所属は従わない。

## 3. 実因子配置と、単純一次因子を infinity にする全変換

実 anticanonical 五次式は定数倍を除き

```
line_conic: h²*L*Q,
 L=x−y−12h,
 Q=36h²+3h(x−y)+2xy;

three_lines: h²*(x+3h)*(y−3h)*(x−y−3h).
```

line_conic の Q は rank3 である。`Q|h=0=2xy` なので h は割線、
`Q|L=0=2(y+6h)²` なので L は接線で、接点は `[h:x:y]=[1:6:−6]`。
`h∩L=[0:1:1]` は Q の外にある。
three_lines では平面商の四直線 h,L1,L2,L3 は三本ずつ独立である。
重複しているのは h のみであり、三本の単純一次式は三次元を張る。

各単純一次因子を新しい H とする全四通りを定数 PGL5 変換で計算した。
新しい H 行は厳密に H で割れ、Euler-vertical 補正で H-Casimir 代表になる。
逆変換で元の全150係数を再構成し、全六磁気係数上の Jacobi も確認した。

全四 chart とも新 modular matrix は rank2、平方零、定数部はその像内である。
従って nil22-image 型のままで、SS20 への直接認識は得られなかった。
modular zero での定数 Pfaffian は、line_conic が `−1/8748`、
three_lines が順に `1/5184,1/5184,−1/5184` であり、いずれも全磁気係数に依らず非零。
この非退化性だけを SS20 の仮定と取り違えない。modular 型が異なる。

既存の「相異なる三 H と既約 Q」「相異なる五 H」「二既約 conic」全源定理は、
ここでの反復 h をそのままの式で許していない。
同じ次数和にまとめ直したり、非既約 quadric/cubic を作っただけでは、
それらの全源閉包への包含は証明されない。

## 4. 同じ anchor ansatz だけでの因子分離の限界

二つの通常の projective vector field を affine 平面へ書けば、二次部分は
ともに radial field の一次式倍である。その二つの二次部分は平行なので、
anchor determinant の次数四部分は恒等的に零となる。
従ってこの ansatz を保つ限り、同じ infinity h に対する実五次式は必ず h² で割れる。
これは任意の次数の jet 実験ではなく、全 projective-action ansatz に対する次数恒等式である。

従って反復 h を**被約に分離する**実族を作るなら、その族はこの固定 h に対する
projective-anchor ansatz の外へ出る必要がある。
反復因子を保ったまま既知成分へ認識する別の実族まで除いたわけではない。
未完了の仕事は、この全二・三磁気変数源を保つ実変形又は別の既知全源認識である。
以前停止した一つの T01 frame とその非適合を、全 T01 非所属へ拡張しない。

## 5. 再利用可能な保存物

- `code/classification_nil22_stable_infinity.py`
- `research/classification/jobs/nil22_stable_infinity.json`
- `reports/CLASSIFICATION_NIL22_STABLE_INFINITY_20261107.json`（0.960 秒）
- `code/classification_nil22_stable_projective_action.py`
- `research/classification/jobs/nil22_stable_projective_action.json`
- `reports/CLASSIFICATION_NIL22_STABLE_PROJECTIVE_ACTION_V2_20261107.json`（0.373 秒）

両 mode は実行前に finite job gate を通した。V1 の projective-action 報告も不変で残し、
V2 では「Borel の既約表現」という不正確な説明ラベルのみを上記の sl2 制限の記述へ直した。
matrix、全係数、正規化、結論の数学は V1 と V2 で同じである。
producer を import しない別作者 replay はまだ未実施であり、成分数は増えない。

## 6. 全不変量上での T39/T40 認識と、three_lines の被覆済み境界

上の多項式正規化後、定数非退化値を固定 Darboux frame に送った。
frame、S61 perspective と全61座標の分母は全て非零定数である。
二・三不変量を全て未定元として、全150係数の復元、全134 S61 方程式、
全350 Jacobi 方程式、T39 と T40 の各45 graph を一度ずつ検査した。
`classification_nil22_stable_graph.py`、同名 job と
`CLASSIFICATION_NIL22_STABLE_GRAPH_20261107.json` に全式を保存した（4.005 秒）。

line_conic では、全二不変量に依らない u3 graph 残差が
T39 で `1/2`、T40 で `12` となる。従ってこの全源は両 graph の外である。
three_lines の T39 残差にも定数 `−7/4` がある。

three_lines について `Rnormal=k0*x²+k1*xy+k2*y²` と書くと、
T40 の45残差のうち非零なのは厳密に次の三本だけである。

```
9*(4*k0+2*k1+k2),  27*(k1+k2),  81*k2/4.
```

従ってその graph membership locus は `k0=k1=k2=0` と等しい。
元の全六磁気源では、上記 `I3=I4=I5=0` の三次元閉部分源が該当する。
正規化はその全パラメータで可逆なので、この**境界全源は T40 に入る**。
それ以外の three_lines 磁気値は T40 graph の外である。

これは境界 locus の全支持認識であり、three_lines 全源を T40 へ入れたという結論ではない。
この境界に T40 のみが通るとも主張しない。
上の負の graph 結果を、他の MPS 成分又は既知非MPS成分への非所属へ拡張しない。
