# W≠0 の実際の交叉を既存218式へ明示的に認識する

2026-11-03。個々の Ac8/Ac7 型の一致から交叉を推測しない。W=1 の実際の
テンソルを一つの GL5 行列で移し、既存の全 target image ideals と照合する。

## 1. 有理的な frame と literal な旧 source

元の座標を (h,x,y,z,w) とし、新座標を

`(x0,x1,x2,x3,x4)=(x,−h,y,z−x,w−h)`

とする。この行列は行列式1で、rank4 の点 (0:1:0:1:0) を[e0]へ移す。
全150物理係数の厳密な恒等式により、像は

`x0²J−x0*x1 ∂y∧∂z+(2x0*x3+x1²+x3²) ∂y∧∂w`。

これは旧 universal critical-β source の
`p=s=r=0, k0=1,k5=1/3, その他のk=0` に文字どおり一致する。potential は
`H=x0*(x1*x4+x3²)+x1²*x3+x3³/3` である。S61 normalizer に通して
全61座標、134関係、150再構成を検査する。primitive β24 も旧 critical
rank2 の基点 B=C=1 に完全一致し、Sp4 orbit 探索や型による類推は不要。

## 2. 旧 finite 表示をこの点へ局所化できる理由

旧49変数の target 座標では、対応点は **v4=v7=1、他の全変数0** となる。
CA、CB の10個の実際の q 座標の fiber 行列式はこの点でも−72。N の9個の
fiber block も可逆である。β の7方向の orbit 単位小行列は β だけに依存し、
同じ critical β なので不変である。従って旧全 source、実際の image kernel、
G7 形式積を、この点で再び使える。二つの正規化の tensor product に置換
してはいない。

旧109+109関係の全218式をこの target 点へ代入して零と確認し、さらに
新しい218×49数値 Jacobian を係数から作る。階数は CA28、CB28、和29。
旧 F を加えた N の階数は30である。これは critical slice で接次元
21、21、20、19、全錐で38、38、37、36を与え、別に行った W1 の実際の
kernel jet 計算とも一致する。

重要なのは、旧 M2 検証が多項式イデアルの大域的等式

`F²∈ID, ID+(F)=IN, ID:F=Q_emb`

を証明している点である。Q_emb は v1,v4,v5,v7 以外の45本の独立線形式。
新しい点はその上にあるので、F の非零加群は局所化しても消えない。
v4,v7 の平行移動はこの線形 prime 自体を変えない。これらの受理済み等式を
再局所化するのであり、接空間の余剰だけから nilradical を推定しない。

## 3. 全完成環と深さ

証明書に保存した218式で v4,v7 をそれぞれ v4+1,v7+1 に置換して原点へ
移す。この49変数の商の完備化に、旧 β-orbit の7方向と symplectic-value/
perspective の受理済み10方向、合計17個の smooth 変数を加えたものが
W1 の **実際の全交叉 D の完成環**である。

この点でも actual odd ideal の regular sequence は9個、N 側は8個であり、
free な k0,k2,k5 の平行移動で失われない。従って N の critical slice は
Ac7 に3個の smooth 方向を加えたもの、深さ4である。nilradical は主な
平方零 ideal (F)、その全加群は critical slice で正則4次元加群となる。
完全列とこの associated prime から D の深さも4。17方向を戻すと

* 錐 D: `(dim,edim,depth,e)=(28,37,21,2)`。
* 射影 D: `(27,36,20,2)`。
* nilradical: 正確な冪零次数2、全加群 `C[[u1,...,u21]]F`。
* associated primes: 被約 N の極小 prime と、滑らかな21次元の embedded
  prime の正確に二つ（射影では embedded 層の次元20）。

環の直和分裂は主張しない。全非零 W が同じ projective orbit にあるという
受理済みの群恒等式は CA、CB、N と実スキーム的交叉を保つから、この全
完成型は W≠0 の全てへ移る。W=0 の元の B18 原点へは拡張しない。
新しい大域 prime や Poisson 成分を数えた結果でもない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_W1_S61_intersection_recognition_20261103_v63.py --output work/W1-S61-recognition
```

60秒・2GB上限。新しい source jet、M2、Gröbner 基底は実行しない。
