# 実 native unit adapter の二次係数が不要となる source 直和

前の一般的な基底変更の自然性より強い、元の比較 source 固有の消去がある。
新しい N2、gauge2、physical jet を求めずに証明できる。

元768比較列の最初640本は `(B48,P8)`、残る128本は
`(Res4(N^t)*Res4(S12_by32(G8)),0_8)` である。後者を `L*S` と書く。
L は48対48の形式単位である。同じ actual Cfree4 の monic と
`(1,tauA,tauB,tauA*tauB)` を使って Res4 を取ることが必要であり、異なる
source frame の行列を混ぜない。

原点では `G8=g*tauA*tauB` で、

```
g=(1/2,25/8,0,0,0,5,-125/2,0).
```

従って native S0 の28 Koszul 部分の QQ rank は7である。非零ベクトル
g による exterior contraction の rank7と、tauA*tauB による C4 上の
乗法の QQ rank1の積だからである。四つの free-zero 部分は別の native
四 target にあり、rank16を持つ。L0 は単位なので native origin 像の
rank は合計23となる。原点の N と逆行列の両積を独立に再生して確認する。

一方、凍結済みの実 X0 において旧640列の rank は4で、固定独立列は
`604,613,622,631`、全768列の rank は27である。したがってこの二つの
origin 像の和は **直和**である。単に全 source の定数部分の和が零という
だけでなく、元の定数 column 結合 c の全56成分について

```
B0*c_old=0,    (L0*S0)*c_native=0
```

が別々に従う。builder は27消去後の全741本、および一次 rank37の消去後の
全620本について、この二つを実 frozen origin 行列で直接確認する。

L0 の可逆性から `S0*c_native=0` であり、全次数の積の二次係数を展開すると

```
(L*S)_2*c_native
 = L0*S2*c_native + L1*S1*c_native + L2*S0*c_native
 = L0*S2*c_native + L1*S1*c_native.
```

よって **L2 全体**が消える。これは N2 を任意に零とする仮定でも、B を
同時に変えたという仮定でもない。元の source 行列自体の積の恒等式と
origin 像の直和を用いている。L2 に含まれる N2 および Res4 の二次補正は、
いずれもこの共通左因子の二次係数として消える。残る `L0*S2` 側に入る
monic、core、実 physical/source の二次変化はこの証明では消えない。

前の m1 消去定理の補正 `Xp1*a1` は既知の一次 source だけを使うので、
ここにも L2 は現れない。従って三つの未決 class の入力に必要な native
二次部分は、既知の L0,L1 と **実 S2** の source coefficient gate に
減る。L2 の求解を目的とした native second-gauge job は不要であるが、
S2 がすべて既知または零だと結論してはならない。

この proof/preflight は新しい二次係数、source solve、jet、rank検索、348
行列を計算しない。元の Q、B、N、旧正則商函数および全 source 列は変更
しない。Tの全長、原点深さ、nil 三乗の決定は引き続き別の gate である。

再現:

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_native_adapter2_cancellation_preflight_20261104_v63.py
```

60秒・2GB・32MiBの固定 source 検査。旧凍結 source は変更しない。
