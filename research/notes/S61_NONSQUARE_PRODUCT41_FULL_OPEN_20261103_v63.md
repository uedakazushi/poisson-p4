# 元の六二次式の product rank は FULLU 全域で41

2026-11-03。元の三 core で

```
F=500a+90e²+3e−24,
B=(50a+18e²−9e−2)d+60e−40,
C=100a+18e²−9e−2,  J=5ad−6e−4
```

とし、**FULLU={adFBCJ!=0}** とする。元の六二次式 q0,...,q5 と、元の七 light 変数
`(v0,v1,v2,v3,v4,v6,v7)` から作る42本の三次式の係数行列は、標数0の FULLU の
**全幾何点で正確に階数41** である。有限例外を残す先行 packet は変更しない。

証明の依存は独立受理済みの次の二つである。

- `S61_NONSQUARE_PAIR_EMPTY_QQ_20261103_v63_REVIEW`。
  実804係数からの全 pair rank2、先行 gcd proportionality、binary reduction、有限Σの
  Nakayama 証明を合わせ、actual RAW24×56 の全点 row rank24 を与える。
  係数行列を写像 R56→R24 と読むと全域で全射であり、その転置の24本の三次式が各点で独立である。
- `S61_NONSQUARE_RAW24_PRODUCT41_BRIDGE_20261103_v63_REVIEW`。
  元六式全係数からの unit row change、実際の ell 座標、(Z,W) 商での四式独立性、
  17 kernel 方向、係数10aの非零 syzygy により、各点で RAW24 rank24 から元42行rank41を導く。

従って元の42行の rank41 欠損集合は空である。全域で階数が一定なので、その image は
84次元の pure cubic 自由束の rank41 subbundle、商は rank43 の局所自由束である。
ここで商の global な自由基底や、具体的な全開 right inverse を構成したとは主張しない。

全開の規約に注意する。新しい pair packet の P は `P=ad` である。
`B0=(50a+18e²−9e−2)P+(60e−40)a=aB` とすると、上の FULLU は
`a*P*F*B0*C*J!=0` と一致する。d=0 は含まない。J=0 は先行 product 欠損有限定理の
範囲であり、この全開定理をそこへ延長しない。

全 cubic row module の57基底、heavy projection rank7、および今回の product41 を
単純に加えて full57 を主張しない。heavy7 との和から直接得る下界は48であり、残る実9×43
criterion は `S61_NONSQUARE_CUBIC_NINE_ROW_CRITERION_20261103_v63` に固定した。
全三次階数57、四次全射、長65、ambient t 冪零性、全 Poisson 成分の尽くしは未解決の別 gate である。
実行済みの constructive pivot 数を増やしたとも主張しない。

検査入口 `code/s61_nonsquare_product41_full_open_20261103_v63.py` は両独立受理 report の
scope と source/evidence hash を検査する。新しい GB、有限体、generic-field solver は不要である。
