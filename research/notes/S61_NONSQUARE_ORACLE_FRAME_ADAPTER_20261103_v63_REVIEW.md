# oracle の二 frame と元の行の独立監査

2026-11-03。対象は `S61_NONSQUARE_ORACLE_FRAME_ADAPTER_20261103_v63`。
実際の C122 の rank や追加の点を計算せず、係数規約の一致を監査する。

## 1. 同じ元135式・同じ34行の graph

既受理の oracle レビューは元の134式と G6 から全135式を独立に再構成し、
同じ固定座標、44座標の順序、全 J/B、選択34行、T/U、重みを確認している。
今回もその入力 hash、旧 native source report の全入力 hash、保存 weights を確認する。
oracle の全 T/U atom は weighted chart に一致し、旧 source の初期 tangent は

```
Tn=T34*D,  D=diag(a,a,a,1,a,B/50,a,a),
Un=ud*U34,  ud=a²FB/25000.
```

D と ud は aFBC 開上の既知単元であり、D は重みを保つ。

## 2. なぜ genuine 行は一つの scalar 倍で結ばれるか

この点を anchor から推測してはいけない。元の二次式を Jx+B(x,x)、同じ選択34式の
graph の ordinary degree d の局所項を Z_d と書く。旧 native の recurrence を読むと
帰納的に

```
W_native,d = ud^(d−1) Z_d,
A_native,d = ud^(d−1) F_d(T34*D*v)
```

となる。d=1 は Tn の定義、d≥2 は二次式の分解次数の和から共通の ud^(d−2) が
出て、Un の ud が一つ加わることによる。残余 A も native の den*Q+JW から同じ指数。

extract は A2 に ud²、A3 に ud を掛ける。従って各**元の行ごと**に

```
A2*ud²+A3*ud=ud³*(F2+F3)(T34*D*v)
```

となる。weight3 行は ordinary degree4以上を持たないので、これがその全 residual である。
selected34 graph の一意性とも一致する。

extract の content 処理と Python の追加処理は既知単元による各行別の除算のみ。
元 full63 は unit-Schur の**入力前**に保存され、元 labels の最初42行は quadratic
products、残り21行は元方程式を一行ずつ取り出した genuine cubic である。
既受理の row-basis レビューは、残る57行がこの正規化済み63行の実際の部分集合で
あることを全成分で確認している。従って genuine16行を product41 の未知の組合せで
変えたものではない。

以上から各 genuine 行について global な一つの既知単元比が**先に**証明される。
今回の16個の非零 anchor はその scalar の値を特定する役割に限る。

## 3. 独立に比較した係数

レビュー source は author の Oracle.coefficient/core_value を用いず、凍結 DAG の
atom と四演算から、目的 monomial の係数を QQ(a,e,d) 上で直接再帰計算する。
六つの quadrics の全29個、計174係数（heavy 線形の零を含む）を比較する。
41 product 行の全91列、計3731成分を元 quadratic と保存63 normalization から再構成する。
16 genuine anchor と全57 source label、全単元因子を照合する。

全148858 genuine cubic core 項を再展開して照合したわけではない。
全行の一致は第2節の元 source に束縛された共変性と既受理の全63再構成に依存する。

## 4. quartic と H² の規約

旧 heavy-square source は Tn の heavy 列で作った native A2 を ud*D_H² で割る。
従って旧 b77,b131 は T34 frame の係数である。この事実は既存 oracle 独立レビューの
元 J/B からの直接 H 一変数計算とも一致する。

encoded frame の全33 genuine quartic には

```
Qhat_k(v)=F_k(D*v)/D_H²,  D_H=B/50
```

を使う。H² の座標倍率 D_H² が共通 row 単元で打ち消されるので、旧 b77,b131 と
それから作った u,v、単元 determinant は文字どおり保持される。全33行は元の
weight4 非選択行のままであり、新しい関係を追加したものではない。

独立係数比較の後だけ public API を呼び、quadratic 一係数、product 一係数、
genuine cubic 一係数、quartic H² 一係数の四 fixture を照合する。
query がメモリ内 DAG に node を追加することを数えるが、保存ファイルは変更しない。

以上の frame adapter を受理する。actual C122 rank64、全 quartic rank、全 cubic57 の
各点下界、ambient 除去、追加の点・体・GB は今回の結論ではない。

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/s61_nonsquare_oracle_frame_adapter_review_20261103_v63.py
```

60秒・2GB上限。全入力 hash、係数数、16単元、API結果、時間とRSSを REVIEW folder に保存する。
