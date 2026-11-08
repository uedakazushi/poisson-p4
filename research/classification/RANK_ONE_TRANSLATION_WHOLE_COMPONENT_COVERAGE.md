# rank一・c非image：全long/short源の成分被覆

2026-11-08、PRIMARY `E_multiple_hyperplane`。純全source合成。
**producer証明。最後のM1/M11行は新しいglobal dimension inventoryのroot受理を待つ。**
このnoteは台帳、旧finite-routingの不変bytes、旧negative quadratic証明書を書き換えない。

## 1. 全入口と有限被覆

modular線形部rank一、定数部cがimageに入らない場合の全affine型はlongとshort。
longでは `MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §5 の

```
X=∂x+x∂y, r=y−x²/2, q=(z,w), d=div U,
U(g)=d*g, (U+1)H=F*g,
deg g≤1, deg U≤1, deg F≤2, deg H≤3,
Pxy=−F−(1+d)r, Pxz=u, Pxw=v,
Pyz=xu+r*g_w+H_w, Pyw=xv−r*g_z−H_z, Pzw=g
```

を全て使う。rは原式導出用であり、非線形座標をPGL所属の代わりに使わない。
`RANK_ONE_TRANSLATION_WHOLE_FINITE_EXHAUSTION.md` と独立FINITE_EXHAUSTION査読により、
全g/U/Hの零・非零、全scalar/Jordan/zero-eigenvalue/translationを保った十七whole源で尽くされる。
Fは全行で任意二次以下の六係数であり、constant Fだけへ縮めない。

## 2. 十七行の全帰属

| 全long source | 全パラメータを保持した閉成分への帰属 |
|---|---|
| g非定数 | long非定数gの新25成分 |
| g非零定数、H3零 | 旧positive-weight NEW25 |
| g非零定数、H3非零、μ=1 | long定数g・μ1の新25成分 |
| 同μ=1/3 | long定数g・μ13の新25成分 |
| g零、M1 | linear-Hの全25像そのものが全成分、global bound theoremによる |
| g零、M2 | 既nil31-third25、全actual permutation/scalar逆による |
| g零、M11 | product-Hの全25像そのものが全成分、global bound theoremによる |
| g零、M3 | 既nil31-half25、全actual permutation/scalar逆による |
| g零、M21 | 既nil31-DDH25、全factor変更とcritical λ=−1/5の閉性による |
| g零、B1 | full25 component support、true50cutの真正quadratic最大根基による |
| g零、B2 | full25 smooth unique成分、同一点D100と全source下界による |
| g零、B3 | full25 component support、true50cutの真正quadratic最大根基による |
| g零、B4 | full25 smooth unique成分、同一点D100と全source下界による |
| g零、B5 | full25 smooth unique成分、同一点D100と全source下界による |
| g零、B6 | full25 smooth unique成分、同じprojective pointの整数化24によるD100 |
| g零、V1 | full25 component support、true50cutの真正quadratic最大根基による |
| g零、SC | scalar binary cubicの新25成分 |

M1/M11について用いるのは
`RANK_FOUR_OUTSIDE_LOW_MODULAR_FLAG_DIMENSION_BOUND.md` §9 の全成分次元論である。
全25像の一般元は唯一重複hのcanonical modular行列がrank二なので、
proper閉旗「ある重複Hのmarked Nがrank一以下」の外にある。
その開集合の全成分は射影次元≤25。一方、両全像は正確25の既約閉集合なので、
それを含む全空間成分と一致する。新しい固定点や高次jetを使わない。

M2/M3/M21の正逆はfinite-routing §§4–5、他の全source/成分性は次の凍結証明を使う：

* `RANK_ONE_LONG_NONCONSTANT_25_COMPONENT.md` と別作者AMBIENT査読。
* `RANK_ONE_LONG_THREE_NEW_25_COMPONENTS.md` と別作者THREE_AMBIENT査読。
* `RANK_ONE_GZERO_FOUR_BINOMIAL_25_COMPONENTS.md` と別作者NINE_AMBIENT査読。
* `RANK_ONE_GZERO_THREE_QUADRATIC_SUPPORT_COMPONENTS.md` と別作者FIVE_QUADRATIC_SUPPORT査読。
* `RANK_ONE_LONG_LINEAR_AND_PRODUCT_ENTRANCE.md` と別作者LINEAR_PRODUCT入口査読。

各nonempty rank四開からの全source閉包への所属を使うので、係数零のboundaryも保持する。
finite-routingのunit係数が零なら該当する別のmonomial/zero sourceへ戻す。
g=H=0はPf零でrank四の外、非零constant Hは(U+1)H=0に反して空である。

## 3. 全short源の復帰

`RANK_ONE_SHORT_TO_LONG_WHOLE_POLYNOMIAL_CURVE.md` と独立
`RANK_ONE_LONG_ORBIT_AND_SHORT_CURVE_SEPARATE_REVIEW.md` のactual polynomial λ族を用いる。
λ零は任意の元short source、λ非零は真のaffine frameで同じg分岐のlong sourceとなる。
全F/Hの変換、六括弧の正逆、λ零の可除性はその証明に保持される。

従ってλ非零の全点は§2の有限閉和に入り、その閉性からλ零も入る。
λ非零におけるframeのpoleをλ零へ代入したり、shortの一点だけをlongへ移したと仮定しない。
全shortのconstant g/cubic resonance、g零、非定数gの全zero/unit complementが同じ議論で戻る。
これでrank一・c非imageの**全long/short rank四点**が有限既知閉和に入る。

## 4. M1/M11の一般軌道と主張の限界

`RANK_ONE_GZERO_VERTICAL_POISSON_FIELDS.md` により、全g零sourceでXはPoisson affine場。
さらにY_ell=ell(q)∂yがPoissonである必要十分条件は(U+1+d)ell=0であり、全F/Hを保持する。
従ってM1の全PGL軌道は23以下、M11はd=−1なのでXと∂yが常に独立で22以下。
既に保存された同じ九点ambient結果のM48/M47は下界23/22を与える。
既約全sourceと非零minor開により、一般軌道はそれぞれ23/22、reduced Rosenlicht商の次元は2/3。
この一般軌道結論に新しい点又はminorは使っていない。

M1/M11の既一点の接次元28/27、quadratic normal coneのplane／二lineは不変である。
両者とB1/B3/V1のcomponent support証明から、smooth、局所unique、正確なnil指数、
全非被約Poisson schemeの構造は従わない。他のsmooth-positive行だけの結論を表の全行へ拡張しない。
rankA一・c=image、rankA零、rankA二以上、reduced/reducibleの他因子型の点別被覆は本合成定理の射程外。
