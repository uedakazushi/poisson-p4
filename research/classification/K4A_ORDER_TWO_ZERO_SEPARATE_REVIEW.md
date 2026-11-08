# K4a actual24 親の order二零点：全大域 operator の独立 review

2026-11-08、PRIMARY `E_reduced_nonnormal`。別作者による **全元係数 Fraction 再構成 PASS**。
K4a の accepted actual24 親は completed T36/T38 のどちらでもない。
元の全一軌道 K4a source も両成分の外にある。
**T37 所属、新二十四次元成分、他の既知親からの相違又は分類完了は証明していない。**

## 1. 全 source に対する依存と今回の有限計算

`K4A_COMPLETE_TORUS_SLICE_SEPARATE_REVIEW.md` と正式 V2 receipt の
完全有限 C* graph、smooth unique24、および全非零 parameter が一つの actual24 orbit
という独立受理済みの定理を入力とする。
元 K4a source の全正規化はその canonical23 orbit に戻り、canonical point はこの親の境界にある。
今回それらの D101、G48、normal、IFT graph 又は rank profile を再計算していない。

同じ literal graph の `t=174`、整数再parameter化の `e=1` を使った。
これは既受理全非零 orbit の一つの表示であり、新しい点・q chart・profile の選別ではない。
独立実装 `classification_k4a_order_two_zero_review.py` は standard-library の
`Fraction`、AST による有限 Laurent 多項式読取、疎な加減乗算と微分だけを使う。
producer、current helper、SymPy、CAS を import・実行しない。
active 有限 contract を先に gate PASS としてから一回実行した。

全百五十 coefficient order を ten skew pairs × fifteen degree二 monomials と照合し、
保存 graph の全 coefficient を `174^degree` 倍して actual tensor を再構成した。
これを producer の全十 actual 括弧と比較し、零 entry も含め完全一致した。

## 2. Projective tensor の intrinsic first jet

任意の斉次二次代表 P と Euler 場 E に対し、T=E∧P の全十 cubic entry を使う。
Euler gauge P→P+E∧L は T を変えない。
q_i が非零の affine chart で q_i=1 とすると、T の i を含む entry は
projective bivector の局所係数であり、残る entry はその局所係数の座標多項式結合である。
従って projective Π の値と一次 jet が零になることは、T の全局所 entry と
全 chart 方向の一次微分が零になることと同値である。

T は次数三なので Euler 公式 Σq_l∂_lT=3T により、
全 ambient first derivatives ∂_lT(q)=0 は T(q)=0 も含む。
逆に chart の値と一次 jet の零から、残る radial derivative も零になる。
従って五十式

`∂_l(E∧P)_ijk(q)=0`, `i<j<k`, `l=0,…,4`

は Π の intrinsic な order≥2 零点条件そのものである。
代表変更と PGL に不変であり、chosen affine chart 上だけの条件ではない。

## 3. 全五十式・全十五 monomial と非零 minor

上記 same-orbit actual P から全十 cubic を定義通り作り、全五十微分を独立に展開した。
全十五 monomial の順序は `q_i q_j`, `0≤i≤j≤4`。
全五十式、各 equation order、50×15 の全二百三十四非零 entry が producer と完全一致した。

保存された row

`0,1,2,4,5,6,7,9,10,11,12,15,16,17,20`

と全十五 column の square minor を QQ Fraction elimination で計算し、

`det = 3411634563591564165120 ≠ 0`

を得た。行列は十五 column なので rank は正確に十五、全 kernel は零である。
実行は **0.01753秒、RSS12,620 KiB**。
正式 receipt は `CLASSIFICATION_K4A_ORDER_TWO_ZERO_SEPARATE_REVIEW_20261108.json`。

もし q が complex projective solution なら非零 vector `(q_i q_j)` が kernel に入る。
少なくとも一つの q_i² は非零なので矛盾する。
同値に、この五十式は degree二形式の全空間を張り、homogeneous ideal は irrelevant ideal の二乗を含む。
従って全 P4 上の projective zero scheme は空である。
QQ 点だけの不在、有限点標本又は chosen q chart の結果ではない。

## 4. T36/T38 の completed full source に必要な閉条件

accepted completed source `MPS_SMOOTHING_COMPLETION_20261103_v62/T36.json` と
`T38.json` の**全括弧**を、p0 と四 smoothing parameter を記号のまま保持して読み直した。
全 coefficient monomial は x0 に非依存、base degree二である。
独立実装は各全 source の十 cubic から全五十微分を作り、
q=[1:0:0:0:0] の値が Laurent parameter ring 内で恒等的に零と確認した。
period 又は smoothing の特定値を代入していない。

よって period denominator の非零域全体の両 source に order≥2 零点がある。
PGL 飽和でも intrinsic 性から同じである。
全 tensor と q の incidence は上の五十式で閉、q の parameter P4 は proper なので、
「order≥2 零点を少なくとも一つ持つ」は ambient tensor 空間の閉条件である。
この議論は representative 空間でも成立し、Euler gauge 不変なので projective tensor 空間へ降りる。

accepted full-family component theorem を用いると、その閉条件は T36/T38 の全 component closure に成立する。
period denominator が零の境界へ有理 source 式を直接代入する必要はない。
full family の一般 source から閉包へ閉条件を移したのであり、
一つの既知 seed point から別の source 全体の帰属を推論したのではない。

## 5. 親と元の canonical source への正確な含意

K4a の全非零 literal curve は一つの dense actual24 orbit である。
その e=1 に上の閉条件が成立しないので、actual24 親は T36/T38 のどちらとも異なる。

さらに canonical K4a point P0 は accepted theorem により smooth unique24 点である。
もし P0 が T36 又は T38 に属せば、そこを通る full component の一意性と
両既知成分の次元二十四から、その成分は唯一の actual24 親と一致する。
すると親の dense orbit にも必要閉条件が成立するはずで、上と矛盾する。
したがって P0 は両成分の外にある。
全一軌道 K4a source は P0 に actual PGL 同値なので、同じ非所属が全 source に戻る。
source の境界全点まで非所属とは主張しない。

T37 は今回の必要条件の適用対象でなく、親候補として残す。
他の全 prior24 比較、T37 の actual 認識、新成分性・相違の全監査は別義務である。
本 review による新成分計数は零。
