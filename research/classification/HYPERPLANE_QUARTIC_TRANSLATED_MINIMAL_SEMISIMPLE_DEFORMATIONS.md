# 最小冪零・長 block translation：残る二軌道の全 MPS 被覆

2026-11-07。課題 `E_reduced_nonnormal`、純冪零型10。
`HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_NORMAL_FORMS.md` の全90係数からの還元と
別作者の全三軌道確認を入力とする。
同稿の C 軌道は既に実曲線から半単純型20へ被覆した。
本稿は残る A,B の**全軌道**を、元の二次 tensor の実多項式 Poisson 曲線で被覆する。
従って型10、単純 H と既約四次 F を持つ全支持は既存 MPS40 の閉和へ入る。
成分下界を追加せず、他の冪零型・全分類の完了も主張しない。

## 1. 全係数を保持する曲線

五変数の順序を(h,x,y,z,w)とし、P の h 行を零とする。
以下の六括弧は全て斉次二次であり、全 t に対して strict Jacobi を満たす。

### A：線形 potential の代表

    Pxy=hx/2+t²x²/4−y²/4,
    Pxz=−yz, Pxw=−h²+yw/2,
    Pyz=−hz−t²xz,
    Pyw=−h²t+hw/2+t²xw/2,
    Pzw=hx+t²x²/2−twz−y²/2.                         (A_t)

divergence は

    Z_t=y∂x+(h+t²x)∂y+t z∂z−t w∂w.

t=0 は三正常形の A と、定数項も含め全六括弧で正確に一致する。

### B：三次 potential の代表

    Pxy=hx/2+t²x²/4−y²/4,
    Pxz=−yz/3, Pxw=−yw/6−3z²,
    Pyz=−hz/3−t²xz/3,
    Pyw=−hw/6−t²xw/6−3tz²,
    Pzw=hx+t²x²/2−twz/3−y²/2.                     (B_t)

divergence は

    Z_t=y∂x+(h+t²x)∂y+(t/3)z∂z−(t/3)w∂w.

t=0 は正常形 B と全六括弧で正確に一致する。
初期 Pfaffian は、それぞれ

    (xh−y²/2)²/2+h³z,
    (xh−y²/2)²/2+hz³.

これは既約四次であり、h 因子は持たない。
両曲線の全 Pfaffian における y⁴ 係数は常に1/8なので、四次式の射影点は全 t で定義される。

## 2. 全非零 t の半単純 cofactor を認識する

二曲線を一様に κ=1又は1/3 と書く。
h=1、t≠0 で modular zero は

    (x,y,z,w)=(−1/t²,0,0,0).

この点の定数交代形式の Pfaffian は両方とも `1/(8t⁴)` で、全非零 t で非退化。
ξ=x+1/t² とすると線形 cofactor は

    A_t=[[0,1,0,0],[t²,0,0,0],[0,0,κt,0],[0,0,0,−κt]].

実固有座標

    (a,b,c,d)=(tξ−y,tξ+y,z,w)

の行列式は2tで、固有値は(−t,t,κt,−κt)である。
どの t≠0 にも分母零の例外はない。
κ=1では tensor を t で割り、固有空間を並べ替えると既済型20の cofactor となる。
定数交代形式は非退化で cofactor に保たれるため、負・正固有空間間の crossing 行列は rank 二。
κ=1/3では scalar 3/t と座標の並べ替えで型7、固有値(−3,−1,1,3)となる。

従って実曲線上の非零 t の点は、完全に処理済みの半単純 cofactor 型へ入る。
同じ固有値の場を持つという観察だけでなく、全 tensor の Jacobi・次数・定数部を確認した適用である。

## 3. 既約四次条件と閉和への全点包含

固定次数の可約多項式の集合は、各次数分割の二つの射影多項式空間からの積写像の
有限和の像であり閉じる。
上の F_0 は既約なので、F_t が既約となる非空 Zariski 開集合が t=0 の近傍にある。
全 t に対する実 Poisson 曲線であるため、平坦性や有限 jet の持上げは仮定しない。

さらに A では座標重み(0,2,1,4,−1)、character1、
B では(0,6,3,4,5)、character3に対して、曲線の全 monomial は同変である。
point を固定する compensated PGL・scalar トーラスは t をそれぞれ t/a、t/a³ へ送る。
複素数体では全非零 t が同じ実軌道になる。
従って非空開集合で確認した rank 四・既約四次の条件は全 t≠0 に成立する。
定数部の上の非零 Pfaffian は rank 四について直接の全域確認にもなっている。

`HYPERPLANE_QUARTIC_SEMISIMPLE_SUPPORT.md` の型7全支持定理と
`HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md` の型20全支持定理を適用すると、
t≠0 の両曲線全体が既存 MPS40 の有限閉和に入る。
その閉性から t=0 も入る。
連結な PGL5 と非零 scalar が各成分を保つので、A,B の全実軌道も同じ閉和に入る。

これは、MPS 点が一つ含まれるだけの任意曲線から他点の所属を推測する論理ではない。
非零 t の全像の包含を全支持定理から得た上で、その真の閉包を取っている。
既済 C の全点被覆と合わせて、全型10の支持が閉じた。

## 4. 有限な探索と数学的証明の区別

実装 `code/classification_translated_minimal_quartic.py` の
`--cofactor-cubic` と `--cofactor-linear` は、元の全斉次二次 monomial のうち
上のトーラス重みに整合する係数を全て置き、divergence と cofactor 不変性を全線形で解いた。
これは「ここまでの jet が解けた」という計算ではない。
得られた (A_t),(B_t) を**全 t の多項式恒等式**として全 strict Jacobi に再代入した。

B の許容 monomial は15係数、t次数は先験的に2以下。
A は43係数、全許容 monomial における最大 t次数は9。
A の実曲線自体は次数2のみである。
初期 job の prose にあった bound8 は、yw の z² 項の相対重み−9を落とした手数えだった。
producer は cutoff を一切使わず、全 monomial を先に列挙して43係数・max9を保存しており、
元の計算にも証明書にも欠落はない。job の説明だけを9へ訂正した。

得られた正例の証明書は

- `CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_CUBIC_COFACTOR_CURVE_20261107.json`
- `CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_LINEAR_COFACTOR_CURVE_20261107.json`

有限契約は同じ prefix の `jobs/..._cubic_cofactor_curve.json` と
`jobs/..._linear_cofactor_curve.json`。実測はそれぞれ0.326秒、0.844秒。
探索途中の不整合を他の変形幾何の非存在や新成分性へ用いていない。

## 5. 別算術による全曲線の確認

`code/classification_translated_minimal_quartic_review.py --cofactor-curves` は
producer・producer engine・CAS を import せず、標準ライブラリ Fraction と
別作者の nil22 independent reviewer の疎多項式演算だけを用いる。
両曲線について独立に次を再構成した。

- 元の全代表との t=0 完全一致、全斉次二次係数。
- 全十 strict Jacobi と全350 projective Jacobi、全 Pfaffian・divergence。
- 全座標 monomial のトーラス同変性。
- 全 t≠0 の modular zero、定数 Pfaffian、実固有座標の行列式2tと全固有値。
- `A_t C+C A_t^T=0` と非退化定数形式から、型20への crossing rank 二の適用。

全て PASS、実測0.00646秒、最大RSS 14,336 KiB。
独立証明書は
`CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_COFACTOR_CURVES_REVIEW_20261107.json`、
契約は `jobs/hyperplane_quartic_translated_minimal_cofactor_curves_review.json`。
外部の MPS theorem と既済半単純全支持証明自体の再実行とは区別する。
この結果の成分数追加は零であり、全型10の既知閉和への帰属だけを主張する。
