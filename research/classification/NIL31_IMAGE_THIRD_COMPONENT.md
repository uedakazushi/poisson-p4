# nil31 image内 r=−1/3：真正二法線支持から得る別の25次元成分

2026-11-07。主課題 `E_multiple_hyperplane`。
`NIL31_IMAGE_ZERO_COMPONENTS.md` §4 に残した全十変数 third 源も、
**25次元の全 Poisson 既約成分**を与える。half の25成分とも既存92成分とも異なる。
この点の全接空間は27次元であり、滑らかさ・局所成分の一意性は主張しない。

## 1. 元の全源と全周囲

全パラメータ表示、旧15源からの PGL/Euler 恒等式、全 strict Jacobi/Pfaffian、
固定点と全十微分は `CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_20261107.json`
の `branches[1]` をそのまま使う。自由変数は A0,…,A5,u,v,t,q の全十個、q≠0。
実退化式は `(q/3)h³a²` であり、h=0 が唯一の三重超平面である。

固定点は A=(1,2,3,5,7,11),u=2,v=3,t=5,q=13、tensor12倍。
全350式・全150係数の微分D、全Euler/PGL/scalar/十源微分Fについて

    rank D=98, rank F=50, DF=0, dim kerD=52.

kerD から実F50に補な真正二法線 N0,N1 を取り、`rank[F N0 N1]=52` を確認した。
Dの98本の独立な標準座標列と二法線を合わせた100次元空間Lは、F50の真の補空間である。
`rank[F L]=150` をQQ上で確認したので、Lは実際の50本の線形cutの共通零点として選べる。
固定した因子模型内の接空間や、kernelの任意の二方向を使った判定ではない。

## 2. 全Dを消す二つのfunctionalと根基原点

全三つの normal 二次 Jacobi ベクトル

    Jac(N0), B(N0,N1), Jac(N1)

に対し、全D150列を消す二本の有理functionalを構成した。
pairing を正規化すると、その二つの二次初期式は正確に

    q1 = s² − A t²,        A=113435/3825851,
    q2 = st + B t²,        B=214061/7651702.

ここで s,t は二つの normal 変数であり、全源のパラメータ t とは別である。

    B²−A = −1690119521019/58548543496804 ≠0,
    t*q1 − (s−Bt)*q2 = (B²−A)t³.

従って t³ はこの二次idealに入り、その根基は(s,t)である。
独立な98式に形式陰関数定理を適用してL内の98方向を消したとき、
二functionalはDの全像を殺すので、この二次式が残idealの真正初期式となる。
その根基が最大idealなので、真の切片自体の被約支持は原点、Krull次元は零である。

50本の実際の線形cutに対するKrull次元不等式は元の150係数空間の局所次元を50以下に抑える。
全実F50からの下界と一致するので、全十源の既約 PGL 飽和像閉包は全空間成分となる。
Euler24方向と射影scalar一方向を引いた次元は25。
全Jacobiスキームの接次元は52−25=27のままであり、精密なnil指数、一般被約性、
形式的source×normal積又はこの点の局所成分一意性は含意しない。
q=0 は元の既約多項式源の閉境界として同じ成分へ含まれる。

## 3. half 成分および既存92成分との相違

half と third はともに一般Dが二つの異なる超平面の3,2型である。
因子次数だけでは両者を区別できない。
唯一の三重超平面を無限遠に取った canonical affine modular 場を使う。

- half の三重超平面は old_h=0。旧全源の modular 場は常に
  `old_y∂old_x+old_z∂old_y` であり、線形部の Jordan 型は(3,1)、定数成分はimage内。
  affine零点を持つ。
- third の三重超平面は h=0（old_z=0）。全十変数で直接 divergence を取ると
  `−(2/3)(c∂b+∂c)`。線形部は rank1、A²=0、定数成分は長さ2鎖のtopに残り、
  affine零点を持たない。

PGL同値は内在的な唯一の三重超平面を保ち、その無限遠chartではaffine同値となる。
modular場は非零scalarとaffine共役を除いて上の型・零点の有無を保つ。
これは両源の一標本だけでなく、非零p,qを保つ**全稠密既約源**で成立する。
もし同じ25次元成分なら、両PGL飽和constructible像はその非空開上で交わるため矛盾する。
従って両25成分は相異なる。

既存92成分との比較はhalfと同様である。次元25以外の成分とは次元から異なる。
全third閉包のrank4部分はproperな掛け算写像の閉条件 `D=H³L²` に入るため、
SNC点を持つT01/T39/T40、滑らかなrank5残二次式を持つpositiveNEW25、
三乗因子を持たない一般独立221型の既存nil31成分と異なる。
一般rank2の成分とも閉rank条件で異なる。これは成分性を得た後の比較である。

## 4. 証明書

- 全源：`CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_20261107.json`
- 元の二次計算：`CLASSIFICATION_NIL31_IMAGE_THIRD_QUADRATIC_SUPPORT_20261107.json`
- 正規化functional・根基恒等式・真のcutの再構成：
  `CLASSIFICATION_NIL31_IMAGE_THIRD_COMPONENT_20261107.json`
- standalone再構成：`code/classification_nil31_image_third_support.py`
- 有限契約：`jobs/nil31_image_third_quadratic_support.json`

全Dを消すfunctional、全三pairing、真の補空間、非零minor、全源divergenceをQQ上で再構成した。
二次計算2.009秒、正規化と再構成1.797秒。独立作者による全周囲・全源reviewは別途行い、
統合担当はそのPASS後にだけhalfとthirdの二成分を計数する。
