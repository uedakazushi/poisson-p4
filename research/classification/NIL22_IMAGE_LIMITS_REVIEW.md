# nil22 image の全九変数二源・全八変数一源：regular limit の独立査読

2026-11-07。`E_multiple_hyperplane`。
`quadratic_B` の全源が既済 half25 成分に、`triple_C` の全源が既済 third25 成分に入ることを、
別作者が厳密な全係数恒等式として確認した。**両全源の被覆を受理する。新成分数は0。**
残る nil22 image 源の帰属はこの結論に含めない。
追記した第5節では `double_D` の全八変数源も既済 CN31_25 へ被覆する。

独立 checker は既存の `code/classification_nil22_image_components_review.py` に加えた
`--limits` mode。新証明書は `reports/CLASSIFICATION_NIL22_IMAGE_LIMITS_REVIEW_20261107.json`、
有限 job は `research/classification/jobs/nil22_image_limits_review.json`。
生産者・helper・CAS・NumPy の import はなく、汎用 Fraction/sparse 算術で再構成した。
実測0.0315秒、最大RSS14848 KiB、全検査 PASS。
直前の二新23成分の数学関数と元証明書 bytes は変更せず、旧 job の共有 runner hash のみを
専用 CLI 拡張理由付きで更新した。

## 1. 対象源と依存

元入力 `CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json` の
case8=`quadratic_B`、case6=`triple_C` を用いる。
両者の全自由係数は `m0,…,m5,n,c,d` の九個で、いずれも反転しない。
全六磁気 F0 と canonical S,V から元の全源を作り、保存された全150係数と一致させた。

適用する既済成分定理は、独立証明済みの
`CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_REVIEW_20261107.json` と
`CLASSIFICATION_NIL31_IMAGE_THIRD_COMPONENT_REVIEW_20261107.json`。
half の元15係数源からの全係数移送と、third の正規化全十源が対象である。
成分の局所一意性や endpoint での接空間階数は使わない。

## 2. quadratic_B：元の無限遠を保つ曲線

元 half15源で `n6=n7=n9=n11=0,n10=−1/2` とし、残る十係数を
生産者報告にある Laurent 式へ代入する。
それらの式を全て独立に再入力して保存値と比較した。
特に `n5=−1/ε` は ε非零で必ず非零となる。

元の四変数を `(X,Y,Z,W)` として実座標変換は

    (q1,q2,x,y)=(Y−W/ε, εX, Z, W), det=−ε.

この可逆 GL4 変換は h を保つ PGL5 変換である。
全源を代入し、全括弧を tensor として押し出した後、全150係数の ε の負冪が
完全に相殺することを検査した。結果は九係数と ε の多項式で、各空間次数は2以下。
ε=0 の全150係数が元の quadratic_B 全九源に正確に一致した。

全 ε と全九係数で strict Jacobi10式と projective Jacobi350係数は恒等零。
斉次 Pf は `−h²x²/2`、affine modular は `(x,y+εq1,0,0)`。
従って ε非零では全て既済 half25 全十源の実 PGL 像であり、rank4である。

## 3. triple_C：内在的な三重超平面を無限遠へ移す曲線

元の `(h,q1,q2,x,y)` に対して permutation `(3,1,2,0,4)` を施し、
元の x を新 h にする。変換後の tensor Q から

    E∧R,  R=(0,h/3,y/3,0,0)

を引いて h-Casimir representative とする。
R は線形であり、h の逆元は最終式にも途中の正当化にも必要ない。
全150係数でこの実 PGL/Euler 恒等式を確認した。

既済 third の正規化十源を `(h,a,b,C,D)` で表示し、
その全係数を独立の明示式と照合してから Laurent 値を代入する。
特に source period `q=−3/ε` は ε非零で必ず非零。
実 GL4 変換は

    (q1,q2,x,y)=(C−D/ε, εb, a, D), det=−ε.

全150係数の負冪相殺を検査した結果、九係数と ε の多項式 tensor が得られた。
ε=0 は上の変更後の triple_C 全源と完全一致する。
全 strict10/projective350 Jacobi は恒等零、斉次 Pf は `h²x²`。
affine modular は `−(2/3)(1,y+εq1,0,0)` である。
従って ε非零では既済 third25 全十源の実 PGL 像に入る。

## 4. 全支持への閉性の適用

任意の九係数を固定しても、上の式はそのまま成立する。
反転したものは punctured curve の ε のみで、磁気係数や n,c,d の例外はない。
最終 tensor の全係数は ε=0 に正則であり、Pf も恒等的に非零なので射影点としても延長する。
ε非零の全点が、それぞれ同じ固定された既済の閉成分 half25／third25 に入る。
従って閉性により ε=0 もその成分に入る。
triple_C について最後に可逆な PGL と Euler representative を戻す。

これで両元の全九変数源と全ての係数零境界を被覆した。
有限 jet の一致、Pf 因子型の一致、特殊 endpoint の共有から帰属を推論する議論ではない。
実 punctured PGL/source 恒等式と、全係数で検証した正則な endpoint による結論である。

## 5. double_D の全八変数源も CN31_25 へ入る

root の追加依頼を受け、同じ独立 checker の `--double-d` mode で検証した。
`CLASSIFICATION_NIL22_IMAGE_DOUBLE_D_REVIEW_20261107.json`、
`jobs/nil22_image_double_d_review.json` が新たな算術記録・有限契約。
0.0184秒、最大RSS14848 KiBで PASS。
先行二源の immutable 報告は保全し、旧 job の runner hash は追加 mode の理由付きで更新した。
generic parameter index へ変更した先行 `--limits` も有限契約内で再生し PASS した。

case3=`double_D` の自由係数は六磁気 m0,…,m5 と a,u の全八個。
canonical `D3=x²y,S=diag(u,0),V0=(0,−1−(a+1/2)u),V1=diag(a,−1−2a)`
から元の全 M と F0 を独立に構成し、全150係数で元 source と一致させた。
元 nil31-image15源の全 Laurent 置換は生産者報告の全15項と一致する。
さらに `MODULAR_NIL31_COMPONENT.md` §1 の全六 Jacobi compatibility を残し、
これらが全て恒等零となることを独立に検査した。
`n7=−2/ε` は puncture 全域で非零であり、既済の **全 n7 開** の定理を適用できる。
その適用先は `CLASSIFICATION_RANK2_NIL_31_COMPONENT_REVIEW_20261107.json` で
独立証明済みの固定成分 CN31_25 である。

実 GL 変換は第2節と同じ determinant −ε の行列。
全150係数で Laurent 負冪が相殺し、全八係数と ε の正則多項式 tensor となる。
ε=0 は元の double_D 全源に正確に一致する。
全 strict10/projective350 Jacobi は恒等零、斉次 Pf は

    −h x²(y+u h/2),

affine modular は `(x,y+εq1,0,0)`。
Pf の x²y 係数が非零なので、u=0 や a=0,−1/2 等でも rank4。
磁気・a・u の例外的な非零条件は一切ない。
従って任意の固定八係数について puncture は同じ CN31_25 内に入り、
閉性から endpoint も入る。全八変数 source と全係数零境界を被覆する。
これも成分数を追加せず、残る六源の帰属には結論を与えない。
