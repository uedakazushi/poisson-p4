# Smooth cubic scroll rank零：修正全支持と既約五次不存在の独立確認

2026-11-08、PRIMARY `E_reduced_nonnormal`。
**判定 PASS：smooth cubic scroll 上で零となる rank四・被約既約 actual quintic 枝は全空。**
producer の旧四源 cover は棄却したまま保存し、正しい L6∪Segre五の全被約支持で置き換える。
可約 divisor を持つ全点の既知成分帰属は本結論に含めない。新成分数は零。

## 1. 入り口と全 source の適用範囲

既受理の `CUBIC_SCROLL_ABELIAN_CONORMAL_SEPARATE_REVIEW.md` を利用する。
全 smooth cubic scroll は一つの PGL5 軌道であり、標準 ideal は

```
q0=x0*x3-x1*x2, q1=x0*x4-x1*x3, q2=x2*x4-x3².
```

Euler 全源補題から、その上で零となる全 quadratic projective bivector は
`Σ qg*Ag`、Ag constant skew、の全三十係数で尽くされる。
被約既約 quintic D の rank零非plane面では、secant-line degree 論により generic multiplicity は二。
conormal ideal が nonabelian aff1 なら Frobenius linearization は
aff1⊕C² へ分裂し rank二となる矛盾から、conormal は abelian。
従って accepted complete30×12 source は全対象に必要である。
この入り口の仮定を、新しい四源算術や dimension 上界だけで代用していない。

旧独立 receipt は全30源、全159線形行、全30×12 kernel、
全350 Jacobi と全126 signed Pfaffian source rows、二十四 generator の元行との対応を確認済み。
今回、新rank・minor・operator・pointは計算せず、その complete source を使う。

## 2. 修正全被約支持を両方向に確認

元順序 r0,…,r11 を

```
(p,a,s,q,v,f,r,h,u,w,b,t), delta=b-a, eta=t-2s
```

とする。全 unit/zero split は次を尽くす。

- eta≠0：元 g19,g1,g9,g20,g24,g3 から L6。
- eta=0,v≠0：g4/g8,g5,g3,g13,g12 から L6。
- eta=v=0,a-2f≠0：g2-g3,g6+6g3,g22,g7,g24,g12 から L6。
- 残る eta=v=0,a=2f：t=2s,b=2f+delta の全 subspace。

L6 は `t=s,b=a,w=-v,f=u=h=0`、自由六係数である。
前三枝の各 forcing polynomial identity を独立 Fraction で確認した。

最後の subspace では `P=p-2u,L=q+h` と置くと、元全二十四式は
`CUBIC_SCROLL_ZERO_TWO_WHOLE_SOURCES_V2.md` §2 の正しい八式 A1,…,A8 と
**定数係数で両方向に同値**。全復元と全逆を直接展開して確認した。
旧誤推論 fu=0 は使わず、fu非零の全補集合を保持する。

```
M=[[s,f,u],[delta,w,h],[P,-L,-r]].
```

七式は M の七 minor、残りは二 minor x=uL-fr,y=-Ph-delta*r の y+2x=0。
全九 minor の符号と八式の対応を別に直接確認した。
M の第一 row が非零なら、上二 row の全 minor 零から row2=k*row1。
既知 rows1,3／cols1,3 minor が零なので y=k*0=0、続いて x=0。
第一 row が零なら x=0、続いて y=0。
従って全 reduced support は rank M≤1、その逆も全 minor 零から従う。
scheme ideal の radical 性や missing minors の literal ideal membership は要求しない。

全 rank≤1 matrix は C 上外積 `(A,B,C)^t*(D,E,F)` で表され、零 matrix も含む。
元 r への全 polynomial map は

```
(CD+2AF,2AE,AD,-CE-BF,0,AE,-CF,BF,AF,BE,2AE+BD,2AD).
```

これと L6 の全六 parameter map が、全十二 source の reduced support を尽くす。
chosen nonzero chart、dimension 一致、又は部分族の generic factor を外挿した cover ではない。

## 3. 元 tensor・全350・全126と実因子の別算術

新 checker は標準ライブラリ Fraction と旧独立 checker の多項式基本演算だけを用いる。
producer、SymPy、CASを import しない。二つの全 map は独立に記述した。
元30×12 matrix と q0,q1,q2 から全十括弧を組み、
表示された全150係数、全350 projective Jacobi、保存 source の全350行、
定義からの signed Euler-Pfaffian と保存全126係数を別に再構成した。
両全源で全て一致する。

L6 の actual D は全 parameter 上で正確に positive physical degree 一×四。
全 Segre五は、外積 parameter の全六変数を保持したまま

```
D=(1/2)*q2*(D*q0-E*q1-F*q2)*H1
```

であり、H1 は producer V2 report に保存された全 linear form。
この全 polynomial 積を直接展開し、元 actual D と全係数一致することを確認した。
各因子の physical homogeneous degree 二・二・一も個別に確認した。
parameter-only scalar の処理も含め、因子の存在を old F/U/S 部分族から外挿していない。

全 source parameter の任意の specialization 後、各 physical factor は零又は元の正次数を保つ。
いずれかが零なら D=0、それ以外なら少なくとも二つの positive-degree factors を持つ。
重複因子の場合も被約既約 quintic にはならない。
outer-product rank零、全 parameter零、全反転補集合を含め、全対象を尽くしている。
PGL5 inverse はこの因子性を保つので、標準 scroll だけの座標依存結論ではない。

## 4. 保存物と検証範囲

- checker：`code/classification_cubic_scroll_segre_support_review.py`
- job：`jobs/cubic_scroll_segre_support_separate_review.json`、実行前 gate PASS。
- receipt：`reports/CLASSIFICATION_CUBIC_SCROLL_SEGRE_SUPPORT_SEPARATE_REVIEW_20261108.json`
- 全1095確認、0.0486秒、RSS15832KiB。

元四源 cover の失敗は `CUBIC_SCROLL_FOUR_SOURCE_COVER_REJECTION.md` と exact first-failure archive に残す。
四部分源の valid factor identities と旧 whole-cover boolean を混同しない。
修正 V2 は別の全六 parameter源と正しい全八式で受理している。

twisted cubic cone、plane、scroll 上の rank二制限、他の導手面、
及び可約 D の全点既知成分帰属はこの確認の対象外。
