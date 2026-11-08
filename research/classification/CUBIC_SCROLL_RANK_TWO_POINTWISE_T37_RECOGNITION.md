# Smooth cubic scroll・rank二の全actual既約quintic源はT37へ

2026-11-08、PRIMARY E、producer全source接続。最終独立receiptは別に固定する。
新operator／anchor／M1／profile／jet／CASを本接続のために実行していない。
新componentの増分は零。

**定理候補。** P4上のrank四Poisson構造πの実退化因子Dが被約既約quinticであり、
非正規支持にsmooth cubic scroll Sがあり、π|Sが一般rank二なら、πは既知閉成分T37に入る。
R=div(π|S)のreduced／irreducible、clean、normalization-flat又はgeneric四full-P2 fibersを仮定しない。

## 1. 全normal-rootの入口を戻す

既受理 `CUBIC_SCROLL_NORMAL_ROOT_DIRECTRIX_REDUCTION.md` と
`CUBIC_SCROLL_CANONICAL_ROOT_EXCLUSION.md` は、この全sourceをA25二型へ戻す。
rank一normal quadraticはprimitive root A=2C0+5fとdiv gamma=C0、
rank二は同じclassの二primitive rowsとdiv gamma=C0。
canonical dq2-root、他primitive root classes、全zero／vertical complementを既排除した後の
**全源還元**であり、oldG像dimensionだけの判定ではない。

`CUBIC_SCROLL_A25_WHOLE_RESIDUE_CONSTANT_REDUCTION.md` は、このwhole入口から
同じflat fiber constant kをLとΩの両divisorで比較した。
node22はk=2とk=−2の矛盾、node31は既排除。
cuspは全Ωmultiple／lambda零を戻してmu=3lambdaを得る。
全actual Aut(S)=Aut(P2,p)共役と許されたpresentation変換の明示逆により
全cusp pairは(N3,diag(1,3,0))へ戻る。
global exceptional-restriction injectivityが元projective πも一意にするので、
**whole A25 actual sourceは一actual PGL orbit又は空**である。

従って残る義務は一つのgenuine cusp-scroll pointのT37帰属であり、
以下は別担当の凍結済みpositive exact sourceを使ってその存在と帰属を与える。

## 2. 既存のgenuine seedは仮定を全て満たす

`PLANE_M2_CONCURRENT_CUSP_SCROLL_AND_MUTUAL_MAP.md` のcanonical A、period5/2を使う。
座標(a,b,c,u,v)、全十括弧は既実行
`CLASSIFICATION_PLANE_M2_CONCURRENT_DOUBLE_WHOLE_NORMALIZER_20261108.json` cases[0]に保存されている。
これは元fullA8 sourceからの実normalizerによるcanonicalで、恣意的sampleではない。
本接続にはその一つのexplicit actual tensorだけで十分であり、他sourceの尽くしを仮定しない。

同じsourceのu=1 chartでactual五次式は

```
F_A=−256/243*[c+9(a−4b)v/16−v³/8]²
    −3/2*[a−2v²/9]³.
```

括弧内の二つは(a,c)を置き換えるregular coordinates。
従ってその二つの零schemeのdense開でnormal transverse typeはordinary A2 cusp、normal quadratic rank一。
全closure S'は2×3 matrix

`[[u,2v,9b],[v,9a,4c]]`

の全minorsで定まるsmooth S(1,2)。matrixは可逆なrow／column尺度とlinear variable順序で
standard scrollへ戻り、独立な五ambient variablesを使う。
全reduced SingFはこのS'と元plane(u=v0)の合併だが、追加planeを消して仮定していない。

intrinsic affine bracketのexact readbackは

`{b,v}|S'=4(9b−v²)/27`。

従ってS'上のrestriction rankは一般二。
actualFはcのprimitive quadraticで、discriminantのv零restrictionは
−512*a³*u⁵/81という非零非平方なので絶対既約。
actual rank四、smooth scroll、generic cusp、restriction rank二の全仮定を満たす。
既canonical-root排除により、このseedもA25 cusp sourceである。

## 3. 名前付きT37への既存actual positive bridge

別担当の凍結 `PLANE_M2_CONCURRENT_DOUBLE_WHOLE_T37_COVERAGE.md` と
`CLASSIFICATION_PLANE_M2_CONCURRENT_ACTUAL_T37_FRAME_20261108.json` を使う。
canonical Aとperiod−1/2のcanonical Bにはactual RAW恒等式

`(a,b,c,u,v)→(−a,b−a,c/3,u,v), M_*Pi_A=Pi_B/3`

がある。正逆と全150係数を確認している。

同じBから作ったexact polynomial actual family P(e)は、e零でB、e非零ではCstarにより
e一へ移る。全350Jacobiは全eで恒等零。
e一のtensorはQQ[theta]/(theta²+118)上の明示actual GL5と非零scalar、Euler補正により
**全correction込みcompleted T37 source**のtensorと全150係数で一致する。
全inverse及び全126Euler-Pfaffianも一致する。

このmembershipの接続には、fullambient smoothness又は親の唯一性を使う必要はない。
e一が閉PGL-stable T37に入り、全e非零がそのactual orbitに入り、
polynomial族の閉性からe零のB、したがってAもT37に入る。
native／数値近似／graph分母chartの陰性を陽性帰属の代わりにしていない。

## 4. Whole smooth-scroll sourceへ戻す

§2のseed Aはwhole A25 cuspの一actual pointであり、§3でT37に入る。
§1の全actual normalizer／presentation逆とglobal extension一意性から、
任意のA25 cusp pointはAのactual PGL orbitに入る。
したがってwhole cusp sourceはT37へ入る。node全空とcanonical-root等の全排除を戻し、
冒頭のsmooth cubic scroll rank二の**全pointwise**結論を得る。

これはseedの既知帰属を未読higherへ延長した推測ではない。
全sourceの実同一軌道への橋と、independentに固定した一pointの陽性帰属を合成した証明である。
rank零smooth scroll、twisted cubic cone、他の導手次数、他factor型又は全Poisson分類を
この定理だけから結論しない。
