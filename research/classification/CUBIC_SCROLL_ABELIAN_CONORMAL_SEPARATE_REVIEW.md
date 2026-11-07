# smooth cubic scroll・abelian conormal 全源の別作者確認

2026-11-07。producer は generic_exhaustion、reviewer は plane_collision_review。
対象は `CUBIC_SCROLL_ABELIAN_CONORMAL_NONDOMINANCE.md` の smooth scroll 部分。
**判定：全 rank零・既約被約 quintic 枝の PGL 飽和像 dimension≤23 として PASS。**
cone、plane、rank二 restriction、全 E 帰属、rank四不存在は結論しない。

## 1. 新しい線形条件は全候補に必要

全滑らか cubic scroll は S(1,2) の単一 PGL5 軌道である。
三つの二次 ideal generator と Hilbert–Burch syzygy は

```
q=(x0x3−x1x2, x0x4−x1x3, x2x4−x3²),
(x3,−x2,x0), (x4,−x3,x1).
```

この列から h0(I_S(2))=3、q(S)=0、linear-normality を得る。
既受理の Euler 全源補題により、S 上で零となる全 projective bivector は
`Π=Σ qg Ag`（Ag は任意の constant skew matrix）で尽くされ、全三十係数である。
非退化 S を点ごとに固定する projective linear field は零なので、Euler 商の隠れた
kernel によりこの source が縮むこともない。

`RANK_ZERO_SURFACE_GENERIC_FIRST_JET.md` §§1–3 を独立に精査した。
非plane S の generic multiplicity≥3 は、二点で各三重に零となる secant line を
五次 D に含ませ、Sec(S)=P4 に矛盾する。従って multiplicity は二。
rank零では四次元 cotangent Lie algebra の derived algebra は二次元 conormal ideal に入り、
normal quadratic Pfaffian の非零性からその Lie algebra は Frobenius になる。
conormal ideal が aff1 なら全 derivation が inner なので quotient の二 lift を調整して
ideal と commute させられる。残りの相互 bracket も center 零により消え、
aff1⊕C² となって rank≤2 になる。この矛盾から conormal は abelian。

したがって generic smooth locus で `{I,I}⊂I²` が必要である。
局所 projective generator qg/ell² に直す補正は、Π の係数自体が I にあるため I² に入る。
この必要条件から三 generator pair の全五 ambient first partial が S 上で零となる。
generic 消滅は integral S 上の polynomial 消滅へ延長する。
上界には逆方向の differential-square 定理を要しない。

## 2. 全159行・全kernelと全Jacobiの独立再構成

独立 checker は標準ライブラリの Fraction と疎多項式だけを用いる。
producer、共有計算 helper、SymPy、CAS を import していない。
dominant parametrization

`x=(s*u,t*u,s²*v,s*t*v,t²*v)`

へ全 first partial を引き戻し、三十列の全159行と label を再構成した。
parametrization の dominance は x2≠0 で
`(s,t,u,v)=(1,x3/x2,x0,x2)` と戻せることからも分かる。
159×30 行列の exact rank18、18 minor、free indices、完全な30×12 kernel 行列が一致した。

元の三十 basis と新しい十二 basis の双方から、Schouten polarization により
全350 projective Jacobi 行と全126実五次係数を直接作り直した。
保存された両行列と全係数が一致し、十二係数での Jacobi 二次 row rank は24であった。
これは一つの parameter point への評価ではない。

## 3. native GB を再実行しない dimension 証明

native が選んだ二十四 generator は、それぞれ元350行の非零定数倍であることを確認した。
続いて、六十二の保存 polynomial について、三百四十二の保存 lift を Fraction で展開し、
元 generator の ideal combination と完全に一致することを確認した。
grevlex leading monomial も全六十二個を独立に再計算した。

leading support の全792通りの五変数集合について、被覆しない monomial を一つずつ保存し、
六変数による被覆を実際に一つ提示した。従って minimum vertex cover は厳密に六。
これらの leading monomial は元の全 Jacobi ideal の initial ideal に属するので、
全固定 Poisson source の dimension は `12−6=6` 以下となる。
この議論は六十二個が完全な Gröbner basis であること、radical 性、native dim の正しさを使わない。
dimension が厳密に六であるという下界も主張しない。

## 4. preserving group と全 PGL 像

base SL2、relative fiber Gm、二つの fiber unipotent parameter は、scroll を保存する
忠実な六次元群を与える。abelian-conormal 条件は intrinsic なので、この群は全必要 source X を保つ。

非零 X に対する map `PGL5×X→P H0(Λ²T)` の fiber には

`(h,λ)↦(g h⁻¹, λ h·Π)`

が入る。第一成分から h が決まり、Π≠0 より λ も決まるので、この部分は七次元。
個々の Π の stabilizer が零であることを要求していない。
従って全 PGL 像は `24+6−6−1=23` 以下。
全必要 source を用いた上界なので、rank二・零又は磁気境界を途中で捨てる必要はない。

これで smooth cubic scroll 上の全 rank零・既約被約 quintic 枝は、
未知25〜27次元 G 成分の一般元から排除される。
旧三十変数 rowspace の負結果や timeout を、今回の正の証明へ読み替えていない。
新しい幾何学的必要条件と、全 ideal の exact dimension 上界が別に必要であった。

## 5. 保存物

- `code/classification_cubic_scroll_abelian_review.py`
- `research/classification/jobs/cubic_scroll_abelian_review.json`
- `reports/CLASSIFICATION_CUBIC_SCROLL_ABELIAN_REVIEW_20261107.json`

事前 finite gate 後、全再構成は0.285秒で終了。
既存 source/native report と producer bytes は変更していない。
twisted cubic cone の別局所化、他の曲面、G 全体の完成は本確認の対象外である。
