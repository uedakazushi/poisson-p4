# cubic scroll の blowup：全 P2-bundle と残る quartic 境界

2026-11-07。主義務 G の次の幾何入口。新しい dimension 上界・空性・成分は主張しない。
以下の incidence と line-bundle 計算は手で証明する。
Poisson 射又は normalization の flat 性を仮定していない。

## 1. blowup は全 base 上の P2-bundle

smooth cubic scroll の三 quadrics を
`q=(x0x3−x1x2,x0x4−x1x3,x2x4−x3²)` とする。
`[a:b:c]` を新しい `P2` の座標とし、incidence

`Y={a x3−b x2+c x0=0, a x4−b x3+c x1=0}⊂P4×P2`

を取る。二行の coefficient matrix は全 `[a:b:c]` で rank 二である。
`c≠0` なら `x0,x1` の pivots、`c=0` でも `(a,b)≠(0,0)` により最後三列の二行は独立。
従って

`0→F→O_P2^5→O_P2(1)^2→0`

の `F` は rank 三 vector bundle、`Y=P(F)`（lines convention）。
これは全 base・全 boundary を含む smooth P2-bundle である。
`c1(F)=−2B`, `c2(F)=3B²`。

Hilbert–Burch の二つの syzygy から、`P4\S` の quadric-net map
`x↦[q0:q1:q2]` の graph は Y の open dense subset。
Y は既約なので、その graph closure は Y 全体に等しい。
従って `Y=Bl_S P4`。Symmetric algebra と Rees algebra の未証明の flat 性を使わず、
全 incidence の既約性と dense graph で同定した。

`H` を P4 hyperplane の pullback、`E` を exceptional divisor とすると

`B=2H−E`, `E=2H−B`, `K_Y=−3H−B`。

実 quintic D の S に沿う一般 multiplicity 二により、strict transform は

`D'=5H−2E=H+2B`。

よって D' は section `s_D:F→O_P2(2)` の零となる **relative hyperplane**。
この section は basepoints を持ってよく、その点では fiber が丸ごと P2 になる。
その補集合を捨てて P1-bundle と呼んではいけない。

## 2. exceptional relative conic と全 discriminant

`c≠0` では `x0=(b x2−a x3)/c`, `x1=(b x3−a x4)/c`。
この plane の中で `E` は `x2x4−x3²=0` という smooth conic。
`c=0` では

`(x2,x3,x4)=t(a²,ab,b²)`

となり、plane coordinates は `(x0,x1,t)`。
全 scroll ideal の制限は `t(bx0−ax1)` の scalar multiples だから、E の fiber は

`{t=0} ∪ {bx0−ax1=0}`。

二つは常に異なる lines。第一は negative directrix `C0` の実像、第二は scroll の ruling line。
従ってこの conic bundle の discriminant は正確に line `L={c=0}`、double-line fiber はない。
`E|_fiber` は degree 二、`D'|_fiber` は degree 一である。

## 3. 残 A=(2,5) の V=C0 と、必須の四つの full fibers

`CUBIC_SCROLL_CANONICAL_ROOT_EXCLUSION.md` 後の二つの未閉型は、どちらも
normal quadratic の vertical divisor `V=C0` を持つ。
従って `D'∩E` は `E_(C0)` を一つの divisor として含む。

二つの projection を区別する。`E_(C0)→C0` は original exceptional P1-bundle、
一方 quadric-net map `p:Y→P2` では `E_(C0)≅C0×L→L` である。
実際 C0 上の `dq0,dq1` は normal fibers を独立に張り、`dq2=0`。
従って generic `y∈L` で、relative line `D'_y` は conic fiber の `C0` component に等しい。

`L≅P1` 上では前節の explicit coordinates から

`F|L=O_L² ⊕ O_L(−2)`。

第一 summand は `x0,x1` という fixed vertex plane-vector-space、最後の embedding は
`t↦t(a²,ab,b²)`。
`s_D|L` は `O_L²` を零にするため、その全情報は

`kappa:O_L(−2)→O_L(2)`,

すなわち **全 binary quartic `kappa(a,b)`** である。
`kappa` が恒等零なら `p^*L⊂D'` となり、D' が vertical divisor を成分に持つ。
実 D の既約性と D' の class `H+2B` に反するので `kappa≠0`。

従って `Z={kappa=0}⊂L` は長さ四の zero scheme。
その全点では `s_D` は F の全 fiber 上で零、従って `D'_y=P2` となる。
squarefree の四点だけではなく、211/22/31/4 の全 multiplicity partitions を含む。
base の L 以外に追加 isolated zero が無いという主張もしていない。

このため残二型を「smooth P1-bundle の hyperplane family」として数える方法は誤りである。
必須の長さ四の full-fiber boundary、そこでの D' の特異点及び conductor を扱う必要がある。

## 4. 未証明の次の矢印

quadric-net map が全 Poisson source に対して Poisson 射となることは本 note では証明していない。
raw homogeneous bracket の simultaneous Euler gauge、全 symbolic-square/ordinary-square の degree 四同一視が
必要であり、単なる q_i の比例から推測しない。
また D' の normalization がこの base 上 flat、全 full fibers が regular 又は小解消を持つとの仮定はしない。

次の候補は、全 quartic partitions における closed inverse-residue の conductor 制約、
又は元の全 source に対する whole zero/nonzero proof である。
これらを完了するまでは smooth scroll rank 二の A25 二型と cubic cone は未解決であり、
旧 generic-minor/GB/jet の反復をこの観察で正当化しない。
