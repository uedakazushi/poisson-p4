# Whole T37 の rank四境界：次数五以下の Poisson 曲面が退化因子を覆う

2026-11-08、PRIMARY E、純手証明。新しい point、operator、minor、normal/profile、jet、GB 又は CAS は実行しない。証明は actual canonical T37 の多項式・既独立検証済み residue 恒等式から直接出発する。canonical index二の必要源を normalization の十分条件として逆用しない。

## 1. Actual canonical pair と次数五の base surface

凍結済み `CANONICAL_INDEX_TWO_NODE_WHOLE_T37_COVERAGE.md` §§2–3 と root 別作者検証を入力とする。座標 (a,b,c,d,e) で

```
q=ac−b², A=c²+bd, B=a²+be,
C=c²e+a²d+2acb+bde−b³,
F=det[[C,A,B],[A,a,b],[B,b,c]],
g1=bB−cA, g2=bA−aB.
```

F は絶対既約 quintic、actual divisor は非零定数倍の F。掲載された全150係数の actual Poisson tensor は、真 completed T37 の全 period 一の tensor と actual PGL/Euler で一致する。同じ true all-one 点を K4a の actual24 親へ移す `K4A_T37_ACTUAL_MAP_SEPARATE_REVIEW.md` と、`K4A_COMPLETE_TORUS_SLICE_SEPARATE_REVIEW.md` の非零parameter全一軌道／C=T37 の定理を合わせると、この canonical tensor 自身の PGL24軌道は whole T37 に稠密である。normalization の converse ではなく、この actual 同値と既成分次元の一致を使う。D の smooth rank二 locus における inverse residue は、非零定数を除いて dlog(g1/g2)。この最後の主張は全六二形式の元 polynomial identity として既に独立検証されている。

全五座標で literal identity

```
bF=g1*g2−q³                                               (1)
```

が成り立つ。F=Cq−cA²+2bAB−aB² と bC=AB−q² を代入すれば直接展開できる。

S の homogeneous ideal を I=(q,g1,g2) とする。これは行次数 (2,1,1) の行列

```
[[A,B],[a,b],[b,c]]
```

の二次 max minors である。q は rank三の既約 quadric で、g1 は q の倍数ではない（q=0,(a,b,c)=(0,0,1) では g1=−1）。従って grade は少なくとも二、determinantal bound により正確に二である。Hilbert–Burch の全 resolution は

```
0 → O(−4)^2 → O(−2)⊕O(−3)^2 → I → 0.
```

よって S は pure codimension二・Cohen–Macaulay・saturated で、Hilbert numerator 1−t²−2t³+2t⁴ の二階微分から degree は五。

b=0 では (q,g1,g2)=(ac,−c³,−a³) なので reduced support は a=c=0 の line。従って surface component は b=0 に含まれない。b unit 上では ac=b² により a,c も unit、g1 の e 係数 b² が unit であり、

```
a*g1+b*g2=−A*q
```

より I=(q,g1)。q で c、g1 で e が一意に消え、b=1 chart は a∈C*,d 自由の smooth integral surface となる。CM の unmixedness と b unit open の稠密性から全 S は integral である。このため以下の [S] は generic multiplicity 一の degree五 cycle。

## 2. 全一般点を通る次数五以下の invariant surface

κ を一般の非零定数とし Gκ=g1−κg2。Gκ は非零 cubic で、既約 F の倍数ではない。従って CI(F,Gκ) は pure surface cycle、degree十五。

S の generic point では b unit、(1) から F∈I²、Gκ∈I。ambient の generic codimension二 regular local ring で二式の intersection multiplicity は少なくとも二。従って

```
[CI(F,Gκ)] ≥ 2[S],     degree([CI]−2[S]) ≤ 5.             (2)
```

特に S 外の全 residual surface component の degree の和は五以下である。exact residual multiplicity や normalization の十分条件を必要としない。

D=V(F) の一般点 p は S 外、g1g2≠0、rank二 locus に取れる。κ=g1(p)/g2(p) とすれば p は CI の S 外の surface component Zκ に属する。dlog(g1/g2) は π_D の annihilator なので、Zκ はその一般 smooth locus 上で π の二次元 leaf であり、Zariski closure の radical ideal は Poisson stable。従って p を通る integral Poisson invariant surface の degree は五以下。

この構成は全一般 p に適用でき、例示した一枚の surface の存在だけではない。もし κ の有限の bad 値を避ける必要があれば、その逆像は proper closed subset なので一般点 p から除ける。各 true T37 の実 PGL 像と scalar に同じ被覆が成り立つ。

## 3. 全 component closure への DVR 延長

**補題。** 正則 Poisson family π_R on P4_R（R は char零 complete DVR）の generic divisor が上の dense T37 軌道の点で、special actual quintic D0 が非零とする。D0 の任意の一般 smooth rank二 point p0 は、次数五以下の integral Poisson invariant surface component に属する。

証明。actual quintic は全五Pfaffianを使う一つの degree五 homogeneous polynomial。special が非零なので relative D_R は relative Cartier divisor、DVR 上 flat。generic quintic は既約であり、D_R は generic fibre の閉包である。

p0 で D0 は smooth だから、formal smoothness により p0 を通る p_R∈D_R を選べる。generic p_K は residue pencil の base/bad locus を避けるよう選ぶ。必要な有限拡大と DVR base change を許す。第2節から p_K を含む integral Poisson surface Z_K、degree≤5 を取る。

この一枚 Z_K の Hilbert polynomial を固定し、projective Hilbert scheme の properness で flat family Z_R⊂P4_R に延長する。全 T37 source に一様な Hilbert polynomial bound を仮定する必要はない。p_R は generic fibre で Z_K に属するので closedness により special p0 も support(Z0) に属する。F_R∈I(Z_R) は generic inclusion と flatness から延び、従って Z0⊂D0。

任意の ambient local Hamiltonian derivation δ_R={f,−}_πR を取る。I_K は δ_K-stable。I_R の section g の δ_R(g) は O_ZR 内で generic に零であり、O_ZR は R-flat で uniformizer torsion がないから全体で零である。よって I0 は全 local Hamiltonian derivation に stable。

char零では stable ideal の radical も stable：g^n∈I に n 回 derivation を作用させ、radical で還元すれば n!(δg)^n=0 となる。reduced quotient の各 minimal prime へ localize すればその局所環は field なので、各 minimal prime 自身も derivation stable。したがって special support の全 irreducible component は Poisson invariant である。

p0 の Poisson rank は二。p0 を含む invariant component は local Hamiltonian flow の二次元 leaf を含むので、dimension は少なくとも二。fixed Hilbert polynomial の degree は二だから dimension は二を超えない。従って p0 は top-dimensional invariant surface component Z に属する。flat degree の総和は deg Z_K≤5 で、effective generic multiplicity を含めても deg Z≤5。embedded point 又は低次元成分だけに marked p0 が残る逃げは、rank二と invariant minimal prime の議論で排除される。証明終。

以上は Chow 上の invariant-cycle locus が閉という未証明入力を用いない。必要なのは一つの Hilbert polynomial に対する properness、flat quotient の torsionfreeness と char零の differential radical だけである。normalization の flatness・canonical 性の境界保存も仮定しない。

## 4. Whole T37 での結論

whole T37 は第1節の24軌道の閉包。rank四点 π0 がそこに属すると仮定すれば、curve selection と normalization／有限 base change により第3節の DVR family を選べる。従って π0 の actual D0 の各 component の一般 smooth rank二点は、degree≤5 の invariant surface component に覆われる。

対象 component 上で一般 leaf closure の degree が六以上、又は一般 leaf が Zariski dense なら、この必要条件に反する。ここではその target 計算を仮定せず、T37 の whole 必要条件のみを証明した。新 component 数、全 Poisson 支持の完成又は他の既知親についての主張は追加しない。
