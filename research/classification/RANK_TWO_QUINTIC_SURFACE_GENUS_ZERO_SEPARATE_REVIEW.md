# 次数五 rank 二・正規化種数零：導手分裂と全射影上界の独立監査

2026-11-07。producer は `root`、reviewer は `log_manuscript`。
`RANK_TWO_QUINTIC_SURFACE_GENUS_ZERO_NONDOMINANCE.md`
（SHA256 `50b850f09d8d6f7f44e9e7de343d761db67b83e811c460b20b3c09d3d7fe5a95`）を
別作者として読み、以下の証明により **全分岐の G 非優越命題として受理**する。
local unit 自身の Galois 不変性、normal な導手、全域で定義された involution は仮定しない。
CAS 又は数値計算は行わない。

**受理する命題。** actual D が既約被約 quintic である rank 四 Poisson tensor について、
次数五 integral 非正規曲面 S 上の一般 restriction rank が二、正規化の sectional genus が零なら、
全 PGL 飽和射影像は24次元以下である。
005、122 と複数の導手 class は空、残る113/023/014 の上界は24/23/22である。
genus 一、他次数、全 G/E、既知成分帰属又は新成分数の主張ではない。

## 1. generic node の二導手 sheet は rank 二なら分裂する

\(K=\mathbb C(S)\) とする。\(\delta=1\) の generic node の normalization conductor algebra は、
二つの K-rational branch 又は一つの quadratic separable residue-field extension の二択である。
前者は二つの大域 prime、各々 S に birational、後者は一つの prime \(\Gamma\) と
\(E=\mathbb C(\Gamma)\)、\([E:K]=2\) を与える。
char 零の quadratic extension なので field involution \(\tau\) がある。
導手 scheme 自身が normal、smooth 又は everywhere étale だという仮定は不要である。

有限分離拡大と完成後、node は \(K'[[u,v]]/(uv)\)、S は \((u,v)\) である。
Hamiltonian derivation が各 minimal prime を保存することから、全 normalization branch に regular bracket が降り、
normal parameter の ideal を保存する。
係数体を Poisson 中心と仮定せず、二つの tangential function \(a,b\in K\) を
separating transcendence basis に選ぶ。\(\pi|_S=H\partial_a\wedge\partial_b\)、\(H\ne0\) とする。

同じ dualizing generator を二 branch へ制限すると
\[
 \Omega_+=g_+\frac{dv}{v}\wedge da\wedge db,\qquad
 \Omega_-=-g_-\frac{du}{u}\wedge da\wedge db,
 \quad g_+|_{v=0}=g_-|_{u=0}=g_0\ne0.
\]
これは ambient regular unit g と \(\operatorname{Res}(du\wedge dv/(uv))\) を使った
node dualizing module の直接計算である。
mixed normal bracket は branch parameter で割れるため、
\(\alpha=i_{\pi_X}\Omega_X\) の residue に寄与しない。
二 residue は共通の非零定数 convention を除いて
\(g_0H,-g_0H\)、従って非零で互いに反対符号となる。
これは dualizing residue の **trace-zero 関係**である。

この関係は座標に依らない。uniformizer を unit 倍しても \(d\log u\) の residue は変わらない。
ambient unit g や u,v が \(\tau\) で変わっても、同じ invariant \(\alpha\) の residue は functorial に交換される。
**g 自身が \(\tau\)-固定だとは置かない。**
one-prime の場合、E 上の residue \(\rho\) は二 embedding の下で反対符号になるので
\(\tau(\rho)=-\rho\) となる。
一方 \(d\alpha=0\) と simple pole から \(d_\Gamma\rho=0\)、従って \(\rho\in\mathbb C\)。
field involution は C を固定するので \(\tau(\rho)=\rho\)、ゆえに \(2\rho=0\)。
rank 二からの非零性と矛盾する。

従って S の generic node 上の conductor は必ず二つの prime に分裂する。
これは global double cover の involution を scheme 全体へ延長する主張ではなく、function field と generic residue の議論である。
ramification、導手の codimension 二の特殊点、他の非正規成分はこの結論を変えない。
rank 零では \(H=0\) となるためこの lemma を流用しない。

## 2. 全正規化と plane conductor の分類

既受理 length 式は \(\delta_S=1\)、唯一の残余 plane P、\(\delta_P=1\) を与える。
rank に依存しない Bertini/Höring の部分から全 minimal scroll 型は113/122/014/023/005。
005 の line vertex の像では二次元非正規支持は plane のみであり、次数五 S は存在しない。
rank 零専用 cone 排除は使わない。

残る四型で standard smooth projective bundle を \(Y=\mathbb P_{\mathbb P^1}(E)\)、tautological class を L とする。
113/122 は Y=X、014/023 は zero-summand section P1 のみを点へ縮める small morphism。
exceptional divisor が無いので Weil/Pic class は \(\mathbb ZL\oplus\mathbb ZF\) で扱え、
全 conductor の strict transform は \(C=3L-3F\) である。
small exceptional curve へ追加の divisorial pole を許していない。
\(L^2C_S=10\)、\(L^2C_P=2\) である。

effective \(A=aL+bF\) の degree 二条件は \(a\ge0\)、\(b=2-5a\)。
\(H^0(\operatorname{Sym}^a E(2-5a))\) の全 summand を見ると、次表になる。

| E | 全 effective \(C_P\) | 最大 a の理由 |
|---|---|---|
|113|\(2F\)、unique \(Q=L-3F\)|a二以上は \(2-2a<0\)|
|122|\(2F\) のみ|a一以上は \(2-3a<0\)|
|014|\(2F,Q+F_b,2Q\)、unique \(Q=L-4F\)|a三以上は \(2-a<0\)|
|023|\(2F\)、unique \(Q=L-3F\)|a二以上は \(2-2a<0\)|

014 の a一は highest coordinate と base 一次式、a二はその coordinate の平方だけ。
\(2F\) の二 fiber 又は二倍 fiber も全て保持する。
各場合で residue rank on P は仮定していない。

一般平面切断は rational quintic。
S が generic cusp なら五つの transverse ordinary cusp を持つが、
degree 八 tangent Wronskian の共通五因子は dual degree 三以下、
その rational dual の bidual degree 四以下を強制し次数五に反する。
これは既受理の全 Gauss 補題であり、cuspidal α に node の留数論を課したのではない。
よって S は generic node、§1 より \(C_S=\Gamma_1+\Gamma_2\)、各 \(L^2\Gamma_i=5\) となる。

## 3. 空となる全 class と残る三型

prime \(\Gamma_i=a_iL+b_iF\) で a零なら prime vertical divisor は F だけで次数一。
従って \(a_i\ge1\)。\(C_P=2F\) なら \(C_S=3L-5F\) だから、交換を除いて
\(\Gamma_1\sim L\)、\(\Gamma_2\sim2L-5F\)。
しかし \(\operatorname{Sym}^2 E(-5)\) の全 section は次の通りである。

* 122 は零。
* 113 は highest coordinate の平方と base 一次式の積。
* 014 は highest coordinate w と \(cv+q_3(s,t)w\) の積。
* 023 は highest coordinate w と \(cv+\ell_1(s,t)w\) の積。

全て prime divisor を与えない。small map も divisor は縮めないので画像で prime へ戻ることは無い。
014 の \(C_P=2Q\) も \(C_S=L+5F\) となり、二 prime の正 L 係数の和が一になるため不可能。
よって残る全場合は113/023 の \(C_P=Q\)、014 の \(C_P=Q+F_b\) だけ。
全て \(C_P\sim L-3F\)、\(C_S=2L\)、二 prime は \(|L|\) に属する。
\(C_P\) が reduced なので P の generic transverse type も node である。

## 4. closed residue 源が全て一次元以下であること

全 codimension 一の node branch bracket は regular。
\(\Omega_Y\) の pole divisor は C だけなので、\(\alpha\) は全 conductor に沿う closed simple-pole form。
Y smooth と exceptional curve の codimension 二性から、他の pole は生じない。
各 prime 上の residue は constant、同じ全 residue を持つ二 closed forms の差は
全 divisor で regular、Hartogs と \(H^0(Y,\Omega_Y^1)=0\) から零となる。

113/023 では residue を \(\lambda_1,\lambda_2,\mu\) とする。
\(L+nF\) は十分大きな n で ample なので、一般二切断曲線上の residue theorem は
\[
 (5+2n)(\lambda_1+\lambda_2)+(2+2n)\mu=0
\]
を全 n で与える。従って \(\mu=0\)、\(\lambda_1+\lambda_2=0\)。
014 では S の二 prime の residue は §1 の局所計算により \(\lambda,-\lambda\)。
P の二 prime Q,Fb も同じ計算で \(\mu,-\mu\)（P が rank 零なら両者零）である。
\((L+nF)^2Q=1+2n\)、\((L+nF)^2F_b=1\) なので \(2n\mu=0\)。
ここでも残る全源は一つの \(\lambda\) だけである。

fixed D に対し \(\pi_D\mapsto i_{\pi_D}\Omega_D\) は generic smooth locus で injective。
ambient restriction も \(H^0(\mathbb P^4,\Lambda^2T(-5))=0\) で injective。
従って同じ projection/image D に対する全 projective Poisson fiber は高々一点。
tensor について選んだ residue ansatz だけを数えたのではなく、全 rational closed form の源の上界である。

## 5. 全 Schubert 源と actual automorphism の自由 fiber

全三型で \(K=H^0(L-C_P)\simeq H^0(3F)\) は四次元、\(H^0(L)\) は八次元。
P の二 linear equations は全 \(C_P\) に沿って消えるので、実 projection 五空間 V は
\(\dim(V\cap K)\ge2\) を満たす。
node は各 branch に一次、cusp は二次以上で消えるという局所事実もこの必要条件と整合する。

intersection 二の stratum は \(\operatorname{Gr}(2,4)\) の四、
quotient の \(\operatorname{Gr}(3,4)\) の三、lift の \(\operatorname{Hom}(\mathbb C^3,\mathbb C^2)\) の六で合計十三。
intersection 三は \(3+4+2=9\)、四は三であり、全 closed Schubert 条件の上界は十三。
finite、birational、実 conductor、Poisson descent はこの源を狭めるだけである。
113/023 は Q unique、014 の Fb だけが一つの追加 parameter を持つ。

actual polarized automorphism subgroup の dimension を全 End(E) summand から数えると
113 は \(11-1+3=13\)、023 は \(12-1+3=14\)、014 は \(14-1+3=16\)。
これらは complete L-image X 上に faithful に降りる。
条件となる unique Q、fiber choice、Schubert condition は群作用に不変である。
target frame を付けた finite birational map を固定する automorphism は function field 上で identity、
よって Y の dense open と全 Y で identity となる。
個別 tensor stabilizer 零を仮定せず、この実群の全 dimension を fiber から引ける。

|型|projection と plane 導手 choice|frame|実群|全 image bound|
|---|---:|---:|---:|---:|
|113|13|24|13|24|
|023|13|24|14|23|
|014|13+1|24|16|22|

fixed-D projective Poisson fiber の零次元性から、同じ上界が全 Poisson 像へ移る。
rank 零 cone theorem の流用、特殊 centre の捨象又は未証明の family flatness は無い。

## 6. 判定と残余

新しい nodal split lemma と全五 scroll/class/射影源を独立に確認した。
producer 命題を degree 五 rank 二・正規化種数零の **全 G 非優越**として受理する。
generic cusp、非分裂 node、plane の全 residue rank、二倍 fiber、small exceptional curve、
全特殊射影を以上の各段階で明示的に扱った。
genus 一、他次数、plane restriction の残枝、全 G/E と成分認識は別義務。
152頁 TeX/PDF、global 台帳、Git は変更しない。
