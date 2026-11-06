# nil22半追加枝：全源の超平面分裂による被覆

2026-11-07。E_multiple_hyperplane。
`MODULAR_NIL22_TOP_SUPPORT.md` の `A=B=−1/2,n8=−4un9` 全源を扱う。
以下の全パラメータ実Poisson pencilにより、この追加枝の全 rank四支持は
既済の三相異超平面＋既約二次因子層の有限閉被覆に入る。
新しい成分を追加する証明ではない。

## 1. 全七変数源の正規化と補集合

u=n10,b=n6,c=n9,s=n5 と置く。s=0 は、既に全点T01へ被覆した半型 n5zero源との交わりである。
従って以下では s≠0 を扱い、補集合を捨てない。

\[
 W=w-2u,\quad Z=z-2uy,\quad
 a=n_0+2ub,\quad e=n_7-4u^2c,\quad
 d=n_2+2ub-24usc
\]

および

\[
 \delta=-2c,\quad e_0=-2b+4cs,\quad j=-(d+2ce)/(2s)
\]

を用いて、行列式一の affine 変更

\[
 X=x+\delta Z+e_0W,\quad Y=y+\delta W,
 \quad T=Z+jW-2e
\]

を行う。元の全整数化源の係数は、a′=a+2ce として

\[
\begin{aligned}
 P^{XY}&=2a',&P^{XT}&=-YT+4sYW-2sT,&P^{XW}&=W(2s-Y),\\
 P^{YT}&=-T,&P^{YW}&=-W,&P^{TW}&=-W^2
\end{aligned}                                                    \tag{1}
\]

となる。全原係数との比較は `CLASSIFICATION_NIL22_HALF_EXTRA_CANONICAL_20261107.json`
で恒等的に確認した（0.455秒）。古い凍結 witness は正確に a′=−13995,s=3 へ移り、
新しい標本を選んでいない。
以下では a′をa、X,Y,T,Wをx,y,z,wと書く。
全二変数 a,s を実際の源として保持し、一軌道とは主張しない。

## 2. 反復因子を分裂する全Poisson pencil

任意の a,s,ε に対して

\[
\begin{aligned}
 P_\epsilon^{xy}&=2a(1+\epsilon y),&
 P_\epsilon^{xz}&=-yz+4syw-2sz,&
 P_\epsilon^{xw}&=w(2s-y),\\
 P_\epsilon^{yz}&=-z(1+\epsilon y),&
 P_\epsilon^{yw}&=-w(1+\epsilon y),&
 P_\epsilon^{zw}&=-w^2-\epsilon zw
\end{aligned}                                                    \tag{2}
\]

と定義する。全係数は元の affine 座標で次数高々二であり、ε=0 は (1) の全源である。
(2) が全パラメータで Poisson であることは、次の簡単な場の計算から分かる。
w≠0 の補助座標 r=z/w において

\[
\begin{aligned}
 U_\epsilon&=(1+\epsilon y)\partial_y+(1+\epsilon r)\partial_r,
 &W&=w\partial_w,\\
 R_\epsilon&=-U_\epsilon\wedge W,
 &Z_\epsilon&=2a(1+\epsilon y)\partial_y
       +4s(y-r)\partial_r+(2s-y)W
\end{aligned}
\]

と置くと、(2) は `∂x∧Zε+Rε` と一致する。
Uε と W は可換なので Rε は Poisson。
さらに `[Zε,W]=0,[Zε,Uε]=(1+εy)W` だから `L_ZεRε=0` となる。
全場は x に依存しないので、これで (2) の全 Jacobi が零となる。
これは補助座標による恒等式の証明であり、r=z/w を PGL 作用とみなしていない。
元の係数が多項式なので、この Poisson 恒等式は w=0 にも延びる。

## 3. 実因子と全源の閉被覆

h で次数二に斉次化して h を Casimir にすると、厳密に

\[
 \operatorname{Pf}(P_\epsilon)=-2w(h+\epsilon y)q_\epsilon,
 \quad q_\epsilon=a h w+(a\epsilon+2s)hz-2syw.                \tag{3}
\]

従って実射影五次因子は

\[
 \{h=0\}+\{w=0\}+\{h+\epsilon y=0\}+\{q_\epsilon=0\}.
\]

qε の `(h,y,z,w)` における対称行列式は `s²(aε+2s)²/4`。
各固定 s≠0 の元源点について、ε≠0かつaε+2s≠0は稠密開であり、
三超平面は相異なり独立、qεはrank四の既約二次式となる。
既済の全3H＋既約Q被覆定理を適用すると、この稠密な実変形は既知の有限閉成分の和に入る。
閉性からε=0も同じ有限和に入り、(1)の全s≠0源が被覆される。
最後に§1のs=0源は既済T01へ戻るので、元の全半追加枝を尽くす。
同じ成分への一点からの移送や、元点のsmoothnessを使っていない。

`CLASSIFICATION_NIL22_HALF_EXTRA_DEFORMATION_20261107.json` は全(a,s,ε)のstrict Jacobi、
全元源との端点一致、次数二の斉次化、Pfaffian、rank四minorと場の保存を直接検算する（0.162秒）。
producerは共有した `code/classification_nil31_long_b_deformation.py --half-extra`、
有限契約は `jobs/nil22_half_extra_deformation.json`。
longB用の旧producer完全bytesは旧longB証明書に保存されており、以前の証明書とjobは変更しない。
例外層担当が補助座標の全括弧、三つのLie括弧、元の六係数とPfaffian、
三超平面の独立性とrank四条件、全支持定理の適用および両補集合を独立に手確認して受理した。

参考として、同じ凍結点で全D100、純作用接像46、実源接像47と真の三法線を保持した
重み検査は、三法線すべて重み零を返した。
`CLASSIFICATION_NIL22_HALF_EXTRA_WEIGHTED_SLICE_20261107.json` はその限定的結果だけを保存する。
本証明はその後に形式次数を増やしたものではなく、全パラメータの実pencil (2) を用いる。
被覆に不要となった三法線二次障害計算は実行していない。
