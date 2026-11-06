# 階数三接配置 T(w)=0：全源の閉包は24次元成分

2026-11-07。主課題 E_quadric3。
`RANK3_TWZERO_NEW_COMPONENT_EXISTENCE.md` の存在結論を強める。
Q=w+u² の T(w)=0 枝を覆う実既約源の PGL5 飽和閉包 Z は、
**全射影 Poisson 空間の24次元既約成分**である。
従前の59成分と相異なるので、個別に構成した成分の下界は60となる。
全 Poisson 空間の網羅性、Z の一般被約性・正規性は主張しない。

## 1. 全支持を覆う実族

`RANK3_TANGENT_ZERO_TW_RESTRICTION.md` は、全六余剰クラスを含む逆形式表示に
全 Pfaffian 条件を課し、この枝の全点を

\[
 \mathcal T=SL_2\times\mathbb G_m\times\mathbb A^7
\]

の正則像で覆っている。具体的な式は同ノート(9)–(10)と
`code/classification_rank3_twzero_source.py` にある。
源の全射性は、対称一次行列 pencil の二つの階数一係数を分解する幾何学的証明であり、
有限個の標本からの推測ではない。全点が実際に Poisson であることは、
全パラメータの閉性・Pfaffian 恒等式により確認されている。

以下では、SL2 の chart における
`a=b=c=1,δ=D=D0=S=2,C=g=T0=T=1` の点 P* を使う。
この点は r=w(2v+z+w) という非零余剰クラスを持つ。
primitive の r=0 境界だけでの接空間検査ではない。

## 2. 真の49方向と三法線

全150個の二次同次係数空間で、元の350式 `F(P)=E∧Jac(P)` の微分を D とする。
Euler-vertical・全GL・scalar方向と、全11個の実族微分を合わせた行列を S とする。
独立な整数係数演算による QQ の厳密結果は

\[
 \operatorname{rank}D=98,\quad
 \dim\ker D=52,\quad
 \operatorname{rank}S=49,\quad DS=0.
\]

実11微分を加える前の群・鉛直・scalar像は47次元である。
従って実源の二方向を落として五法線を作る必要はなく、
全実族に対する残りの法線は三方向である。

`CLASSIFICATION_RANK3_TWZERO_NORMAL_QUADRICS_20261107.json` に保存した三つの
整数二ベクトル R0,R1,R2 は

\[
 DR_i=0,\qquad
 \ker D=\operatorname{im}S\oplus\langle R_0,R_1,R_2\rangle
\]

を満たす。次に D の独立な98列を U として取り、その可逆な98行 block を M とする。
U の元の150係数方向は `im S⊕span(Ri)` と独立であり、全150方向を張る。
従って実族方向49本を零と置く**実際の101次元線形切片**を選べる。

この切片内で98個の従属座標を形式陰関数定理で消去し、真の局所環を

\[
 A=\mathbb C[[n_0,n_1,n_2]]/I
\]

で表す。残る252本の元方程式は全て保持する。
群作用が実族の二方向まで動かすという仮定や、族方向に関する形式積構造は使わない。

## 3. 三つの真正初期式で支持が零次元になる

`R(n)=n0 R0+n1 R1+n2 R2` とし、F(R(n)) の六つの二次単項式係数を350×6行列 Q にする。
上記の選択98行に添字 sel を付ければ、形式消去後の二次係数は厳密に

\[
 B=Q-U M^{-1}Q_{\mathrm{sel}}.                 \tag{1}
\]

これは全方程式の余核を取ったものであり、一部の方程式だけの近似ではない。
独立な正規化で得た QQ 行空間の基底は

\[
\begin{split}
 &n_0^2,\qquad n_0n_1,\qquad
 n_0n_2-\frac{28}{3477}n_2^2,\\
 &n_1^2,\qquad n_1n_2-\frac{7}{13908}n_2^2.
\end{split}                                                     \tag{2}
\]

各式は真の形式残余式の最初の非零斉次項であり、`in(I)` に属する。
全初期 ideal を生成するという仮定は不要である。

支持の判定には(2)のうち三式だけでよい。
`n0²=n1²=0` から被約支持では n0=n1=0。
第三式の非零係数から n2=0 も従う。
あるいは `c=28/3477`、`f=n0n2−c n2²` と置けば

\[
 n_2^3=c^{-2}\{n_2(n_0^2)-(n_0+cn_2)f\}
\]

という所属恒等式でも確認できる。
従ってこれら三初期式の根基は `(n0,n1,n2)`、`dim A=0` である。
全 Gröbner 基底、全初期 ideal、非被約長の計算を次元結論の前提にしない。

## 4. Krull次元上界と全成分性

元の全150係数局所環を49本の線形式で切ると、上の零次元切片になる。
Krull の次元不等式により、元の局所次元は高々49である。
一方、S は実際の Poisson 族と群作用の微分で階数49なので、
その既約像閉包の次元は少なくとも49。
従ってその閉包を含む全周囲成分と像閉包は同じ49次元で一致する。
24次元の Euler-vertical 商と一つの射影尺度を除けば

\[
 \dim Z=49-24-1=24.
\]

これで Z が全射影 Poisson 空間の既約成分であることが分かった。
全 T(w)=0 支持は第1節の全射な実源とその PGL 飽和により Z に含まれる。

この論法だけで P* を通る他の局所成分がないとは結論しない。
別の枝が同じ実49接空間を持つ可能性まで排除する必要はなく、
全源の閉包自身が成分であることが支持被覆には十分である。

## 5. 相異性と独立査読の範囲

次元24が確定したため、相異性には以下の短い議論で足りる。
従前59成分のうち次元24でないものとは、次元だけで相異なる。
同じ次元の既知一般 rank4 成分は MPS の8成分と CPLUS_EXT である。

全源の退化五次式は `LQ²`、rank(Q)=3である。
射影的な掛け算写像による `LQ², rank(Q)≤3` の像は閉であり、
零を含む affine cone を引き戻すことで rank2 境界も保持できる。
従って Z の**全閉包**はこの条件を満たし、rank4 部分の退化因子は非被約である。
各 MPS 成分は被約 SNC 退化因子を持つ点を含み、CPLUS_EXT の一般元は
相異なる因子の和 `H+Q1+Q2` である。従って Z はこれらの成分と一致しない。
既知 rank2 成分との比較は、必要なら閉条件 `π∧π=0` でも直接排除できる。

この24次元成分定理の相異性には、S61 graph 認識や高次元成分の因子 incidence は不要である。
`RANK3_TWZERO_NEW_COMPONENT_EXISTENCE.md` と同 REVIEW は、正確な次元を使わずに
未知成分の存在だけを示す独立な弱い証明として保持する。
そこに記録した graph 診断を、本定理の追加依存にはしない。

形式切片の論理は generic_structure 担当が独立に点検した。
特に群による積分解に依存せず、49本の線形切片と Krull の不等式を使うこと、
局所成分の一意性を主張しないことを確認した。
新しい有限二次係数は anisotropic_review 担当が別の Fraction 疎多項式実装で再構成した。
`reports/CLASSIFICATION_RANK3_TWZERO_NORMAL_QUADRICS_REVIEW_20261107.json` は、
全 D の階数98、実源49、実源＋法線52、全65列の D による零化、
源の階数上界に必要な13本の QQ 関係を確認する。
さらに五つの余核線形汎関数を独自に抽出し、各々について全150列での消滅と、
六つの二次係数への評価が(2)と一致することを直接検証した。
根基が最大 ideal であることも明示的な所属恒等式で独立確認した。
独立実行は約1.54秒で全 PASS、入力・実装 hash の一致も確認済みである。

主実行は約1.57秒、事前 gate と60秒・1.5GiB制限の下で完了した。

```sh
python3 scripts/classification.py gate research/classification/jobs/rank3_twzero_normal_quadrics.json
OPENBLAS_NUM_THREADS=1 .venv/bin/python code/classification_rank3_twzero_normal_quadrics.py --output work/rank3-twzero-normal-quadrics-replay.json
```

原 report は `reports/CLASSIFICATION_RANK3_TWZERO_NORMAL_QUADRICS_20261107.json`。
三つの物理二ベクトル、全実源49方向の選択、98方向の消去 block、全五初期式、
全係数行列の hash と入力 hash を保存している。
これにより、全法線支持の必要な上界を再現できる。
