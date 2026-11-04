# Actual N-normalizer pairの独立監査：regular rootsとall-order hypersurface

2026-11-04。author runnerを実行せず、独自のLaurent係数辞書と17単項式quotient jetを用いる。
actual E4二次係数は既存Conn独立レビューを依存gateとし、新しい高次source計算は行わない。

## 1. 物理sourceとcached arithmetic

A,B,gと全16H係数を独立symbolとして、h→(h-Bx)/Aの変換をLaurent多項式で再構築する。
元の5x5 antisymmetric tensorを実際にpushし、Euler verticalを加え、CB表示の全25entry・全150quadratic係数と照合する。
反転するのはAだけ。B18の交換済みsourceではA=-1がunitであり、conductor、B、gは反転しない。

N28のcached source2を独自の疎17-monomial環へ入れ、全28unit frame式と全150physical outputを再評価する。
同じR0九関数のmixed2jetを使って、各族のnuと全12even/odd係数をpullbackする。
原点で消えるconductor係数をunit扱いする演算は無い。

row18のconductor pullbackの一次項は、CAで2X108、CBで-X108/2。
nuのp_N²、X_j*p_N成分と物理numeratorから、符号を含めた根の一次項を独立に得る。

```
tauA=-X112/8+X117/6-p_N/2+O(m²),
tauB=-X112/3+8X117/3+2p_N+O(m²).
```

全12outputでmixed係数が一致し、source-node/R0 bindingsも元pure48 sourceと照合した。
rank8のminorは独自QQ消去で-1/2、唯一のlinear normalは

```
5X112-20X117+24tauA+6tauB.
```

## 2. 零weight N射とregularityの論理

N tailのsourceからzero-weight出力だけを取る射は、有限-w polynomial模型と非負weightのtorus作用で正当化される。
lambda=0への作用はambient上の多項式mapとして延長し、閉でGm-stableなNのidealを保つ。
したがって全formal parameterにわたり同じNzero→CAzero,CBzero射である。
完成環にw=1を直接代入してはいない。

非零のb18はnormal source B_N内の非零函数なので、まずそのfraction fieldで
eta=(q18-a18)/b18を定義できる。この開ではconductorを反転した成分と正規化が一致する。
eta²=nuというmonic relationがあるため、etaはB_N上integral。
B_Nはregular8、従って整閉であり、etaは実際にB_Nのregular functionとなる。
非優越N→CA/CBに正規化の普遍性を直接適用した議論でも、b18を中心でunit扱いした議論でもない。

全12のphysical equalityは非零b18の開で成立し、両辺がregularなB_N内にあるため全体へ延長する。
このsame-source extensionにより二つのmonic式も全output差も全次数で零。

## 3. All-order L因子と四変数coreへの正確な帰結

monic二式を既存unit blockで消したambient S9→B_Nはrank8。
formal IFTによりこのmapは全射で、regular8 imageのkernelは一つのregular parameter L_formalとなる。
その一次項が上で計算したnormalである。
したがって全12のactual equationは全次数でL_formalに可除である。二次factorだけからの外挿ではない。

既存のfull quadratic initialから、quotientのlinear spanは五つの独立座標。
S9のformal座標を(u5,core4)へ取り直せば、全idealは

```
I_E4=L_formal*(u1,...,u5,h1(core4),...,h7(core4))
```

と書ける。五つのL_formal*u_i関係はE4の中に残す。
四変数core ideal、関数係数、socle、E4全typeや元D0のdepthはまだ未計算。

独立code `code/ambient_B18_N_normalizer_pair_review_20261104_v63.py`。
Report/folder `AMBIENT_B18_N_NORMALIZER_PAIR_PREFLIGHT_20261104_v63_REVIEW`。
