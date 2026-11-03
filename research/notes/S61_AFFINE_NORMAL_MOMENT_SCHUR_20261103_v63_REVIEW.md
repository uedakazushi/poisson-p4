# Controlled moment Schurの532/536構成と135式の独立監査

2026-11-03。凍結した `S61_AFFINE_NORMAL_MOMENT_SCHUR` と、新しい `...SHARP_MOMENT_SCHUR` の両方を承認する。
新しい532構成は旧536構成を変更せず改善する。**一つのscalar unit判定は正確に未知残余Bの空性と同値**である。
実際のS61 unitやdeltaの巨大展開は計算していない。

## 1. 各未知fiberでの全sourceとnorm

先行独立reviewが検査した実5×670表示とR'→Bの全射を使う。
未知幾何点pではC_p=k[h]/G、A tensor k(p)=kなので、全134 remainderのidealはJ=(h-h0)である。
G=(h-h0)g4とし、r_i=(h-h0)s_i（deg s_i≤3）と書く。
D=k[h]/g4からJへのf↦(h-h0)fは同型であり、s_iはDのunit idealを生成する。
GやDが非被約でもこの議論は有効である。

r(a)=sum a^i r_i、S(a,h)=sum a^i s_iを取る。
image basis `(h-h0)h^j,j=0..3` からbottom四係数への変換はtriangular determinant1である。
独立コードは一般g4,h0について5次・4次companion行列のintertwiningを全20entryで確認した。
よって新四列r,hr,h²r,h³rのminorは厳密に

    delta0(a)=Norm_D(S(a,h)).

旧四列(h-a)r,(h-a)hr,...のminorは `g4(a)*Norm_D(S)` である。
これはdeterminantの乗法性による普遍恒等式であり、fixtureからの外挿ではない。
独立fixtureはproducerのremainder/resultant算法に対し、companion powersと4次multiplication normを使った。
simple、多重根、完全ramifiedの全行列が一致した。G=h^5、r0=h²、r133=hでは新delta0=a^532、旧delta=a^536となる。

## 2. 全未知点で非零であること

Artin代数Dの各local factorでは、s_iの少なくとも一つがunitである。
そのresidueにおけるS(a)は非零polynomial（次数≤133）である。
generic aでSは全local factorのunitとなり、normは非零になる。
4次normのa次数は≤4*133=532。旧式ではg4(a)により更に4を加えて536である。
monic degree4 g4がmultiplication rのannihilatorを張るので、degree≤3 domainは常に補空間であり、
新式に(h-a)を付ける必要はない。

全未知幾何点でdelta0(a,p)が非零polynomialだから、その全係数が生成するidealはBでunitである。
これは自由R'上でdeltaがunitという意味ではない。
primitive polynomialはAss(B[a])の各p[a]を避けるため非零因子であり、既存のAssとnilpotent submoduleを局所化で失わない。

さらに **B→B(a)[delta0^-1]は忠実平坦**である。
polynomial extensionと二つのlocalizationは平坦であり、各p∈Spec B上のfiberは
`k(p)[a]` を非零K[a]と非零delta0_p(a)で局所化した非零環である。
従って全fiberが非空で、忠実平坦性が従う。
これはB(a)の全点を一chartが覆う主張ではない。B(a)→そのdelta局所化は一般に忠実平坦ではなく、
delta=0にある閉点を失い得る。しかし元のBに関する空性・ideal・module情報はこの忠実平坦なmapで検出できる。

known E5では全Phi entryが零なのでdelta0=0である。従ってA(a)_delta0=B(a)_delta0となり、
この一つのscalar chartが零環であることとB=0は正確に同値である。

## 3. 全点coverと有限sourceへの対応

新四列は各jで `sum_i a^i Phi_column(5i+j)`、旧四列は
`sum_i a^i(Phi_column(5i+j+1)-a Phi_column(5i+j))` であり、全134sourceを保持する。
全indexと有限sparse recipesを別証拠として保存した。actual deltaを展開したとは主張しない。

三つの独立generic aを使えば、全三deltaが零になるincidenceはB上にfinite fibersを持つので次元≤2であり、
三次元constant空間のgeneric点を含まない。従って全Bを三chartで覆う。
また新式では0,...,532の533固定QQ値、旧式では537値で全点を覆う。
Vandermonde補間により係数idealと全評価値idealは同じなので、nilpotentsを含むschemeの開被覆である。
固定値coverにはdim Bの上界を必要としない。

## 4. 135本の実algebra表示

M_rはCでのmultiplication rの5×5行列、Pは先頭四列のbottom block、delta0=det Pとする。
delta0を逆にして先頭四relationを消すと、target 1,h,...,h⁴の後ろ四つは1の倍数となる。
特に `h=tau=U/delta0`、`U=-(toprow*adj(P))_0` であり、最後のrelationはdet M_r/delta0である。
1が実unit generatorなので、moduleとしてだけでなくunital algebraとして

    (C/(r))_delta0 = R'_delta0/(N),  N=det M_r

が成立する。r自体が元134の結合なので、元Aはこの環へ全r_i(tau)を加えた商である。
各r_iのh次数≤4だからdelta0^4で分母を払えば、Nと134本のclearした式の計135本になる。
全phi entryのformal次数≤2より、N≤10、Uとdelta0≤8、clearした式≤2+4*8=34である。
先行670本degree≤10の選択肢も同じ実schemeを表す。

独立coefficient-ring fixture G=h^5+c*h+b,r=h-xでは、companionから直接delta0=1、tau=x、
N=-(x^5+c*x+b)を確認した。実S61 scalar unitの標本ではない。

このpacketは空性判定の有限recipeとそのscheme的正確さを証明する。
actual unit、scalar cofactor、全w飽和、一般五点尽くし自体は未完了であり、新しいmodel rank・GB・Macaulayは実行していない。
review code: `code/s61_affine_normal_moment_schur_review_20261103_v63.py`。
receipt: `reports/S61_AFFINE_NORMAL_MOMENT_SCHUR_20261103_v63_REVIEW/`。
