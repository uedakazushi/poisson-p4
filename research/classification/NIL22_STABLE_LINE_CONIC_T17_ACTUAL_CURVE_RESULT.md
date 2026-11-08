# Stable line-conic 全六磁気源の T17 actual confluence：有限恒等式結果

2026-11-08、PRIMARY E_multiple_hyperplane。
純証明 `NIL22_STABLE_LINE_CONIC_T17_ACTUAL_CURVE.md` の bytes を変更せず、
同一全変数式の ONE finite readback 結果を記録する。
状態は **producer 全恒等式 PASS、root 独立手査読 PASS、正式独立 receipt 待ち**。
新 component count は零。

## 1. 完全全源の検証結果

job `jobs/nil22_stable_line_conic_t17_actual_curve.json` は実行直前 gate PASS。
60 秒・1 GiB の上限で、標準ライブラリ Fraction sparse Laurent の一括全変数計算を実行した。
実測 0.192733649 秒、最大 RSS 13,952 KiB。
ε と全 m0,…,m5 は同時に保持し、source/helper import、別点、Jacobi operator/rank/minor、
profile、jet、local ideal の計算は無い。実行 source 全 bytes と hash は report に埋込み済み。

全て PASS：

- true completed T17 全三 period・二 edge の正逆と k=s+r、t1=nu*k の全積関係。
- raw completed representative と H-Casimir gauge の全 Euler 差。
- original completion の t0*t1/p1 を formal kappa で文字通り計算し、分子内で相殺後に kappa=k を代入。
- 全150係数の forward/inverse と元 raw completed representative の厳密な復元。
- actual GL5 の両 inverse、det=-1/(9ε²)、tensor scale ε。
- 全負 ε 冪の相殺、全括弧の ε 次数高々二。
- 全任意六磁気端点が公開元 stable_line_conic 全150係数に一致。
- 全六磁気と四 true period derivatives＋二有限 shift の polynomial inverse、det=-2916。
- 全 strict10 Jacobi、全 projective350 Jacobi 係数の恒等零。
- signed Euler Pfaffian 全126係数の変換則、元端点実 D の厳密な一致。

実端点 D は `-h²*(x-y-12h)*(36h²+3h(x-y)+2xy)/2`。
磁気 parameter は一つも反転しない。全零、magnetic invariant 零、係数零境界も含む。
formal kappa は完全に消去され、実 curve は QQ[ε,m0,…,m5,h,Q1,Q2,x,y] 上の多項式。

## 2. 正しい assignment と whole membership

本 curve は F→smooth Q、G→hL の assignment である。
以前の F→hL,G→Q 固定 plane 診断は別の class の陰性で、本曲線への障害ではない。
中心 coisotropy の全式と true period derivative の E=εE1+ε²E2 が、
full tensor の regularity と任意 full magnetic endpoint を同時に与える。

各任意に固定した全六 m に対し、k(0)=2、t(0)=81/2、nu(0)=-1/2、l=1。
ε!=0 と必要 source open は非空稠密。
そこで true completed T17 source の actual PGL 像と非零 tensor scale なので T17 内。
全多項式 tensor は非零 rank4 の元 endpoint へ特殊化するから、named component の閉性により
**元 stable_line_conic 全六磁気源、全 PGL 飽和とその閉包が T17**。
local smoothness、unique component、ambient dimension、形式 higher support は不要。

root は実行前に coisotropy、全 intermediate 六括弧、E1/E2、四磁気 map と全逆、
二有限 shift、physical F/Q assignment を独立手査読した。
本 producer 結果と root の独立手査読は区別し、正式 receipt は別ファイルとして保持する。
既 three-lines 全六磁気→T03 と合わせれば、stable 二源の帰属義務は両方陽性になる。

## 3. 再生と metadata 履歴

report の `executed_source_bytes` は job の inline source SHA と一致する。
source は出力 report を exclusive create するため、再生時は既存 bytes を上書きしない別作業領域を使う。
source parameter/frame と全多項式 family/D は report に保存した。
monomial 順は `(ε,m0,m1,m2,m3,m4,m5,h,Q1,Q2,x,y,kappa)`。
最後の kappa 指数は全項零であり、temporary parameter は target に残らない。

初回 contract の hypotheses を list で書いた metadata gate failure は、
`inputs/nil22_stable_line_conic_t17_contract_metadata_history_20261108.json` に exact bytes 保存した。
その時点で mathematical run は一度も始まらず、nonempty string への訂正後 gate PASS、
root GO、fresh gate PASS の後に本同一全変数 batch を一度だけ実行した。
数学的失敗を隠した再点計算ではない。
