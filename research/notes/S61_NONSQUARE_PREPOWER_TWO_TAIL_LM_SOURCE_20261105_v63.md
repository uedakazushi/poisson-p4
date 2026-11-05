# 保存支持による二tail非割切り：反転因子を六つに固定したSOURCE

2026-11-05。保存済み4500項N0と2317項N2の支持だけを読み、新しいtarget積・因子探索・GCDを行わない。
\(Q_0=N_0/D=T_0/(a^4C^4)\)、\(Q_2=N_2/D=T_1/(a^3C^4B^4)\) で、D197は同じ。

## まず完全な局所化リストを固定する

受理済み `GLOBAL_FULL57_SOURCE/unit-inventory.json`（SHA `75c20c39…`）と独立review6640…は、
base ringを文字通り
\[
R=\mathbf Q[a,e,d,1/(adFBCJ)]
\]
と定め、factor orderを **a,d,F,B,C,Jだけ** としている。
元 `RAW_H_B_35_MINORS/certificate.json` の FULLU_product49係数も、この六因子の積と全係数一致させる。
同ファイルの古いrank/GB主張は今回取り込まない。

Delta407,Saux,N6,N7bar,P,h22,p,R434,κ3,A_anchorはinventoryで明示的にNOT_inverted。
T,D12は後続十分開のcorollary localizerだけで、今回の環には反転しない。
**本比較をR[1/D12]や旧Ω環へ拡張しない。** さらに反転すれば以下のstrictnessが消える可能性がある。

F,Cはaに関する非零定数先頭の一次式、Jはeに関する定数先頭一次式。
Bはdに関するprimitive一次式（Rの因子3e−2に対するb1の剰余50a≠0）。
a,dと合わせ、六因子は既約で相互に非associate。literal scalar normal formも相異なる。
同じ議論はCへの基礎拡大でも成り立つ。

## 全反転因子に対するN2の非可除性

a,dは保存支持の最小指数0から。F,Jは独立受理済みvF(Q1)=vJ(Q1)=0を
\(Q_1=B^4N_2/D\) に適用する。Cも正確なvC(T1)=4と
\(T_1=a^3C^4Q_1\)、C∤aからvC(Q1)=0、従ってC∤N2。
Bは正確なvB(Q1)=4と全B⁴乗算戻し恒等式からvB(N2)=0。
こうして反転される六つをすべて覆う。D12や他のchart因子の非可除性を仮定しない。

## 小さなleading-monomial証明

固定lex a>e>dで保存支持から
\[
\operatorname{LM}(N_0)=a^{16}e^{10}d^{18},\qquad
\operatorname{LM}(N_2)=a^{12}e^{11}d^{14}.
\]
商指数のe成分は−1なので、N2はN0を通常の多項式環で割らない。
もし上記の六因子環RでN0=N2 H/sなら、sは六因子の積であり、N2は各因子と互いに素。
UFDで消去すれば通常のN2|N0となり矛盾する。
従ってこのRでは \((Q_2)\subsetneq(Q_0,Q_2)\)。

これは二tail scheme idealの比較だけである。radicalや開集合が厳密に変わるとは言えず、
unit ideal、actual bad Fitting等号、補集合消滅、全FULLU57/239も主張しない。
