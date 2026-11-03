# actual CA/CB odd U³ classの独立検証

`AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63` を受理する。両族でactual Iの \(J/\mathfrak m_RJ\)、\(J=I+(\nu)\)、への像は正確に10次元であり、\(\mu(I)\ge10\) である。Iの生成数上界や全ideal型は決めない。

## full sourceの再構成

既受理の全critical idealの七つのmonic消去式と \(u=U,v=w=0\) を同時に解き、全十source係数を独立に求めた。二次・三次項を捨てたcritical tangentではなく、全monic sectionである。旧critical idealをこのsectionへ制限したidealが正確に \((U^3)\) となることも確認した。

raw CA/CBの全20変数sourceを直接canonical化し、**全generic sourceを微分した後**にsectionを代入した。両族でratioはsection上0であるが、CBのratio微分列はvariable fを含んだままである。fの第三係数とHの二次・三次補正も全て残した。全150 physical出力、variable fiber20列、全flag9列、150×29 Jacobianが証明書に一致する。

## 独立なinverse三次展開

共通actual28行でのblockを \(T_0+T_1+T_2+T_3\) とし、非可換な行列の順序を保って

\[
R_0=T_0^{-1},\quad R_1=-R_0T_1R_0,
\]

\[
R_2=R_0T_1R_0T_1R_0-R_0T_2R_0,
\]

\[
R_3=-R_0T_3R_0+R_0T_1R_0T_2R_0
 +R_0T_2R_0T_1R_0-R_0T_1R_0T_1R_0T_1R_0
\]

を構成した。積が次数3まで単位行列になることを検証し、これをfree列に掛けた全四次数のvertical vectorをproducerの再帰結果と比較した。

全150 target微分は次数0、1、2で零である。選択28行とbranch row61は三次も零で、残る全150三次係数は保存値と一致した。producerを別の出力先で一度実行した証明書もbyte一致した。

## traceとactual cubic class

branch微分が次数3まで零で、\(Dt(0)=1\) だから、tracefree tはsection上で次数4以上となる。rank2表示 \(F_i=a_i(R)+b_i(R)t\) をDで微分すると、even部分は三次に寄与しない。また各b_iのsection上の次数0、1、2が零なので、Dtの非定数部分も三次に寄与しない。従って全Schur三次係数はactual b_iの三次制限そのものである。

基底1,pの変更によるtracefree tの単位倍も中心値1なので、このsectionで次数3から始まるodd係数を変えない。source ramification idealの生成元をactual odd係数に置き換える仮定をしていない。

七つのactual一次基底への全定数結合を引き、さらに既受理のWU/WV二次基底row46、95への全定数結合を引いた。全150結果を独立に照合すると、非零なのはraw row63だけで、

\[
\mathrm{CA}:\quad 12U^3/125,\qquad
\mathrm{CB}:\quad-4U^3/243
\]

となる。

## J/mJの十番目のclassと正確な範囲

このsectionでは全actual b_iの次数が3以上、\(\nu=t^2\) は次数8以上である。従って全Jの元は次数3以上となり、\(\mathfrak m_RJ\) の元は次数4以上となる。U³係数はJ/mJ上の線形functionalを与える。

新しいresidualはfull R firstjetが零であり、critical三変数tangent上の二次classも零になるように引算してある。このclassが以前の七一次＋二二次のspanに入るなら、まず一次部分の独立性により七つの係数が零、次いでWU/WVの独立性により残る二係数も零となる。しかしU³ functionalは非零なので、そのclassは零ではない。従ってactual Iの像は10次元以上である。

既知の全ideal Jは、七一次、三二次、一本の三次classを持ち、\(\dim J/mJ=11\) である。前の独立検証は \(\nu\notin I\) を全次数で証明した。もしIの像が全11次元なら、NakayamaによりI=Jとなり矛盾する。よって像は **正確に10次元** である。

自然な全射 \(I/\mathfrak m_RI\to\operatorname{image}(I\to J/mJ)\) から \(\mu(I)\ge10\) が従う。一方、\(I\cap\mathfrak m_RJ=\mathfrak m_RI\) は未証明なので、\(\mu(I)=10\) や全十生成表示を結論しない。cubic classを含むleading dataを全formal idealと同一視しない。

本稿はdepth、normalizationのsubintegrality、actual D、あるいはU-axisのanalytic cusp型を決めない。R19\Nを個々のCA/CBのconductor枝として使う議論も含まない。

独立codeは `code/ambient_B18_CA_CB_odd_cubic_review_20261103_v63.py`、receiptは `reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63_REVIEW.json` である。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/ambient_B18_CA_CB_odd_cubic_review_20261103_v63.py \
  --replay work/B18-CA-CB-cubic-producer-replay \
  --output work/B18-CA-CB-cubic-independent-review.json
```
