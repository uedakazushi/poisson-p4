# 停止済み prefix の出力予測監査と、別方式の数学的条件

2026-11-04。これは source/forecast の読取りだけである。凍結した24 pivot と
124840項の212×59残余を変更しない。停止した行23・列65の行更新も、別の pivot
も、形式級数計算も行わない。

## 旧37.83MiB予測の内訳

旧 worker は各 exponent に符号込み22文字を一律計上し、有理数の和の分母を積で
見積もっていた。実残余の content 分母は最大92bit、primitive 整数は最大89bit、
実 exponent は H≤9,L≤4,s0/q0/t0≤5 である。次の pivot は QQ 定数2048/27。
この metadata から同じ旧予測37830400 bytesと係数146462項上界を再構築する。
旧予測が誤りなのではなく、非常に保守的だったことを区別する。

## 新しい、証明付き byte 上界

各既存多項式を content*c の形 c*Z、Z primitive 整数多項式として読む。
factor cf*Zf、正規化 pivot entry cp*Zp に対して、係数を一つも掛けずに
product の共通 scalar content cf*cp は正確に求まる。その整数部分の係数は
min(#Zf,#Zp)*max|Zf|*max|Zp| で抑えられる。

old=ca*Za との和には

```
D=lcm(den(ca),den(cf*cp))
M=|num(ca)|*(D/den(ca))*max|Za|
  +|num(cf*cp)|*(D/den(cf*cp))*min(#Zf,#Zp)*max|Zf|*max|Zp|
```

を使う。Dは scalar整数のLCMだけで、未知の多項式 GCD や field 除算ではない。
全結果係数は整数 H_i/D、|H_i|≤M と表せる。g=gcd(H_i) とすれば

```
digits(g)+sum_i digits(H_i/g) ≤ n*digits(M)+1
```

である。従って g 自体や新しい係数を計算せず、元行列と同じ primitive-content
serialization の全数値文字数を厳密に抑えられる。分母は一度だけ記録し、
各 exponent の文字数は support 和集合から実際に数える。残余の非更新成分、
factor、正規化行についても既存値と scalar rescaling の文字数を用いる。
gzipでRAW breachを救う議論ではない。

同じ停止候補について旧・新 bound の内訳を証拠へ保存する。他の定数／単項式
候補も項積上界順の最初の少数だけを support/height 予測で比較する。選択順を
変えて計算したわけではなく、予測の一覧は新しい rank/minor の証拠でもない。

## 全 worker の係数と scratch

byte bound を改善しても、150000項上界を緩めない。次の候補の146462項予測は
旧実装と同じ全 pivot support forecast であり、正規化を行う前に全更新を確認する。
L+2k+3b+16には product、shiftを揃えるコピー、結果、back-addの全所有係数を含む。
新しい整数 LCM/height bound は一成分ずつ求めて破棄し、全更新の大整数配列を
常駐させない。保存するのは項数・bit数・byte数だけとする。

仮に継続 worker が承認される場合も、入力は現在の残余だけを streaming で読み、
旧49039項行列や24段の値コピーを同時保持しない。各次 pivot 全体を事前予測し、
初めて不適合となる前に停止する帰納的 guard が、worker 全体の上界になる。
この予測は残る59 pivot の完了や全途中状態の小ささを保証しない。

## state-reference receipt

行操作を

```
pivot=(sourceRow,sourceCol,QQcoefficient,parameterExponent)
normalized row[j] = previous[row,j]/pivot
new[i,j] = previous[i,j] − previous[i,col]*normalized row[j]
```

として段番号付き previous-state 参照で記録できる。初期 state は凍結残余の
hashで指定し、全 active original row/column IDs と target IDsを保持する。
factor と正規化行の同じ多項式を二重に保存せず、独立 replay は初期行列から
全更新・最終残余hashを照合できる。これは参照形式の提案のみで、まだ新しい
prefix worker や新しい row-operation coefficient を作っていない。

## 形式 Laurent 級数方式について

QQ[H,L,s0,q0,t0]を一つのTの正整数冪へ送る写像は一般に単射ではない。
ただし前段の全 specialized denominator が非零と証明済みなら、specialized
pivotの非零係数は実K5 pivotの非零性を証明する。specialized零から実 entryを
捨ててはいけない。実行列の行操作は全五変数の同じ previous-state program として
保持し、級数側の零性で項を省かない必要がある。

区間 f=f_known+O(T^N) を扱い、denominatorのvaluation vが既知なら inverse の
保証精度はN−2v。積の誤差精度はmin(Nf+vg,Ng+vf)、商は
min(Nf−vg,Ng+vf−2vg)で伝搬させる。未確認の leading term を零や単元として
扱わず、全59 pivotについて保証精度内の非零係数が得られた場合だけ、実lazy
programの具体的 pivot profile が閉じる。

実残余の一成分の各変数次数上界は(9,4,5,5,5)。59minorの各次数は高々
(531,236,295,295,295)。例えば基数532のKronecker重みならこの範囲の多項式
は衝突せず、既知の非零59minorも非零に保つ。しかし必要になり得る絶対T次数は
約532^5であり、これをそのままdense精度の予算とはできない。

一方、全212×59状態で各entry高々8個の保証済み係数を保存する設計なら最大
100064個であり、入力の各多項式をstreamingして破棄し、少量の単元逆・畳込み
scratchを含め150000以内に設計し得る。ただし8個の相対係数で59段を証明できる
保証は現在ない。精度損失が起きれば停止する必要があり、追加精度や別重みの
自動 retry は正当化されない。この方式の source/全状態数学条件だけを評価し、
級数 worker は実行しない。
