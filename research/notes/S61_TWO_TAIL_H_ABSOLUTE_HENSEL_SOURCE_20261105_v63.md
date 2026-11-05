# 確定Hの絶対既約性：同じd1・primeで六根seriesを使う有限証拠のsource

2026-11-05。SOURCE inventoryのみ。因子、根、Hensel/Newton係数、recombination値は実行しない。
確定した94項H以外のphysical/beta sourceを使わず、既存d=1、prime32003を固定する。
後述のnative/pure-field workerは、freshコード/custody監査と一回GOの前には実行しない。

## 特殊化とcontentの別々の条件

F(a,e)=H(a,e,1)はa次数6、e次数高々9で、a先頭係数は
`alpha*L(e)`、alpha=17578125、L=55+66e。両方のd1・mod32003制限でa次数6を保持する。
H自身のa-content unitは既受理だが、それだけでFのmodp contentがunitとは言わない。
Fのcontentはa先頭係数の線形因子Lを割るので、a定数係数P(e)を6e+5で
formal pseudo-divisionし、非零かつmod32003で非零の定数余りをsourceで保存する。
これは新しいphysical sampleではなく既知多項式の係数contentの恒等式である。
これによりFp[e]でも代数閉包でもFはaについてprimitiveとなる。

alphaをfield単元としてFbar=F/alphaとし、a先頭係数を正確にLとする。
e=0のmonic六次は、既受理Rabin/微分報告のf0=Fbar(a,0)/55と同じもの。
既存証拠からf0はsquarefreeで、四つのFp線形因子と一つの既約二次因子を持つ。

## 六根と一つの有限体

f0から既知のdegree4 gcdをexact divisionして二次qを作る。
四次の四根を固定Fp内部の有限根探索で取り出す方式なら、全32003 field要素に
固定四次Hornerを適用しても128012 base productsで完了する。
これは固定多項式の内部因数分解であって、新しいe/d/primeやbeta/physical oracleを呼ばない。
代替のinstalled nmod_poly_factor APIは利用可能だが、本有限設計はその完全性契約に依存しない。
四線形因子とqの全product、相異性、qのdiscriminant非平方を確認してseedを完全に検証する。

K=Fp[w]/q(w)を使い、四Fp根とw,-q1-wの六根を固定する。qの既約性によりKは体。
H_axis四次商環とは別の証明用fieldであり、混同しない。
六つのFbar_a(r_i(0),0)を非零unitとして逆証明を保存し、n=1..9について

```
r_i,n = - coefficient_e^n Fbar(r_i,<n(e),e) / Fbar_a(r_i(0),0)
```

で六根seriesを構成する。逆にするのはこの六つの固定field要素だけ。
六つのroot identityと `L*product_i(a-r_i(e))=Fbar mod e^10` を全係数で保存する。

## 62組が絶対因子を尽くす理由

代数閉包上の非自明な因子Gもa-contentがなく、正のa次数を持つ。
そのa先頭係数は定数倍の1またはLしかない。定数倍を調整し、先頭をL^b (b=0,1)、
補因子の先頭をL^(1-b)にする。モノミック化した因子は六根seriesの一意な部分集合積である。
従って係数はKに属する。deg_e(G),deg_e(Q)はそれぞれ高々9なので、e^10までのseriesから
元の多項式を完全に復元できる。formal seriesの存在だけで多項式因子だとは推論しない。

seed0を含むproper subsetは31個。その各々についてb=0,1の二leading allocation、計62組を作る。
G=trunc_e<10(L^b*subsetproduct)、Q=trunc_e<10(L^(1-b)*complementproduct)とし、
e次数18、a次数6までのFULL polynomial product GQをFbarと比較する。
一つも一致しなければKより大きな代数閉包にも因子がないのでFのmodp像は絶対既約。
一組でも一致した場合は、固定特殊化が可約という結果に留め、char0因子とは呼ばない。
別prime/d値や無制限factor searchを追加しない。

## char0と全Hへの帰還

Fbarのmodp絶対既約性からchar0の絶対既約性を得る際は、a-content unitとa次数6保持を使う。
もし代数的数係数の二因子があれば、p上の付値環でGauss primitiveに取り直せる。
a先頭係数の積が還元後も非零なので、両因子の正のa次数も還元で残り、絶対既約性に矛盾する。
単なる一変数Fp既約性とgeometric既約性を取り違えない。

Hのglobal a-content unitもgeometricである。先頭係数はd^3*(80-25d+66ed)で、
後者は代数閉包上でもprimitive linear polynomial、a定数係数の既知二つの非零判定が
dとこのlinear因子を除く。従ってHの代数閉包上の因子も正のa次数を持ち、d1で次数を保つ。
Fがchar0絶対既約ならHも絶対既約となる。

## 可視資源と証拠形式のgate

二次field積はKaratsubaと二次関係の還元を用いbase products五回以内。
六根Newtonの全n<=9評価80640 base products、六root最終check21000、
全64subset積378000、62full products496000、leading allocation29760、
固定四次root探索128012、残りseed/逆/normalize/checkを含め全1.3百万以下。
fresh workerではhard2百万base products、field一時値30 bit、owned24000、RAW4MiB、
45秒soft/60秒hard/2GBを候補とする。64subsetを全てdense70 K-slotで保持しても8960base cells。
一度に一組のG/Q/fullproductだけを作り、62組分を同時保持しない。
guard prefixには全根/subset・current candidate・loop cursor/pending field productを保存する。

independent checkerは全seed product/二次非平方、六固定inverse、全root jets、
subset identitiesと62個の失敗coefficientを有限積だけで検証できる。
「factor APIが一因子を返した」ことを証拠の代わりにしない。
ここではコード実行前sourceまでであり、実装の完全なowned/cursor/custodyを別に固定してからGOを要求する。
