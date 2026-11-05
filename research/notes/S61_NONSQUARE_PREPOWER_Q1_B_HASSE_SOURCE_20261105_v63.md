# 残る既知単元 B：未知の b1 chart を導入しない SOURCE

2026-11-05。Jの独立d-Horner受理後のSOURCEのみ。Bの層・worker・商・GCDは未計算。

\[
B=b_1d-R,\quad b_1=50a+18e^2-9e-2=50a+(3e-2)(6e+1),\quad R=40-60e.
\]

**b1がFULLUの単元とは仮定しない。** 係数b1とRのgcdは1である。
Rの唯一の非定数因子3e−2に対しb1の剰余は50aで零でない。
Bはdに関してprimitive、従って普通のQQ[a,e,d]で既約である。
QQ(a,e)[d]でのd根R/b1のHasse付値は、その普通のB付値に一致する。
これは全係数の因子重複度の主張で、b1=0のFULLU点を捨てる操作ではない。

さらに全係数恒等式
\[
20b_1+R(6e+1)=1000a
\]
がある。整数行列
\[
M=\begin{pmatrix}b_1&-R\\6e+1&20\end{pmatrix}
\]
のdetは既知unit1000a。
adjugate置換 \(x=20X+RY,\ y=-(6e+1)X+b_1Y\) では
\(b_1x-Ry=1000aX\)。従ってFULLU基底上のhomogeneous因子表示でも未知のb1反転は不要。

## 小さい固定 root cache

Q1の原4811整数I、共通D197を保持する。
\[
S_h=R^h b_1^{18-h},\qquad
P_k=\sum_{l\ge k}I_{ijl}{l\choose k}S_{l-k}a^ie^j,
\quad [z^k]Q_1(d_0+z)=P_k/(D b_1^{18}).
\]

S_hの支持はa≤18−h、e≤36−h−2aで、上限19(19−h)。全19本でも3610項。
R cache190と一つのb1累乗を使って各S_hへ直接加算し、全S完成後にR/current workspaceを解放する。
Hasse層はa≤34、e≤69、全次数≤86−k、上限2297項。ここでe≤69−kとは言わない。
input/cache/一accumulator、係数とkey、metadata16000を含む上限37436 slots。

入力L1<2³⁰¹、相関した二項和 \(79^{18-l}179^l\le179^{18}\) により全整数partialは436 bit以下。
hard518を維持する。最初の非零層の逆順確認・再生を含む全source product上限12768337。
CPU完了時間は未測定なので、将来は45秒soft・60秒hard・2 GBと正確prefixを持つ一回だけのworkerを提案する。
全状態prefix上限4642060 bytes、提案枠5 MB、準備入力/コード≤2.8 MB、receipt0.2 MBでRAW8 MB。
実worker準備後の新しいROOT code/source監査が必要である。

## 正の付値だった場合だけの固定除去法

Bの整数係数は−9と−2を含み、9−8=1なのでprimitive。
GaussによりI/Bˢは整数多項式で原D197を保てる。
multivariate lex a>e>dなら先頭単項式は50adで、固定整数LC50による六項の長除算が使える。
未知多項式b1を反転する必要はない。ただしsの値・actual支持から別資源gateを作るまでは実行しない。

既知unit除去は将来の共通κ／二tail補集合入力を小さくするためであり、
補集合の空性、互いに素、実bad Fitting ideal、全FULLU57/239は主張しない。
