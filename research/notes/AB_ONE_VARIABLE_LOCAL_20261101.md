# AB の一変数・非零トレース軌道に沿う正規性と接錐

2026年11月1日。対象は既約十次元の被約多様体 \(B=\overline{AB}\)。

## 結論

[incidence ノート](AB_INCIDENCE_20261101.md)の特異成分
\(\Sigma_{\mathrm{one}}\) のうち、トレース非零の部分では B は正規である。
この部分は一つの \(\mathrm{Sp}_4\) 軌道で、代表を

\[
 Q_*=x^2\partial_x\wedge\partial_y
\]

と取れる。そこでの局所次元は十、埋め込み次元は十六、Hilbert–Samuel 重複度は七。
接錐は

\[
 \mathbb A^6\times\operatorname{Cone}\bigl(S(2,2,3)\subset\mathbb P^9\bigr)
\]

である。\(S(2,2,3)\) は三次元の有理正規 scroll。
完成局所環そのものが接錐の完成と同型だとは主張しない。
全 B の正規性、特に二つの特異成分の零トレース交叉は未決である。

この証明は [被約ファイバーだけによる誤った近道への反例](REDUCED_FIBERS_NORMALITY_COUNTEREXAMPLE_20261101.md)
を避け、正規化環に必要な一次生成元を実際の Q 座標がすべて与えることを確認する。

## 1. 例外曲線の二 chart

proper incidence \(f:\widetilde B\to B\) の \(Q_*\) 上のファイバーは被約な
\(C\simeq\mathbb P^1\)。トレース非零なので全ファイバーの近傍で \(\widetilde B\) は滑らかである。
\(L_0=\langle\partial_y,\partial_w\rangle\) に対する LG の chart を

\[
 B=\begin{pmatrix}b_0&b_1\\b_1&t\end{pmatrix},\qquad
 S_B=\begin{pmatrix}I&B\\0&I\end{pmatrix}
\]

と書く。行列表示は位置・運動量の順 \((x,z;y,w)\) を用いる。
ファイバーは \(b_0=b_1=0\)、t がそのアフィン座標である。
z,w のシンプレクティック交換を用いたもう一つの chart との関係は

\[
 B'=\begin{pmatrix}b_0-b_1^2/t&b_1/t\\b_1/t&-1/t\end{pmatrix}.
\]

固定 L の \(\mathfrak{sl}_3\) core は \(N_{01}\ne0\) の chart で

\[
N=\begin{pmatrix}1\\a\\b\end{pmatrix}
  \begin{pmatrix}-ad-be&d&e\end{pmatrix}.
\]

代表では \(a=b=e=0,d=-1/24\)。三つの垂直二次係数を \(g_0,g_1,g_2\) とし、
曲線の九法線座標を

\[
(b_0,b_1,a,b,d+1/24,e,g_0,g_1,g_2)
\]

に取る。全テンソルの座標変換を一次まで計算すると、法線遷移は対角成分

\[
 (1,t^{-1},t,t^2,1,t^{-1},t,t^2,t^3)
\]

と、成分 \((2,1)=2\)、\((6,5)=12\) の二つだけを持つ（添字は零始まり）。
すなわち二つの非自明な二次ブロックは

\[
 \begin{pmatrix}t^{-1}&0\\k&t\end{pmatrix},\qquad k=2,12.
\]

これらはいずれも自明束である。実際、二 chart の正則可逆なフレーム

\[
 P_0=\begin{pmatrix}1&t\\0&-k\end{pmatrix},\qquad
 P_\infty=\begin{pmatrix}t^{-1}&1\\k&0\end{pmatrix}
\]

について \(TP_0=P_\infty\)、両行列式は \(-k\ne0\)。
従って正確に

\[
 N_{C/\widetilde B}\simeq
 \mathcal O^6\oplus\mathcal O(-2)^2\oplus\mathcal O(-3).
\]

これは行列式次数だけからの推測ではなく、全遷移と明示的なフレームの計算である。

## 2. 実際の Q 座標が全余法線切断を生成する

Q の六十係数を、上の九法線方向について微分する。
曲線の全 t に沿って得られる \(60\times9\) 行列の係数列を集めると、rank は正確に十六。
証明書には非零な十六次 minor と、全係数行列を保存した。
compatibility の四十五座標と六十係数は線形埋め込みで結ばれているので、
これは S61 の Q 座標についても同じ生成性を意味する。

一方、上の法線束から

\[
h^0(C,N^*)=6+2\cdot3+4=16.
\]

従って、Q の座標差をファイバー ideal に引き戻した写像

\[
 \mathfrak m_{B,Q_*}\longrightarrow
 H^0(C,I_C/I_C^2)=H^0(C,N^*)
\]

は全射である。単一点 t の微分 rank 九だけでは、この全射性は得られない。
曲線全体を残した係数列の rank 十六が必要である。

## 3. 完成環の生成補題

ここで用いる一般的な補題を明記する。
\(f:X\to\operatorname{Spec}A\) を Noetherian 局所整域 A 上の proper birational 射、
X を正規とする。閉ファイバー C は連結・被約で、X の滑らかな領域に入る滑らかな射影多様体とする。
\(I=\mathfrak m_A\mathcal O_X\) は C の ideal、\(E=I/I^2\) とおく。
さらに

1. \(H^0(C,\mathcal O_C)=k\)、\(H^1(C,\operatorname{Sym}^nE)=0\)（すべての \(n\ge0\)）、
2. \(\bigoplus_nH^0(C,\operatorname{Sym}^nE)\) は次数一で生成される、
3. \(\mathfrak m_A\to H^0(C,E)\) が全射、

を仮定する。このとき A は正規で、

\[
 \operatorname{gr}_{\mathfrak m_A}A
 \simeq\bigoplus_{n\ge0}H^0(C,\operatorname{Sym}^nE).
\]

**証明。** Stein 環 \(D=H^0(X,\mathcal O_X)\) は有限 A 代数で、X の正規性と
birational 性から A の正規化である。properness とこの有限代数の同定には
[Stacks Project, Theorem 37.53.4](https://stacks.math.columbia.edu/tag/03H0)を使う。
C が連結なのでその上の有限正規化点は一つ。

C は regular immersion だから \(I^n/I^{n+1}=\operatorname{Sym}^nE\)。
仮定1と、逐次短完全列の cohomology から、各 infinitesimal neighborhood の
大域切断の制限写像は全射になる。
形式関数定理により、それらの逆極限は \(\widehat D\) である。
ここでは任意のファイバーへの \(H^0\) の base change を仮定しておらず、
proper な場合の完成と逆極限の同定を使っている。
[Stacks Project, Section 30.20](https://stacks.math.columbia.edu/tag/02O7)の
Theorem 30.20.5 と Lemma 30.20.3 はこの同定と filtration の cofinality を与える。

この逆極限に infinitesimal neighborhood への核 filtration \(J_n\) を入れると、
その associated graded は \(G=\bigoplus H^0(\operatorname{Sym}^nE)\)。
仮定2、3により、\(\widehat D\) の任意の元について、次数ごとに A の座標の
斉次多項式を選び、残差を \(J_{n+1}\) へ送れる。
この逐次近似は \(\mathfrak m_A\)-進位相で収束するため
\(\widehat A\to\widehat D\) は全射である。
有限 A 加群 \(D/A\) の完成が零になり、完成の忠実平坦性から \(D=A\)。
同じ逐次近似を次数 n から始めれば \(J_n=\widehat{\mathfrak m_A^n}\) となり、
上の associated graded 同型を得る。\(\square\)

仮定3を落とすと、被約・連結なファイバーでもこの結論は偽である。
別ノートの A1 錐の非正規部分環の例では、余法線の三つの一次切断のうち
二つしか元の座標が生成しない。

## 4. 適用と Hilbert 系列

今は

\[
 E=\mathcal O^6\oplus\mathcal O(2)^2\oplus\mathcal O(3)
\]

なので全対称冪の \(H^1\) は零。
各直和成分の二元単項式を積に分解することで、対称冪の大域切断環は次数一で生成される。
前節の rank 十六が補題の仮定3を満たす。従って Q* で正規性と接錐の記述を得る。

非自明な三直和部分の次数 n の切断数は

\[
 \sum_{a+b+c=n}(2a+2b+3c+1)
 =\frac{(n+1)(n+2)(7n+3)}6.
\]

その Hilbert 系列は \((1+6T)/(1-T)^4\)。六つの自明方向を加えた局所接錐は

\[
 \frac{1+6T}{(1-T)^{10}}
\]

を持つ。従って埋め込み次元十六、重複度七。
正の三直和部分の切断環は \(\mathcal O(2)\oplus\mathcal O(2)\oplus\mathcal O(3)\)
の tautological 埋め込み、すなわち \(S(2,2,3)\) の affine cone の環である。

\(\Sigma_{\mathrm{one}}\) の非零トレース点は、\(\ell,w\) の symplectic frame を選ぶことで
すべて Q* に移せる。したがってこれは同開部分全体での結論である。
\(\Sigma_{\mathrm{one}}\cap\Sigma_{\mathrm{plane}}\) と原点については、
incidence の源がファイバー上で特異になるため、同じ regular-embedding 補題をそのまま適用できない。

## 再現

```
.venv/bin/python code/ab_one_variable_local_20261101.py --output /tmp/ab-one-variable-review.json
```

[証明書](../../reports/AB_ONE_VARIABLE_LOCAL_20261101.json)は全 chart 変換、法線遷移、
二つの自明化、全係数微分と十六次 minor を保存する。
完成環の補題は上記の幾何学的証明と、本文を確認した形式関数・Stein 定理を使う。
