# AB の全域正規性：独立概念監査

2026-11-02。[主証明](AB_ORIGIN_NORMALITY_20261102.md)を受理する。
被約 AB は原点を含め正規であり、既知 incidence f の正の高次順像は消える。
これは四次までの寸法一致を全次数へ外挿する議論ではない。
局所 Tor の次数上界から正規化の S 加群生成次数を制限することが要点である。

## 1. 局所模型と独立の有限確認

既存 incidence は Z=LG(2,4) の Zariski chart 上で Y×A³_g×U であり、
Y は sl3 の rank≤1 錐である。N と g は全て次数1。
補平面の変更で g に N の線形式を加えるだけなので、この同型は grading を保つ。
従って相対 Tor を計算する際に、大域的な parabolic splitting は不要である。

独立コードは主 producer を import せず、二元二次 s、二元四次 H、二元二次 g から
F=q s/4+J2 grad H、A=DF と Q を直接再構成する。
六 skew 成分の60係数を11変数で微分し、保存された45座標を選んだ行列の
非零11次 minor を求める。これにより45一次生成元から fiber の11一次座標への全射を
別実装で確認した。局所的にその kernel34 を分裂させ、三つの自由 g 変数を除けば、
Tor は Y の8変数分解と34個の一次 Koszul 因子の積となる。

Y の最小分解の ranks=(1,9,16,9,1)、shifts=(0,2,3,4,6) は主の QQ 分解証明書を入力とする。
その exactness 自体を本レビューで再実装したとはしない。
Koszul 畳み込みは独立に行い、i=0,1,2,3 の内部次数を照合した。
これらの Tor は S_+ に殺され、scheme 論的零 fiber Z 上の局所自由層である。
零 fiber が Z 自身であることも、fiber 代数が11一次座標で生成されることから従う。

## 2. 導来議論と全次数 Bott の監査

proper coherent higher images により K=RΓ(X,O_X) の cohomology は有限 graded S 加群であり、
Z の次元3から cohomological amplitude は0から3にある。
S は45変数多項式環なので K は bounded graded perfect complex で表せる。
同次数の単位成分を消して、全 differential の成分が S_+ に入る最小自由複体 F を取れる。
この場合 F^q⊗k は導来 fiber の q 次 cohomology に等しい。

k が S-perfect なので projection formula を適用でき、

    H^j(Z,Tor_i^S(O_X,k)) ⇒ H^(j−i)(K⊗^L k)

が使える。j≤3 であるため、total degree0/1に寄与する i はそれぞれ3/2以下。
内部次数を保存する spectral sequence の subquotient を取っても上界は増えない。
独立列挙で得た F^q の生成次数上界は

    q=0,1,2,3 に対し 4,3,2,0

である。spectral sequence の退化や Tor 束の大域直和分解は仮定していない。

全 n の Bott 論証も妥当である。core の中心重みは2、vertical は4なので、
vertical 次数gの部分の最高重みは a≥b、a+b=2(n−g)+4g=2n+2g≥0。
ρ=(2,1)を加えた u=a+2,v=b+1 は u−v>0,u+v>0 を満たす。
従って u>|v| であり、v>0なら既に dominant、v<0なら第二座標の符号反転一回、
v=0なら singular である。非零 cohomology の Weyl 長は0または1だけとなる。
有限 filtration に戻しても H^{≥2}(Z,E_n)=0 が保たれる。

この消滅により F の最高の正次数項を取って graded Nakayama を適用し、F^{≥2}=0 を得る。
H¹(K)=coker(F⁰→F¹) は次数3以下で生成される。
既知の H¹(E_n)=0 (n=0,1,2,3)、および E_n=0 (n<0) から H¹(K)=0。
再び minimality により F¹=0 となる。
従って D=coker(F^−1→F⁰) は S 加群として次数4以下で生成される。
これらの次数は既に B に含まれるので、S·B⊂B から D=B が従う。

## 3. 全 n Bott に依存しない第二の有限論証

上の主証明とは別に、低次数 cohomology だけで正次数項を順に消すこともできる。
これは有限一致を無根拠に延長する方法ではなく、同じ導来生成次数上界を使う。

まず F^q=0 (q>3) であり、F³ の生成次数は0だけである。
H³(K)=coker(F²→F³) は次数0で生成されるが、H³(E0)=H³(O_Z)=0 なので零。
最小性により F³=0。

次に H²(K)=coker(F¹→F²) は次数2以下で生成される。
既存の E0,E1,E2 の cohomology 計算から、それらの次数の H² は全て零であり、
負内部次数も存在しない。従って H²(K)=0、さらに F²=0。
最後に F¹ の上界3と既知 H¹(E0..E3)=0 を用いて、前節と同じく F¹=0 を得る。
ここから D の次数4以下生成と全域正規性が従う。
この第二論証は主証明の全 n Bott 部分と独立であり、局所 Tor の低次数と既存低次数の
cohomology を使う有限の確認経路を与える。

## 4. 出典、再現、範囲

導来議論の基礎命題は本文を確認した
[proper coherent higher images](https://stacks.math.columbia.edu/tag/02O3)、
[regular 環上の perfect complexes](https://stacks.math.columbia.edu/tag/066Z)、
[unit cancellation](https://stacks.math.columbia.edu/tag/00MR)、
[perfect 因子の projection formula](https://stacks.math.columbia.edu/tag/01E6) に対応する。
これらを組み合わせた AB への有限生成次数の論証は本研究の推論であり、
出典に AB の結論がそのまま記載されているという意味ではない。

    .venv/bin/python code/ab_origin_normality_review_20261102.py --output reports/AB_ORIGIN_NORMALITY_REVIEW_20261102.json

独立有限検算は0.076秒で PASS。主 proof/source の監査と、別構成の fiber rank、
Tor 畳み込み、内部次数上界を記録した。
Y の QQ 分解 exactness、次数0–4の B_n=D_n、既存低次 cohomology は凍結済みの依存である。
主 producer の再実行や、導来圏の証明をプログラムのみで検証したという主張ではない。
実際の T01 零 fiber の非被約構造、全 S61 の尽くし、成分数の上界は別問題のままである。
