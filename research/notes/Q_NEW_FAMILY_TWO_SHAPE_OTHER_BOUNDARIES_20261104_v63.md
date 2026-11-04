# J・Fn・先頭係数境界の有限な必要条件

2026-11-04。r次数19の保存多項式をprimitiveな一般次数と呼ぶ前に、二shape入力の境界を調べた。**J境界に新しい必要因子25r−1を見つけた。** Fn境界のdeterminant・norm条件と先頭係数零の分岐も保存した。一般factorization、新しいfull resultant、主二式のPRS、点・素数probe、source26は行っていない。

## J境界：leading caseを含む厳密な証明

記号は以前の二shape packetと同じで、

\[
 J=(10r-1)s-10r,\quad T=(r-1)R-2r,
\]

\[
 F_n=(r-1)N_q^2+16rT^2(s^2-4rs+4r).
\]

まず

\[
 (10r-1)(s-1)-J=1
\]

なので、J=0とそのs-leading coefficient 10r−1=0は同時に起きない。従ってJ=0の解析でleading caseを黙って捨てる必要がない。

全係数の恒等式として

\[
 \boxed{(10r-1)^2F_n\bigl(s=10r/(10r-1)\bigr)
       =1024r^3T^2}
\]

を検証した。r,Tが既知単元である元のchartではJ=0からFn≠0となる。cleared E₁,E₂はJ=0でZ₀·C₁分子·Fn、Z₀·C₂分子·Fn²へ落ちるので、Z₀も単元ならC₁,C₂の分子が同時に零になる必要がある。

これら一次分子をJ上へ多項式homogenizationで制限し、QQ(R)[r]上のgcdとBézout係数を求めた。gcdはmonic規約で

\[
 \boxed{g=rT(25r-1)/(25(R-2))}.
\]

Bézout係数の分母はrには依存せず、共通R分母R−2を払った**全係数の多項式Bézout恒等式**を保存した。従って一般R（R≠2を表示）で、既知r,T単元を除いたJ境界の必要条件は

\[
 \boxed{25r-1=0}.
\]

対応するJ上のshapeはr=1/25、s=−2/3。これはcleared二式のbase-point候補である。保存Q₁₉をこの因子で割る計算は実行していない。特別R=2の全解析や、一般sourceへのprimitiveモデルの確定を主張しない。

## Fn境界：pseudo-remainderとnorm

Fn=a(r,R)s²+b(r,R)s+c(r,R)と書く。A₁,A₂のprimitive分子（33項、212項）をsでpseudo-divisionし、

\[
 a^{\delta_i}A_i=Q_iF_n+(u_i s+v_i),
 \qquad (\delta_1,\delta_2)=(1,3)
\]

を全係数で検証した。余りは59、352項。aを未知r多項式のまま全域で反転していない。

a≠0の分岐では共通根の必要条件として

\[
 D=u_1v_2-u_2v_1=0,
\]

\[
 K_i=a v_i^2-bu_iv_i+c u_i^2=0\quad(i=1,2)
\]

を保存した。normの根拠となる普遍identityは

\[
 u^2F_n-[a(us+v)+bu-2av](us+v)
 =av^2-buv+cu^2
\]

で、一般scalar変数の係数恒等式として検証した。保存された実係数uᵢ,vᵢをこの式へ代入する。D,K₁,K₂は非零で、(r,R)次数と項数はそれぞれ(29,14)/319、(21,10)/176、(47,22)/828。これらのgcd・零点・Q₁₉との共通因子は調べていない。

leading coefficient零を別に保持する。

* a=0、b≠0なら、b^deg(Aᵢ)Aᵢ(s=−c/b)=0が必要。二つのhomogenized多項式は(r,R)次数(15,6)、(30,12)、79、289項。
* a=b=0なら、c=0がFnに根を持つために必要。c≠0なら根は無い。c=0でのAᵢ条件は未解決として残す。

従ってa=0の可能性を取り落とさず、未知r多項式を全域で反転していない。これらは**必要条件の有限な入力**であり、Fn境界が全て除去されたという結果ではない。

## Hnとs=∞

J·Fnと他の名前付き分母を避けると、二つの正規化写像はHn=0上にもregularに延長する。Hn=0は二shape平面のproper curveなので、その像は独立な二target Z₁,Z₂を支配できない。pole/base locusとの交点は上のJ/Fn問題に含め、無条件には捨てない。新しいHn resultantやsaturationは計算していない。

s=∞については、二式のs次数3,5に対する独立target三係数ずつ、全6個の先頭係数多項式を保存した。一般の独立targetでも残る固定leading rootにはこれら全ての零が必要である。leading content gcdや∞消去は未実行。

## 実装停止と、成功済み証明を再利用した保存

最初の実装は約6.19秒でJ gcd/Bézoutと両pseudo-remainderを保存した後、局所変数が変数tupleをshadowしたため停止した。次のcontinuationは保存済み証明を読み戻し、determinant/normを約0.26秒で保存した後、source recordに`variables`欄が無いというparserの誤りで停止した。

両失敗receiptと元codeは保持した。J gcd・pseudo-division・normを繰り返さず、明示的な変数順で残るleading caseと∞係数だけを約0.77秒で完成させた。**合計の数学処理は約7.22秒で、同じ60秒budget以内**。これは大きい算法の再試行ではなく、保存済み部分からの実装修正である。

最終成功packetは
`reports/Q_NEW_FAMILY_TWO_SHAPE_OTHER_BOUNDARIES_FINISH_20261104_v63/`。
全J証明、norm、leading case、∞係数、入力hashを含む。依存するprefixとcontinuationは同じstemの`OTHER_BOUNDARIES`、`OTHER_BOUNDARIES_CONTINUATION`に、`all_passed=false`の停止記録を残す。対応する三つのPython codeは同名。最終certificate約0.22MB。

Q₁₉のprimitive性、一般r次数、actual target index、NEW方程式は未解決である。次には係数contentを用いる別方式を検討できるが、このpacketでは開始していない。
