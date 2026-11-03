# 四次quadratic-product商の反射性

2026-11-03。これは全四次Jacobi moduleの全射を証明するものではなく、
そのpure quadratic-productsを割った周囲加群の性質を決定する補題である。
係数環は実三パラメータの正則integral環
`R=QQ[a,e,d,1/(a*d*F*B*C*J)]`とする。必要な既知定数の反転は省略した。

六quadのquadratic cofactor写像を `M:R168→R210`、商をQとする。
全22 syzygy定理から、kernelは点別独立な22本の張るsubbundleであり、
像Eはlocally free rank146である。従って

`0→E→R210→Q→0`

が完全で、Qのprojective dimensionは高々1、generic rankは64である。
ここでEを周囲R210のsubbundleとは呼ばない。特殊fiberでの単射性は
未証明であり、それを使う議論ではない。

既証の全aPJ有限欠損定理をFULLUへ制限すると、Qは有限個の閉点を除き
locally free rank64である。各残る閉点の正則局所環R_pは次元3。
上の完全列とdepth lemmaから `depth Q_p≥2`。Qのtorsionが非零なら
その支持は上の有限集合に入り、そのassociated primeは閉点である。
部分加群のassociated primeはQのassociated primeでもあるので
`depth Q_p=0`となり矛盾する。従ってQはtorsion-freeである。

高さ1、2の点ではQは自由、高さ3ではdepth≥2だからQはS2を満たす。
正則環はnormalなので、有限加群のtorsion-freeかつS2という判定により
**Qはrank64のreflexive加群**となる。また上の長さ1の局所自由resolutionから
`Ext^i_R(Q,R)=0 (i≥2)`、`Ext^1_R(Q,R)`は有限閉点にのみ支持される。

三次元正則環上のreflexive加群には孤立した非自由点があり得る。
従って本補題は有限例外でのrank146、Qの自由性・flatness、122×64全四次
残余の全射を意味しない。残余四次方程式はこの反射加群Qへの実際の写像と
して扱えるが、そのcokernelを零とするには追加の証明が必要である。

依存資料は全22 syzygyの独立レビューと、全quad78有限欠損からの
light146有限欠損の独立レビューである。新しい数値計算・標本・GBは使わない。
