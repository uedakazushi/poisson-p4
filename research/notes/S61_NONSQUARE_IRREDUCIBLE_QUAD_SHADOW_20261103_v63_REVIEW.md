# 既約二次式によるproduct rank37下界：独立監査

`S61_BINARY_NONSQUARE_IRREDUCIBLE_QUAD_SHADOW_20261103_v62` を、標数0の全幾何点、開

\[
 U=\{adFBCJ\ne0\},\qquad J=5ad-6e-4
\]

における **元の42本のquadratic-linear productの階数≥37** として受理する。独立な小計算は0.085秒で終了した。全product rank41、full cubic57、quartic、ambient除去、J=0への適用は含まない。

## 1. 元の六light二次式からのrank3証明

公開mod-ell座標記録のq0を全14単項式から独立に再構成した。七つ目のheavy指数が全て0で、各light次数が2、係数に元のdが含まれないことを確認した。clear multiplierも1であり、再構成した多項式は保存された **unscaled q0そのもの** と一致した。従ってP=adへの変更でa-clear shiftは不要である。

主側のsparse permutation determinantを呼ばず、SymPyで実際のq0のw0,w1,w2 Hessianを取った。その行列式は

\[
 \det\operatorname{Hess}_{012}(q_0)=125a^2/32.
\]

対称係数行列はHessianの1/2なので、その3×3 minorは正確に `125a²/256` である。a≠0の全点でrank≥3となる。標数が2でないとき二つの一次式の積の対称行列rankは高々2なので、q0は代数閉包上でも既約である。

四つのmod-ell二次式の全点での独立性には、先行の `M L=−a²P²J I4` のexact polynomial certificateと独立reviewを依存として保持した。今回、その非空一般点だけを全点へ読み替えていない。

## 2. 3平面の存在を与えるincidence

6次元空間Vと4次元quadric空間Wを各幾何点で固定する。P(W)は3次元であり、rank≤2 locusはq0を含まないproper closed subsetなので次元≤2である。

rank rの非零quadric qについて、radical次元m=6−r、L∩rad(q)の次元k、非退化商への像次元s=3−kとする。対応するisotropic3-planeのstratumは

\[
 \operatorname{Gr}(k,m),\quad \operatorname{OGr}(s,r),\quad
 \operatorname{Hom}(S,\operatorname{rad}(q)/(L\cap\operatorname{rad}(q)))
\]

を順に選ぶため、次元は

\[
 k(m-k)+s(r-s)-s(s+1)/2+s(m-k)
\]

である。非退化商のisotropic Grassmannianの式は、Gr(s,r)に対するs(s+1)/2個の対称双線形条件から得られる。char0の非退化形式ではその微分は全射であるため、この次元評価が使える。最大isotropicの場合に二成分となっても各成分の次元は同じで、上界に影響しない。

独立codeではsを走らせて全実行可能stratumを列挙し、rank1,…,6の最大次元

\[
 (6,6,4,4,3,3)
\]

を確認した。rank0はP(W)の非零quadricには無い。

従って `([q],L)` with q|L=0 のincidenceは、rank≤2側で次元≤2+6=8、rank≥3側で≤3+4=7であり、全体でも8以下である。Gr(3,6)は9次元なので、投影像は全体にならない。よってW→Sym²L*が単射となるLが存在する。これは **各core点ごとの存在** であり、一つのLが全coreで一様に使えるという主張ではない。

## 3. 8+12と、元の17 kernel方向

三変数の六個の二次単項式から四つを選ぶ全15組を、主側の指数vector列挙と異なるmultiset表現で検査した。三次shadowの大きさの分布は、8が9組、9が3組、10が3組である。四つの独立な二次式の基底を行簡約してleading monomialを異ならせれば、異なるshadow monomialに対応する積は一次独立なので、制限したlinear productの階数は8以上となる。

Lのidealの一次基底をell1,ell2,ell3とする。12本ellj qiの定数関係をI_L²で割ると、自由なconormal基底ell1,ell2,ell3の各成分に、独立なqi|Lの関係が生じる。その係数は全て0である。従ってrestriction kernelには12次元があり、imageの8次元と合わせて六light product階数は20以上である。

先行pair coverは、元のell,Z,Wが全開で独立であり、q0,q1,q3,q5が(Z,W)を法として独立であることを証明している。元の七light空間で

\[
 \ell(Z,W)S_1
\]

は13次元である。二つの独立一次式の14本の一次倍に唯一の定数関係ZW−WZがあるためで、ellを掛ける操作は多項式環で単射である。さらにell qiの四方向はこの13次元を法として独立である。ellを消去した関係が、(Z,W)を法とするqiの関係になるからである。

以上の17方向は全てell制限のkernelにあり、imageには上の20方向がある。元六二次式のunit row changeも含めて

\[
 \operatorname{rank}(S_1\langle q_0,\ldots,q_5\rangle)
 \ge17+20=37
\]

を全点で得る。これは抽象的なrow moduleの自由rankではなく、実際の多項式係数行列のpointwise rank下界である。

## 4. 第33pivot後の非構成的帰結

独立受理済みの33回のproduct pivot後、残る8product行の階数は全点で4以上となる。左核のprojective fiberはP7内で次元≤3なので、3次元coreとのzero-combination incidenceの次元は≤6である。P7への **像の閉包** も次元≤6であり、proper closed subsetとなる。この投影のpropernessは必要としていない。

その補集合にはQQ点lambdaが存在する。するとlambdaによる一定row combinationのentryはU上のどの幾何点でも同時に零とならない。localized affine coordinate ringへのNullstellensatzから、そのentry idealは単位idealである。少なくとも一つのlambda成分は非零QQ定数なので、このrow combinationは可逆な定数行変換で導入できる。

これは次のunimodular scalar rowが **存在する** という定理であり、lambda、Bézout vector、あるいは次のnative出力を構成したものではない。**実行済みpivot数は33のまま** とする。全41product方向を尽くしたとはしない。

独立codeは `code/s61_nonsquare_irreducible_quad_shadow_review_20261103_v63.py`、小certificateとreportは `reports/S61_NONSQUARE_IRREDUCIBLE_QUAD_SHADOW_20261103_v63_REVIEW/`。主packet、rank4右逆の既存review、pair cover、第33pivotの全指定input hashを照合した。大行列solverは起動していない。
