# 四core二次近似の節点四次曲線：独立監査

既に独立受理された二つの実二次初期式だけを読み、別のQQ行列・多項式演算で検査する。source jet、元の幾何builder、新しい点や素数は使わない。補助conic上の既存の有理方向は、betaの新しい試験点ではない。

二次形式の行列A,Bはrank4,3で、det(A+tB)=(25t²+5t+3)/691200。有限二根は単純なので、それらのpencilのvertexはbase curveの特異点にならない。Bのvertex v=(7/5,2/5,7/24,1)は両式上にあり、滑らかなA上の接平面をvで割った二次形式は非退化である。従ってここが唯一の特異点で、ordinary nodeである。

滑らかなAのquadric上で曲線は(2,2)、算術種数1。可約なら各真の成分はbidegreeの少なくとも一方が2未満なので正規化の種数は0である。唯一のnodeのdelta=1と連結性による種数公式は成分数1を強制する。従って幾何学的に既約である。

さらにvからb=0平面への射影で、Q(w)=w転置Bw=0という滑らかなconicを得る。T=v転置Aw、F=w転置Awとすると、Phi(w)=2T(w)w−F(w)v。Q=T=Fに共通射影零点が無いことをbinary resultantと無限遠点で確認し、二つのquadricの恒等式も全係数で検査する。Phiは全conicで定義され、v以外では射影の逆になるので有限双有理、すなわち正規化である。既存の方向(8,1,7/6,0)がconic上のQQ点であるため、正規化はQQ上のP1となる。Tとの交点は相異なる二点であり、nodeの二枝も保たれる。

この二次idealはprimeなregular sequenceである。従って実最初二つのcore函数についても標準filtered regular-sequence論からassociated gradedがこの二次商となり、二次元・重複度4の整なCM完全交叉となる。高次の項によって実surfaceの正規性は変わり得るので、二次近似から非正規とは推論しない。

これは七函数全部のinitial idealでも、full core tangent coneでもない。残る五函数の高次・socle・E4/D0の深さは未解決である。

再現は `code/ambient_B18_core_nodal_quartic_review_20261104_v63.py --output <新しい出力先>`。全入力hashと行列・射影式・resultantを同名reportに保存する。
