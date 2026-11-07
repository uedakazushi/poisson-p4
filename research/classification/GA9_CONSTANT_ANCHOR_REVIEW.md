# Ga9非零定数c全十源の23次元成分性と既存97との差：独立レビュー

2026-11-07。主義務 `E_reduced_nonnormal`。
**c=1の完全十parameter源のPGL飽和閉包は、既存97と相異なる23次元の全Poisson成分である。**
同じ既存b=2 anchorで、真の102次元切片の全方程式から被約支持が原点だけであることを確認した。
全c=1源に対する結論であり、固定したbの部分族だけを成分と呼んでいない。
局所成分一意性、滑らかさ、正確な非被約厚さ、他のGa9 c型の全帰属は主張しない。
文献上の新規性は未確認で、既存台帳97に対する相対的な新成分性である。

独立実装は `code/classification_ga9_constant_anchor_review.py`、
証明書は `reports/CLASSIFICATION_GA9_CONSTANT_ANCHOR_REVIEW_20261107.json`、
有限契約は `jobs/ga9_constant_anchor_review.json`。
生産者・そのhelper・CAS・Singular・SymPy・NumPy・legacy算術をimport・実行しない。
別作者が既に作成した標準Fraction疎多項式、Laurent parser、独自有理行消去だけを再利用した。
最終実行4.304秒、最大RSS23,372KiBで全検査PASS。gateは入力契約の検査と区別する。
TeX・global台帳・Gitを変更していない。

## 1. 完全十源をpotentialから独立に作る

パラメータは `r∈C*`, `(b,α,p,q,A,D,G,H,k)∈A9`、従って全域は既約である。

    B=b−1/4,
    C=−(1−b)(b−1/4)/(2r),
    E=α−[(1−3b)A+q(b−1/4)]/r,
    φ=r y4/4+A y³+B y²w+D y²+E yw+C w²+G y+H w,
    u=x−yz, Vy=−α−by, Vw=p+qy+bw+r y²,
    A2=φ+Vy φy+Vw φw−k。

A2は二次以下となる。補助座標t=zでの `∂t∧V+Jac(u+φ)` を元四変数へ戻し、
六括弧を独立に

    Pxy=φw+yVy, Pxz=u+A2−zVy, Pxw=−φy+yVw+z,
    Pyz=−Vy, Pyw=1, Pzw=Vw

と作った。全十parameterを保った全六式は元の完全源報告と一致し、
全四strict Jacobi、実Pf=`k−u−φ`、divP=`y∂x+∂z`を独立に確認した。
全十源がGa9 c=1の全支持を尽くすこと自体は、既済別作者
`CLASSIFICATION_NIL211_SHORT_TRANSLATION_REVIEW_20261107.json` の全source還元を入力とする。
今回それを選んだ六磁気方向などへ制限せず、全十実微分を保持した。

hで斉次化したquartic Fはx一次で係数が−h³、定数項は−r y4/4を含む。
r非零なのでgcd(h³,定数項)=1。Gaussの補題から全源上でFは絶対既約であり、hはFを割らない。
従って全実五次式は正確に単純hと既約quarticの積、rank四で、H=(h=0)は唯一の一次因子である。

全パラメータ上で実五次式が `(h,y)^3` に入ることも、全monomialについて検査した。
これは一点のPfaffian形の観察ではない。

## 2. 同じ既存anchor、全周囲D98と実F48、真の102cut

既存r=1,b=2,k=1、他零の点を `x'=x−h` でk=0へ送る全150係数等式を独立確認した。
新しい標本点の探索ではない。このk=0点の八倍整数代表は

    Pxy=14hw−2y², Pxz=8hx+8yz+35w², Pxw=8hz−12wy,
    Pyz=16hy, Pyw=8h², Pzw=16hw+8y²。

ここから全350のprojective Jacobiについて、全150方向の微分行列Dを独立に作った。
全有理消去と非零minorからrankD=98。
全25Euler列V、全25GL5実微分、scalarを合わせたG、さらに上の**全十源**実微分を含むFを別に作ると

    rankV=24, rankG=47, rankF=48, DF=0。

十源の係数微分もrank10を保つ。群方向だけのrank47を全sourceのrank48と混同しない。

保存された四法線N_iを全150係数へ復号し、全Dで消えることとrank[F,N]=52を確認した。
98本のliteral標準座標Uと四法線を合わせたLについて、独立有理minorから

    rankL=102, rank[F,L]=150

を得る。従ってP+Lは48個の独立な実線形cutで定まる真の切片である。
選択した元Jacobi98行とU98列のminorもQQの非零定数として独立に再計算した。
法線だけの支持を周囲へ外挿する証明ではない。

## 3. 全orders有限性は実行前のC* gradingから従う

ambient座標weightは `(h,x,y,z,w)=(0,4,1,3,2)`、点のtensor characterは−3である。
全150列、全350行のweightを各monomialから別に作り、Dの全非零成分のweight一致を確認した。
Fの各実列とLの各列もweight homogeneousである。
点を固定する対角PGLとtensor scaleの組合せを逆向きにすれば、四法線の正weightは

    (1,2,2,3)。

全ambient implicit座標の相対weightは−4から10、元Jacobi全行は−4から16。
これらは全表現から得た先験的な上限であり、jetの次数を任意に切ったものではない。

非零98次minorで、真のL内のU98座標を四法線の形式函数として一意に解ける。
DN=0なので補正は通常次数二以上から始まる。
選択した行・列はweight stableなので、一意性からその形式IFTはC*同変である。
全法線weightが正のため、weight dのimplicit座標はweight dの有限多項式となる。
負weightとweight零の補正は零、正weightは10以下。
従って全IFTはweight10以下で全ordersが完成し、全残式もweight16以下で尽くされる。
QQ(b)の一様なformalモデルや別のparameter開を要しない。

## 4. 全78implicit多項式を元350式へ直接戻す

normal t0,…,t3について、元150係数のtensorを

    P(t)=Panchor+Σt_i N_i+Σ_{j∈U} f_j(t)e_j

として別算術で構成した。
全98implicit座標のうち78個のf_jが非零であり、全monomialが指定weightを持ち、通常次数二以上である。
この完全tensorから全Schouten/Euler外積を直接展開して、元350式の全係数を取り出した。

* 選択された元98式は多項式として恒等零。
* 残る非零179式は保存された全179式と全係数で一致。
* 他の全零行、負weight行、weight零行も省かず零と確認。
* 全残monomialのweightは元の各行のweightに一致し、最大16。

全補正は法線項を含め808 monomial、全残式は5,412 monomialである。
非零minorと形式解の一意性から、この検証済み多項式解は真の完全IFTそのものである。
全350式を課した真の切片の完成局所環は、これら179式の四変数形式商に一致する。
高次残差や未表示の相補行を捨てていない。

## 5. 原179式からの四つの明示liftで根基原点を証明

原179式をR_jとする。証明書の四つの多項式G_iと全係数L_ijについて

    G_i=Σ_j L_ij R_j

をFraction疎多項式として直接検査した。
必要な非零lift係数の個数は2,11,11,11であり、各係数は通常の多項式、法線の逆元はない。
これらの多項式は

    G1=t1⁴, G2=t2⁴, G3=t3³,
    G0 mod(t1,t2,t3)=C t0⁵,
    C=514114689789605014659184048028741872800≠0

を満たす。各leading monomialも独立に取り直し、t0⁵,t1⁴,t2⁴,t3³となる。
三つのliteral pure powerから全共通零点でt1=t2=t3=0、続くG0からt0=0。
元の全式は原点で零だから、根基は正確に `(t0,t1,t2,t3)` である。
原nativeの26本の完全Groebner基底、最小基底、正確なnil指数を仮定せずに支持を閉じた。

## 6. 一点の真のcutから全c=1源の成分性へ

真の48線形cutの完成局所環は零次元なので、Krullの不等式により全raw周囲局所次元は48以下。
実際の全十源・群・Euler・scalar写像の微分rank48が像次元の下界48を与える。
全十源のdomainはC*×A9、PGL/GL群・scalar・Euler追加も既約なので、像閉包Yは既約である。
anchorを含むYの次元は48以下で、上の下界と一致する。
Yを含む周囲成分もanchorで次元48以下だから、Y自体が全空間成分となる。
Euler24とprojective scalar一を除いた成分次元は23。

結論はb=2の部分族だけでなく、全b・全r非零・他の全九係数と、その像閉包の全境界を含む。
scalar/Ga9の他のc型を全て被覆したという意味ではない。
接次元は27であり、滑らかさ・局所唯一成分・完成局所環の積分解は主張しない。

## 7. 唯一Hとcofactor flag：旧四H+quartic成分を直接識別

Hが実五次式の唯一の一次因子なら、h-Casimir homogeneous liftは一意である。
差E∧Lもh-Casimirなら hLj−xjLh=0、linear LについてL=cE、従って差は零となる。
よってその斉次発散Nは共役とtensorの非零scaleを除き内在的に決まり、Hとのdual-image flagも不変である。

全新源の実発散を全括弧から再構成すると、homogeneousには

    Nnew=y∂x+h∂z, rankNnew=2, Nnew²=0,
    Nnew* (dz)=dh、従ってh∈imNnew*。

旧ESSENTIAL5、M、A0、DPMは各印刷された六実括弧を**全て別に再入力し、発散を直接計算した**。

| 旧23成分 | 直接得た実cofactor | rank | Hとのflag |
|---|---|---:|---|
| H_QUARTIC_ESSENTIAL5 | z∂y+x∂w |2|h∉imN*|
| MINIMAL_M | 6x∂y |1|h∉imN*|
| A0 | z∂y+x∂w |2|h∉imN*|
| DPM | z∂y+x∂w |2|h∉imN*|

DPMは共通normal formへu=1/6,v=−1/36,c=ρ/18,d=1/2+ρ/18,a=−ρ/6,b=−ρ/216、
p=x²+2xz+3z²を代入し、ρ=√3の係数を保ったまま微分した。
ρは発散から恒等的に相殺し、浮動小数やCAS、ρ²の数値代入は使わない。
旧全5×5行列のh-columnが零であること、rank、dual imageにhを加えたrank増加も独立計算した。
従ってA0/DPMのJordan型を名前から仮定していない。

旧四成分のcanonical PGL軌道は既済証明により23次元で稠密開である。
新しい全c=1飽和像はconstructibleかつその新23成分に稠密で、全像で上の反対flagを持つ。
もし成分が等しければ両稠密constructible像は非空開で交わり、flagが矛盾する。
このflagを全退化上の閉条件と仮定する議論は使わない。

## 8. 全97との相違と凍結projectionの監査

`CLASSIFICATION_GA9_CONSTANT_PRIOR97_20261107.json` の全97行を走査し、各IDへ個別の理由を保存した。
rank二はrank、23次元でない成分は全成分次元で識別する。
同じ23次元のMPS/SNC成分とは、全新sourceに成立する「ある平面に沿って五次式が三重以上」という閉条件で区別する。
この条件はGr(3,5)上の H0(I_plane³(5)) の射影束からのproper像で閉であり、退化式零もaffine coneに含める。
SNC因子は一般codimension二点で高々二枝しか交わらないため、この条件を満たさない。

他の同次元非MPS成分は、第7節の旧四H+quartic以外、一般因子次数・重複度が1+4と異なる。
そのproper因子積写像の閉像に今回の既約quarticを持つ点は入らない。
全97に未処理の同次元H+quartic例外がないことも明示検査した。
これで新成分は全既存97と相異なる。

凍結projectionのmetadata headが0f0cb285であるのは、独立受理済みのローカル97台帳を
97 commit前にfreezeしたためである。immutable inputのmetadataは変更していない。
レビューはmutableな現行台帳を依存にせず、Git blob
`741c35604ee905eba799fd2dffa6c001a5a20b10` を直接読み、
SHA256 `f77a7c6e455f079f7c5de18628feadd4cfc7288f1cee1d67b5e825b7335c9f94` と一致、
全97 IDとprojection全fieldがそのblobと一致することを確認した。
このblobは公開済み97 commit `d27d787f77b6d7e91b5b97dac2fc8fa343d45107` の台帳である。
後続の台帳追加によってこのレビューが再実行不能になる依存を作っていない。

以上により全c=1源は一つの新23成分へ被覆され、既存97に一成分を追加できる。
実際の成分数・台帳の昇格はroot integratorの担当である。
