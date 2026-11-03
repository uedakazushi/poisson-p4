# Independent audit: B18/W uniform Darboux map, orbit and leaf closures

2026-11-03。指定された算術・幾何の主張を**受理**する。producerは再実行せず、元CA/CBのraw model、coefficient convolution、三角vector fieldの積分、projective degreeを別に構築して検査した。局所完成型そのものや全B18 strataの分類はこの監査の対象外である。

## 元CA/CB曲線とcanonical150

CAのraw modelへf=−h²、g=W、Hamiltonian=x²z−hxw+z³/3を入れる。CBへ同じf,Hとg=−Wf=W h²を入れるのでbinary pencil条件も満たす。これらをcanonical divergence projection後にh,xを交換すると、両方から

    P_W=∂y∧(x²∂h+hx∂z+(h²+z²)∂w)+W x²∂z∧∂w

を得る。この実sourceを、保存された全W sourceおよび独立受理済みodd-quadratic sourceのU=V=0に結びつけた。canonical150係数を別のmonomial order engineで作り、指数による微分・畳込みから全350 strict Jacobi係数と25 divergence係数が零となることを検査した。W標本による外挿ではなくQQ[W]の恒等式である。

Euler/Pfaffianからdegeneracy quinticはW*x^5となる。従ってW≠0のrank4 locusは正確にx≠0である。

## 積分による独立Darboux導出

x=1でW=0のbivectorは∂y∧D、D=∂h+h∂z+(h²+z²)∂wである。D(h)=1からZ=z−h²/2を得る。z=Z+h²/2を使い

    D(w)=Z²+(1+Z)h²+h⁴/4

をhで積分するとF=hZ²+(1+Z)h³/3+h⁵/20となる。T=w−F(h,Z)はD-invariantである。これが作者のforward mapを独立に与え、逆はz=Z+H²/2、w=T+F(H,Z)である。

全六bracketと両合成をQQ[W]で検査した。変換はWに依存せず、W=0も含めてpolynomial automorphism、Jacobian1である。W≠0ではconstant symplectic A4、W=0ではPoisson centerが正確にQQ[Z,T]となる。後者の全affine leavesはZ=c,T=dのsmooth connected fibersで、H,Yを自由座標とするA2である。

## Cremonaとprojective orbit

独立に得たaffine mapsをdegree5へ斉次化し、両合成がx^24倍の恒等写像であることを全五座標で確認した。x-coordinateがx^5であり、最後のcoordinateが非零h^5項を持つため、両方の全coordinate formsのgcdは1である。従ってbidegreeは正確に(5,5)である。base locusはh=x=0のP2であり、射影空間全体のregular automorphismとは主張しない。

全inverse-coordinate substitutionを含むtensor pushforwardを別に行い、g_r=diag(1,1,r,1,1)についてg_(r*)P_W=rP_(W/r)を確認した。r≠0であるから、全W≠0はprojectivized bivector spaceの単一PGL orbitとなる。W=0への非可逆な代入をorbit証明に使っていない。

## 全Casimir値でのleaf cone degree

逆Darboux mapへZ=c,T=dを入れ、hをP1 parameter A/Sとした四sectionをdegree5へ斉次化する。S≠0ならx=S^5が非零、S=0ならw=A^5/20が非零なので、**全有限c,d**でbasepoint-freeである。h/x=A/Sは有理逆parameterなのでmap degreeは1、pullback O(1)=O_P1(5)からimage curveのdegreeは5となる。特殊c,dでdegreeが落ちる余地はない。

yが独立な自由座標なので、leaf closureはこのirreducible rational quinticとvertex[0:0:1:0:0]のjoin、すなわちrational coneである。vertexを避けるhyperplane y=0のsectionがそのquinticだからcone surfaceもdegree5である。この議論はnaive homogenized equationsをcomplete intersectionと誤認していない。

## 保存と範囲

`code/ambient_B18_W_uniform_darboux_review_20261103_v63.py`、`reports/AMBIENT_B18_W_UNIFORM_DARBOUX_20261103_v63_REVIEW/`。元producer/dataのhashを保存し、作者のfilesは変更していない。全Wのlocal completionやglobal component classificationの新しい主張は追加しない。
