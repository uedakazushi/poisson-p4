# Compact normal-direction graphの全42chart：独立監査

`S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63` を、指定beta-only開上の **compact incidenceのnormal方向P41への閉埋込み** と、その全42標準chartの有限多項式recipeとして受理する。独立検算は約1.5秒で完了した。アフィン残余の空性、QQ上のB_i冪零性、一般五点尽くしは含めない。

## 1. Compact対象の正確な定義と全134式

元のdegree2以下の134式を、先行有限5生成packetの記法で五式F_lと129式G=A(z,n)nへ移す。cross inverseによる可逆な定数行操作を行うと、Gは

\[
 z_kn=P_k(n)=T_kn+Q_k(n)\quad(k=0,1,2)
\]

の126式と、三式H=H1+H2となる。五式は同じcross rewriteを引いて `f_l(z)+K_l(n)` にする。これらはq次数2以下で、行変換はqに依存しないため、wでdegree2に斉次化しても全134式との等式が保たれる。affine idealのprojective closureへのw-saturationを済ませたと仮定していない。

対象YはP45_[w:z:n]×P41_[v]内の **divided normal incidence**

\[
 F_0^h=\cdots=F_4^h=0,\quad A^h(w,z,n)v=0,\quad
 n_iv_j-n_jv_i=0
\]

である。単に未除算の134斉次式とn∧vだけを課した対象ではない。後者はn=0の既知五点上にP41を残してしまう。

方向chart v_i=1でn=rvと置く。`G(w,z,rv)=r A^h(w,z,rv)v` は全多項式恒等式なので、**ambient polynomial ring内で**共通rを消せる。これによりincidenceのcross式は、r=0も含めて

\[
 z_kv=wT_kv+rQ_k(v)
\]

となる。quotient ringで未証明のr正則性を先に仮定した操作ではない。

## 2. r=0とnormal-zero無限遠を別々に排除

Yのr=0上では `z_kv=wT_kv` である。w≠0ならvはT0,T1の共通固有vectorとなるので、可逆なC0=[T0,T1]がvを消し、projective vと矛盾する。w=0ならz_kv=0からz=0となり、[w:z:n]がprojective点でなくなる。従ってYはr≠0に入り、schemeとしてこのopenを使ってよい。

一方、**元のraw134射影scheme** のn=0部分は既知planeとの交叉である。先行の五quad最高部rank5・その一次倍rank10により、zの全三次単項式がその最高二次idealへ入る。従ってw=0,n=0にはprojective点が無い。w≠0,n=0のfiberは先行の実際のE5である。

この第二の事実により、raw無限遠の全点はn≠0であり、[v]=[n]を付けることでYのw=0部分に移る。既知五点の除去と、無限遠を失わないことを混同しない。

## 3. 全42chartに共通するQQ恒等式

lambda_k=(T_kv)_i、x=v_iとし、U_k=lambda_k v−xT_kvとする。独立codeでは、任意の行列文字A,Bと可換scalar x,l0,l1の自由wordを展開し、整数係数で

\[
 (xB-l_1I)(l_0I-xA)-(xA-l_0I)(l_1I-xB)
       =x^2(AB-BA)
\]

を検査した。A=T0、B=T1、l_k=lambda_kと代入し、vを掛けると

\[
 (v_iT_1-\lambda_1I)U_0-(v_iT_0-\lambda_0I)U_1
     =v_i^2 C_0v
\]

である。row_i(C0^-1)を掛ければv_i³となる。これは **全iについて標数0で成立する一つの形式的恒等式** であり、主packetの三chartのmodp検査から全42へ推測したものではない。

C0のQQ上の非零性・p単元性は先行のconstant commutator証明と独立reviewを保持する。逆行列はadjugate/determinantでbetaだけを分母とする有限recipeである。今回も実際のmodp C0と保存逆行列の両側積を全42×42 entryで確認し、exact QQ inverseの還元との一致を保証した。

v_i=1とするとu0U0+u1U1=1である。V_k=Q_ki(v)v−Q_k(v)、B_i=u0V0+u1V1を定めれば、全mixed compatibility `sU+V=0` は

\[
 s=-B_i,\quad \zeta_k=-B_i\lambda_k+Q_{ki},\quad
 V_{kj}-B_iU_{kj}=0\ (j\ne i)
\]

と等価となる。実際、compatibilityのu結合がs+B_iを与え、各式からU(s+B_i)を引けば残る式になるので、これはidealの等号であり根基だけの議論ではない。

三つ目のk、三本のH、五本のFを全て保持する。したがってchartのidealは123本のcompatibility、三本H2−B_iH1、五本の

\[
 f_l^h(\zeta,s)+sK_{l,1}(v)+K_{l,2}(v)
\]

をs,ζの上の式で置換したものとなる。次数上界はそれぞれ6,5,10、B_iは4以下、ζは5以下である。一般のplane原点でf_lの定数項が非零なら、f_l^hにそのs²項も含める。

## 4. 係数と元sourceへの独立binding

主producerや`load_fiber`をimportせず、保存されたexact QQ出力の全119070 P係数をmodpへ還元し、compact modelと全て照合した。先行の独立有限5生成reviewは、この同じ出力を整数分母消去で検査し、全134式・129 factorizationに結び付けている。

今回さらに、公開cross入力から三本のHを別に再構成した。五本Fについては、元61変数sourceの単項式にT01点とaffine45座標変換を直接代入し、plane部・normal部・mixed部を分類して全Kを作った。Pだけでなく、全H、全K、五plane式が保存modular modelと一致した。従って第三mixed blockや残り三Hを落とした小模型を検査しているわけではない。

全42chartの環同型は同じintrinsic Yの上で構成されるので、重なりでも一致する。各chartでs,ζを多項式として解き、残る方向idealとの互いに逆な環準同型がある。よってY→P41は **schemeとしてclosed immersion** であり、未知のnilpotent構造も保持する。

## 5. 無限遠とアフィン部分の範囲

元134式のQQ最高q部分を直接調べ、qのunit axes8,16,24,32が凍結中心で全134式を満たすことを確認した。対応するu11,u22,u33,u45の平方項は、順に一つの元式の `2u11²u36`、`2u22²u27`、`2u33²u13`、`2u45²u5` だけであり、これらbeta係数が実際の中心で零となることを確認した。対応する全P二次係数も **QQで厳密に零** であり、これら四点はs=ζ=0のgraph点となる。一般betaで同じ四軸が残るとは主張しない。

別の一般性の根拠として、既存BOUNDARY_DIMENSIONが証明する **full135 raw無限遠の正確な次元3** を依存として保持する。この3次元族は、n=0無限遠が無いことによりcompact134-incidence Yの中に入る。これで一般betaのYが非空であることが分かる。今回、Y134全体の次元上界3やpure dimension3を新たに受理したわけではない。

実アフィン残余はs≠0、すなわちB_i≠0であり、

\[
 n=-v/B_i,\qquad z=\lambda-Q_i/B_i
\]

で回復する。B_iを反転するのはこの **アフィンopenを指定する最後の段階** だけで、compact graphの構成では反転しない。Yのbeta上のpropernessをこのopenへ移してはならない。

modpでのB_i冪零性だけからQQ冪零性は従わない。properかつflatな点の族でも、固定hyperplane内の特殊点が一般fiberでは外へ動き得るため、flatnessだけでも十分ではない。QQ上の実cofactor等式、またはこの関数のsaturationを制御する別証明が必要である。compact Y全体の空性を目標にせず、s≠0部分を保った問いとして扱う必要がある。

独立codeは `code/s61_affine_normal_direction_graph_review_20261103_v63.py`、reportと小さい自由word・source-axis certificateは `reports/S61_AFFINE_NORMAL_DIRECTION_GRAPH_20261103_v63_REVIEW/`。全指定input/evidence hashを照合した。42変数GBや一般QQ chart idealの巨大展開は行っていない。
