# v=0境界：全Ga stabilizerと旧209 sourceの独立監査

2026-11-04。凍結 `Q_NEW_FAMILY_V_ZERO_BOUNDARY_20261104_v63` を対象とする。
領域は全u≠0,v=0。producerを呼ばず、元全Pと旧sourceのencoded方程式から検査する。

## 全stabilizerとparameter同値

元Pからdivergenceとhomogeneous Nを再構築した。N²=0、ImNのminorは4u²/3。
全quadricのdetは32768u⁴/531441で、Q|Hはrank3・vertexはz方向。
従ってH、vertex、ImNを使う全affine旗の議論はこの境界でも有効である。
compatible constant pencilが一意だと仮定していない。

source u、target Uとscalar κを最初から分ける。四つの低bracket一次係数は
eta=κfに対しeta*b=64c,eta*a=d,eta*d=a,64eta*c=bを与え、divergenceはeta*u=U。
affine blockが可逆なのでeta²=1であり、U=±uが必要。元full tauの全六係数も再照合したため
U=−uも実現される。これはv≠0のfirstjet逆写像を境界へ延長する論証ではない。

固定uではeta=1。modular covarianceのz,w係数からb=c=0、z,x係数からj=a³。
従ってf=a²、κ=a⁻²。y'のe,gは未知のまま保持する。
full Pzwの一次w係数は−2g/3なのでg=0。
full Pyzのw²係数a³−a²−64agからa=1、xw係数−a(128e+ug)からe=0。
constant条件はtx=tw=i=0を強制する。u²=256等を除く必要はない。

残るH,K,Zの三constant系を別のQQ(u,T,H,K,Z)係数環から抽出した。
detは−256u/243で、唯一の解はH=3T/4,K=0,Z=T/(64u)。
元Pをhomogenizeし、全150 raw係数の一致を独立検査した。
全行列はI+3TN/(4u)であり、N²=0なのでGaそのもの。
従って**各u≠0の全projective stabilizerは正確にGa**。
v≠0での有限Eがこの境界にも残るとは仮定しない。

parameter間同定は±だけで、各projective orbitは24−1=23次元。
PGL5×Gm_uからのsourceの一般fiberはGaの有限個のcosetなので、飽和の閉包は24次元。
元三parameter族のv→0特殊化ゆえ新25成分の境界へ入るが、新しいPoisson成分とは呼ばない。

## 旧99 source／209条件のactual binding

旧exact99 graphを(n0,n1,n2)=(u,1,0)で特殊化し、独立なQQ[u]疎演算で元のencoded
209方程式全てを評価して零と確認した。保存された209 pullback文字列の零だけを
読み直したのではない。最初の90係数を旧positive Pへ戻した全六式も元Pと一致する。
さらにactual Pから140 Jacobi係数、24 compact-center係数、25 N²係数、20 constant-pencil
係数を別に再構築し、同じ209条件を確認した。

九個のauxiliary座標は全て零。markはcenter=z方向の無限遠点、pencil=∂y∧∂zのrank2。
従って**閉じた209-incidenceへの所属**は証明される。
有限center・非退化pencilの開部分の閉包へ入ることを意味せず、N²=0だけからの
MPS N成分帰属でもない。

## 既知成分の排除とその根拠

受理済みglobal proper CA/CB sourceは全成分を覆う。
その非零quinticはそれぞれg h³ f(h,x)、a b h F(h,x)²であり、幾何的には五つの
linear factorの積。これはproper Chow像として閉条件で、zero divisorも別に含める。
今回の非零3H+smooth irreducible Qはその型でないのでCA、CBとN⊂CAから排除される。

既知Hquartic成分はH+irreducible quarticを持つ23次元の単一orbitの閉包であることが
旧独立reviewと成分定理で受理済み。今回も各orbitは23次元だがdivisor型が異なるため
そのdense orbitではない。境界は次元≤22なので境界にも入れない。
ここでもclosed209-incidenceへの所属とHquartic成分への所属を混同しない。

その他の既知成分、u=0/p=0境界、新成分の全境界の網羅は主張しない。
新しい点・prime・GB・solver探索は使わない。

再現：`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_v_zero_boundary_review_20261104_v63.py`
（60秒/2GB上限）。
