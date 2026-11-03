# CA self-pair上の全CB fiberとsmoothness：独立監査

2026-11-03。指定liftのrank29だけでtargetの滑らかさを推論せず、全proper sourceを先に尽くす。

実P(a)=∂y∧D(a)から75×5の定covector contractionを再構成し、正確なrank5 minor −864a⁴/625を確認した。
著者の展開式とは別に、各4成分のdet[E,v,∂y,D]から175×5の全projective constant-factor行列を構成した。
全係数が著者証拠と一致し、rank4 minorは3a、核は<∂y>である。
従って以下はいずれも**全a≠0**で成立する。

既存の独立監査済み全CB proper模型では、rank2のpreimageは二つのratio端点a_CB=0、b_CB=0で尽くされ、F=0は両端点に含まれる。
a_CB=0の全raw族は二つの独立な定covector Casimirを持つ。全3binaryと16Hを保持した物理式でも確認した。
これはGL不変な必要条件なので、上のrank5によりこの端点全体を除外できる。
残るb_CB=0ではLがprojective constant factorだから、必ずL=<∂y>となる。
ここで一つのflag chartだけへsourceを縮めてはいない。

商(h,x,z,w)におけるmarked line ellはDarboux一次式である。
w係数が非零なら一般cofactorのw²,hw,zw,h²係数が矛盾するので、ell∈span(h,x,z)。
ellがxでなければUFDの因数性からcofactorはxの倍数になる。
正確なeigenbasisのdetは−6a³であり、全candidateは

    x（5D eigenvalue −3a）、x+ah−2a²z（2a）、x+3ah（12a）

だけである。全一般Darboux係数とeigenbasisを独立に照合した。

CBの全標準sourceから

    D(ell)=(2/5)f_k ell,  D(k)=−f+(2/5)f_k k

を確認した。ell≠xではf_kが非零のx倍なのでWがxを含み、k=xを選べる。
上のeigenvalueをlambdaとするとf=lambda x²/4+alpha ell²となる。
D(x)=−3a x²/5との比較はlambda=12a、alpha=0だけを許す。
これは既に明示されたCB flag H=<x+3ah>, W=<h,x>である。
ell=xではD(k)−9a xk/10∈<x²>が必要であるが、全9係数行列はrank3、核<x>だけであり、独立なkが存在しない。

非零fによりratioは[1:0]で固定され、固定flagの19係数写像は単射である。
よって全proper CB fiberは、各a≠0で**正確に一つの幾何学的点**である。
別の独立監査済み全29微分のminor53747712a¹⁰/15625により、そのfiber tangentは零。
一点支持のArtin fiberは被約length1である。

propernessから各点の近傍で全sourceはfiniteとなる。reduced閉像の局所環Aはsource局所環Bへ単射であり、
B/m_A B=CだからB=A+m_A B。有限加群B/AへのNakayamaによりB=A。
従って個別reduced CBは**全a≠0でcone smooth29、射影smooth28**である。
global birationality/normalizationを先に仮定せず、a=0の型・全bad locus・full Poisson schemeを結論しない。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CA_self_pair_CB_fiber_review_20261103_v63.py --output reports/AMBIENT_CA_SELF_PAIR_CB_WHOLE_FIBER_20261103_v63_REVIEW.json --replay work/v62-review/CA-self-pair-CB-whole-fiber-replay

全行列、厳密minor、一般source恒等式、fresh certificate byte一致を検証する。全proper端点被覆は既存独立監査を明記して再使用する。
