# CA self-pair族の全CB旗：一つだけで、CBはsmooth

2026-11-03。三つの不変一次式の一つから第二CB旗を探したが、実際には既存のCB旗だけが可能だった。
単一candidateの失敗から全fiberを推測せず、**proper CB模型の全階数二端点を覆う内在的判定**を行う。
結論は全a≠0で、全CB fiberは被約一点、CBはそこでsmoothである。
CAの二旗をCBの二旗に読み替えない。

## 入力と実際の既存CB旗

独立監査済み実テンソルはP(a)=∂y∧C(a)で、

    5C=(5x²+12a*hx)∂h−3a*x²∂x
        +(5hx+2a*xz)∂z+(5h²+5z²−8a*xw)∂w.

射影Poisson rankは全aで2。
既存の明示CB表示は変換h_new=−x−3ah, x_new=xによるものであり、元の物理座標での旗は

    L=<e_y>, W=span(h,x), marked H=span(x+3ah).

従ってx+3ahを再発見しても第二旗にはならない。
指定CB liftの全29 source微分は、別の独立監査済みpacketで全a≠0に対しrank29と分かっている。

## 全proper端点の絞り込み

全CB proper模型の階数二preimageは、a_CB=0またはb_CB=0、またはF=0で覆われる。
F=0は両端点に含まれる。
これを一つの旗chartの近傍だけへ制限しない。

a_CB=0端点の標準raw式はrawCB(0,g(h,x),H)である。
canonical化してもh,xの二つの独立な定covector Casimirを持つことを、全3binary+16H係数のまま確認した。
この性質はGL移動で保たれる。
一方P(a)の全75×5定covector contraction行列はrank5で、minorは

    −864*a^4/625.

従って全a≠0で、この端点とH-only端点には入らない。

残るb_CB=0端点は、flagの定ベクトルLを常にprojective constant factorとして持つ。
全175係数の(Euler∧v∧P)行列はrank4、minor3a、核は厳密に<e_y>である。
従って任意のCB表示のLはe_yである。
全ての残る旗は、商四変数(h,x,z,w)の二次元Wと、そのmarked line ellで判定できる。

## 全不変一次式の分類

b_CB=0のcanonical tracefree場に対し、marked ellと第二W生成元kは

    C(ell)=(2/5)f_k*ell,
    C(k)=−f+(2/5)f_k*k

を満たす。Hの全16係数を含む標準raw式からこの恒等式も確認した。
特にellはCのDarboux一次式である。

ell=A*h+B*x+Cz*z+D*w、cofactor=rh*h+rx*x+rz*z+rw*wとし、
(5C)(ell)=cofactor*ellを全二次係数で比較する。
D≠0ならw²係数からrw=0、hwとzwからrh=rz=0、最後にh²係数からD=0となり矛盾する。
従ってellはspan(h,x,z)に属する。
この部分では(5C)(ell)=x*A_matrix*ellで、ellがxの倍でなければUFDの因数性によりcofactorはxの定数倍である。
a≠0での三つの相異なるeigenvalueと一次式は正確に

    x:                  −3a,
    H1=x+ah−2a²z:       2a,
    H2=x+3ah:           12a.

三eigenform行列のdetは−6a³である。
従って他のmarked ellを落としていない。

## 三candidateの全W判定

ellがH1またはH2なら、f_kがxの非零倍なのでx∈Wとなり、k=xを選べる。
eigenvalueをlambdaとすると

    f=lambda*x²/4+alpha*ell²,
    C(x)=−alpha*ell²−lambda*x²/20.

実際のC(x)=−3a*x²/5と比較すると、ell²とx²の独立性からlambda=12a、alpha=0が必要。
H1は不可能で、H2だけが残る。
この場合f=3a*x²、W=span(H2,x)=span(h,x)で、既存の実CB表示と全く同じ旗である。

ell=xの場合は

    f=alpha*x²−(3a/2)*x*k,
    C(k)−(9a/10)*x*k ∈ <x²>

が必要である。
k=A*h+B*x+Cz*z+D*wとし、右辺のx²成分を除いた全9×4係数行列を計算するとrank3、核は<x>だけである。
具体的にはh²/z²でD=0、xzでCz=0、hxでA=0。
従ってellと独立なkが存在せず、このcandidateも不可能である。

以上から全a≠0で、全proper CB geometric fiberは既存の**一つの旗だけ**である。
非零fによりratioは[1:0]、固定旗での19係数行列はinjectiveなので、同じ旗にsource係数の余分な解もない。

## 被約fiberとsmoothness

全proper fiberが有限一点なので、そのtarget近傍へ縮めるとsource mapはfiniteである。
独立監査済み29次source微分のminor53747712*a^10/15625が非零なのでfiber tangentは零。
一点支持のArtin fiberはNakayamaによりlength1である。
finite local source代数は1で生成され、実reduced imageの局所環からの写像はsurjectiveかつinjectiveとなる。
従ってtarget CB局所環はsmooth29 sourceと同型である。射影ではsmooth28。
ここでglobal birationalityやglobal normalizationを先に仮定しない。

この族ではCB self-pairを構成できないことが全proper fiberとして分かった。
CB原点の全conductorや全bad locus、他のCB self-pair族、full Poisson schemeについては結論しない。
CA側の二smooth branch・non-CM性と、CB側のsmoothnessは両立する。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/ambient_CA_self_pair_CB_whole_fiber_20261103_v63.py --output <new-output-directory>

全係数・linear rank witness:
`reports/AMBIENT_CA_SELF_PAIR_CB_WHOLE_FIBER_20261103_v63/`。
実行約0.3秒、制限120秒・2.2GB、GBなし。
