# 一つのcontrolled generic Schur chartによる未知残余の正確な空性判定

2026-11-03。新global monogenic quinticから、**一つの一般補助定数aと一つの4×4 minorによる局所化で、未知残余Bの空性を正確に判定する**。minorはB上の非零因子で、全associated/nilpotent submodulesが局所化後も残る。単なる「非空cyclic openがある」という仮定ではない。scalar unitはまだ計算していない。

## 1. 実sourceと一変数fiber

R'=K[n42,S1,S2,H]、C=R'[h]/G(h)、deg G=5とする。全134remainder r_iからA=C/(r_i)を得る。known-five splitとnormal閉埋込みより、B上ではR'→Bが全射である。

Bの幾何点pではA tensor k(p)=k(p)で、hは一つの値h0を取る。従ってk[h]/G内のideal(r_i)は正確に(h−h0)で、次元4である。Gが重根を持ってもよい。各r_iの次数≤4なのでr_i=(h−h0)s_iと一意に書け、s_iはD=k[h]/g4、g4=G/(h−h0)のunit idealを生成する。

## 2. moment polynomialと一つのpivot

    r(a)=sum_(i=0)^133 a^i*r_i,
    b_j(a)=(T_h−aI)T_h^j*r(a), j=0..3

とし、四本のb_jの固定bottom四行のdeterminantをdelta(a)とする。b_jは元5×670 matrixのcolumnsの明示的な定数結合で、全sourceを保持する。各entryのformal parameter次数は≤2である。

fiber上ではr(a)=(h−h0)S(a,h)である。ideal(h−h0)のbasis (h−h0)h^j(j=0..3)からbottom四係数への変換はdeterminant1を持つ。従って厳密に

    delta(a)=g4(a)*Norm_D(S(a,h))

となる。Dは非被約でもよい。s_iがDを生成するので、generic aでSはunitとなる。具体的には他の根での値と、h0が重根ならh0での一次微分の、最大四つの非零linear条件を避ければよい。それぞれaの次数≤133だからbad-generator値は≤532。g4(a)がさらに≤4個を除き、**全ての未知幾何点pでdelta(a,p)は非零多項式、次数≤536**である。

## 3. 非零因子と全nilpotentの保持

delta(a)の係数が生成するBのidealはunitである。もし真のidealなら、そのmaximal ideal上の幾何点で全係数が零となり、第2節に反する。

従ってdelta(a)はB[a]の全associated prime p[a]を避け、非零因子である（同値にMcCoyのprimitive polynomial criterionを使える）。faithful field extension K(a)後も非零因子であり、

    B -> B tensor K(a) -> B(a)[delta^-1]

は単射となる。したがってembedded/nilpotent submoduleも消えない。delta=0にある閉点をchartが全て覆うとは言わないが、schemeの情報はschematically denseな局所化へ単射で残り、idealはdelta-saturatedである。

known E5では全formal parametersが零でPhi entriesも零なのでdelta=0である。よってA(a)_delta=B(a)_deltaであり、

    B=0 iff A(a)_delta=0

が正確に成立する。この空性同値はdim B≤2を必要としない。

## 4. 一つのscalar unit target

四pivot columnを元670column matrixへ追加してもimageは変わらない。deltaを反転したSchur eliminationから、

    g_l(a)=det[b_0(a),b_1(a),b_2(a),b_3(a),Phi_column_l],
    A(a)_delta=R'_K(a)[delta^-1]/(g_0,...,g_669)

である。deltaのformal parameter次数≤8、g_l≤10。a次数は両方≤536である。従って全generic-five問題はこの**一つのcontrolled scalar unit**に同値である。実unit証明、QQ cofactor、nilpotenceは未計算である。

modpでunitまたはnilpotenceを見ただけではgeneric QQに移さない。今回証明したのは局所化の非零因子性と判定の同値であり、空性そのものではない。

## 5. 全点のcoverが必要な場合

既存dim B≤2を使えば、三つの独立generic constants a0,a1,a2で三chartがB全体を覆う。全三deltaが零となるincidenceはB上にfinite fibersを持つので次元≤2、A3のgeneric pointを含まない。nilpotentを含む生成はNakayamaで保たれる。

また0,1,...,536という**537個の固定QQ constants**を取れば、次数≤536の非零多項式はその全てで零にはならないため、全Bを覆う。こちらはdimension boundを不要とする。特定の三つのrational constantsが検証済みとは主張しない。

## 6. 五次元の正確なfixture

producerはsimple、multiple-root、およびG=h^5、h0=0、r0=h²、r133=hの三例で、全5×4 columnsと上のresultant/norm identityをQQ[a]で検査する。最後はdelta(a)=a^536となり、次数536の上界が実際に達する。これらは小さいalgebraic implementation checksで、元S61のunit標本ではない。

保存: `code/s61_affine_normal_moment_schur_20261103_v63.py` と `reports/S61_AFFINE_NORMAL_MOMENT_SCHUR_20261103_v63/`。全134sourceと既証明閉埋込みをpinし、finite companion-column recipeを指定する。未知残余の実空性は未完了である。
