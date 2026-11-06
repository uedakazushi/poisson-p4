# 忠実な一般法線chartと特殊底準備の独立確認

2026-11-06 / v64-classification-plan。G/E0の論理を独立に確認した。元S61のunit、一般五点尽くし、特殊底の全被覆を確認したという主張ではない。

## G：忠実平坦性

`GENERIC_NORMAL_CHART.md`の補題は正しい。Bを可換K代数、(n0,...,n41)B=Bとする。tauを超越元とすれば、各p∈Spec B上のfiberは

    kappa(p)[tau] localized at K[tau]\{0} and ell_p,
    ell_p=sum tau^j*n_j(p).

少なくとも一つのn_j(p)が非零なのでell_pは非零多項式。この整域の局所化は零環でない。多項式拡大と局所化は平坦で、全pでfiber非空だから、**B→B tensor_K K(tau)[1/ell]**は忠実平坦となる。nilpotentのないことや有限次元性を仮定しない。

B_K(tau)からの最後の局所化まで忠実平坦だとは主張していない点が重要である。固定tau値は特定の残余点を失うことがあり、この補題の代わりにならない。

既存splitter A=E5×B、(n)A=(1-e)Aを仮定すれば、ellはE5上で零、B側のnは単位idealを生成する。従って

    B=0 iff (元134式,1-u*ell)=(1) in K(tau)[q45,u]

は正しいexact criterionである。既知E5が既に開閉因子なので、B_redの空性とB=0は同値であり、余分なscheme厚みの計算を要求していない。

## G：同じ実K5上の259式

旧fixed-centerの129残差・三本pure-normal模型と、今回の実K5の130残差・四次元pure-normal空間を区別した。後者の130は元129残差と元ideal内のG6一本である。入力は

`reports/S61_SPARSE_GENERIC5_NORMAL42_EXACT_PRESENTATION_20261104_v63/typed-normal42-operator.json.gz`

の24501節点・P126・H130であり、選択五行は[50,10,52,51,15]。

`R=X(z_h n_j)+b`、`LX=I`、`P=-Lb`、`H=b+XP`を仮定すると、R=0からz_h n_j=P_hjとH=0を得る。逆に全crossとHからR=X(cross-P)+H=0が戻る。Gram逆元はparameter fieldだけにある。K=QQ(H,L,s0,q0,t0)のformal realityからdet(XᵗX)≠0を示し、複素係数へ拡大してもこの非零元は残る、という既存の議論を依存として用いる。複素未知nに正値性を課していない。

M_h=Σtau^j P_hjと置くと、259式は次の通り。

    H130,
    ell*P_hj-n_j*M_h   (h=0,1,2; j=1,...,41),
    ell² F_r(M/ell,n)  (r=1,...,5; 多項式として先に分母を払う),
    1-u*ell.

本数は130+123+5+1=259。省いたj=0はΣtau^j(ell P_hj-n_j M_h)=0から導ける。次数は2、3、4、2以下である。

逆写像z_h=u M_hでは、compatibilityとu ell=1から全126crossが戻り、Hから全130残差が戻る。五つのcleared Fとell²の可逆性から全選択五式も戻る。順写像ではΣtau^j crossによりM_h=ell z_hなので逆合成は恒等。従って全schemeの同型で、元のnormal²項やnilpotentを捨てていない。**259式は新しいunit証明書ではなく、元の空性問題の厳密な別表示である。**

## 実装の確認範囲

`classification_generic_chart.py`のprepareを読解した。PとHのrootは同じtyped operatorへ外部参照し、selected Fは元の選択行番号を保持する。`cleared_full_selected_F`は有限recipeの演算節点であり、実際の五式の全係数を展開・評価したものではない。1081単項式のloopが確認するのはz次数0/1/2の分母払いと次数上界であり、巨大F係数の新たな一致検査ではない。本文はこの射程を正しく区別している。

入力束縛はtyped operatorのhash、source_context_filesの各hash、既存producer/review reportのhashを含む。保存した`CLASSIFICATION_GENERIC_CHART_20261106.json`の全context/proof hashを独立に現物と照合し、一致を確認した。

probeは保存済みQQ中心を法32003へ落とした診断である。tauを超越のまま保つが、betaと標数は固定されている。仮にunitが出ても標数零の一般unitを結論しない、という制限は適切である。今回の独立確認ではprobe・GB・新しいsigned係数を実行していない。

## E0：底分母に必要な因子

`EXCEPTIONAL_BASE.md`の必要条件は正しい。Dを既知三優越成分の和とし、全X上でpi^-1(D(f))⊂Dと仮定する。他の既約成分C_iについてfがその閉像Z_iで恒等零でなければ、C_i内の非空開pi^-1(D(f))と非空開C_i\Dが交わって矛盾する。従ってf∈I(Z_i)。NEW25の複素既約像がV(G28)であるという既存定理を用いるとG28|fが従う。

これはG28を反転すれば全未知枝が消えるという十分条件ではない。49閉像の和が全特殊底であるとも主張しない。NEW25像のDelta非零点があるという既存証拠から、特殊底をV(Delta)だけに限定できない。

DAG graphのelimination idealが閉像idealとなる証明も、源が既約、採用openが非空dense、DPMのsqrt3埋込みを一つ固定、という記した仮定の下で正しい。beta祖先以外の節点を除くと定義域がdense openへ広がり得るが、閉像は変わらない。源分母を底分母と同一視しない点も適切である。

今回読解したE0実装は源の逆元・祖先・hash・元134式のbeta=0制限を保存するものであり、像のimplicitizationや全fiber被覆を行ったと誤表示していない。52件の既存像次元・成分性の証明は引き続き既存台帳と独立reviewへの依存である。


## 主原稿と小行列証明書

主原稿`poisson_p4_classification_ja.tex`の忠実chart、259式の逆写像、特殊像の分母必要条件を本ノートと照合し、数学的な内容の相違は見つからなかった。主分類未証明・既存splitter依存・unit未計算・有限体診断の限界が明示されている。座標番号の1始まりと0始まりは原稿全体で統一する必要がある。

保存したC0証明書の16×16橋と二つの16角接行列について、runnerのdomain-ge法とは別にBareiss法で行列式を再計算した。順に36/625、8428695240、137438953472/996451875で一致した。橋の両側逆積と証明書内の全入力・runner hashも一致した。これは小行列証明書の点検であり、全134式や分類の再検証ではない。
