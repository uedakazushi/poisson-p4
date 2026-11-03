# 実self-pair族の全CA/N旗：全a≠0への強化

2026-11-03。既存のgeneric finite-open結論を上書きせず、新しい全旗判定により**全a≠0**へ強化する。
対象族は先行QQ[a]実族P(a)。CB全旗packetと共有する75×5Casimir、175×5constant factor、三Darboux一次式の証明を入力とする。

## CA全proper端点の判定

CAの全階数二preimageはg=0またはf=0で覆われる。
f=0端点はcanonical化しても二つの独立な定covector Casimirを持つ。
P(a)の定covector contraction行列は全a≠0でrank5なので、この端点は不可能である。

残るg=0の全旗について、唯一のintrinsic constant factorからL=e_y。
四変数tracefree場Cの全不変一次式は、先行の全係数判定により正確に

    ell=x, H1=x+ah−2a²z, H2=x+3ah

である。5Cのeigenvalueはそれぞれ−3a,2a,12a。
CAのcanonical raw式から、全16Hamiltonian係数を残したまま

    C(ell)=f_k*ell/5,
    C(k)=−f+f_k*k/5

が成立することを確認した。

ell=H1またはH2ならx∈Wであり、k=xを取れる。
eigenvalue lambdaに対してf=lambda*x²/2+alpha*ell²となるので、
C(x)=−3a*x²/5との比較からlambda=2a、alpha=0だけが可能である。
従ってH1が第二CA旗を与え、H2は不可能。

ell=xの場合はf=alpha*x²−3a*x*kであり、

    C(k)−(12a/5)*x*k ∈ <x²>

が必要である。x²成分を除く全9×4行列のrankは2、核は正確にspan(h,x)。
従ってW=span(h,x)だけが可能で、これは第一CA旗である。

以上から、任意のa≠0でCA全geometric fiberの旗は厳密に次の二つだけである。

    F1: L=e_y, W=span(h,x), H=span(x),
    F2: L=e_y, W=span(x,h−2az), H=span(x+ah−2a²z).

実二sourceは既存の全150係数多項式恒等式で構成済みで、各固定旗でのlinear fiber係数は一意。
両全29 source微分のminorはaの単項式で、全a≠0でimmersionである。
したがって全CA fiberはschemeとしても二つの被約点である。
この結論には、特定のa=0近傍へ縮めたfinite-open仮説はもう不要である。

## Nの全fiberとsmoothness

N sourceは同じ旗base上でCA fixed-fiberのC[2]を除く19係数subspaceであり、全parabolic作用で保たれる。
従って任意のN sourceは、同じ旗を持つCA sourceへ入る。
二つの既知CA sourceの欠けた係数は、それぞれ0と−a²である。
ゆえにNの全geometric fiberは第一旗だけである。

Nの実19係数と全flag9から全150×28 source微分を独立に再計算した。
固定19次単元blockを除いた9次minorは

    71663616*a^9/78125.

従って全a≠0でrank28、全fiberの接空間は零。
propernessとfinite fiber、およびlength1のNakayama同型により、Nはこの全parameter族でsmooth cone28（projective27）である。
第二CA旗を第二N旗と呼ぶことはできない。

## CAの局所帰結

全CA fiberが二被約点で両sourceがsmooth29 immersionなので、先行の二branch証明は全a≠0に適用できる。
実reduced CAのcompleted ringは二つのsmooth29 branchの実imageであり、

    A=T/(I1∩I2),  T/Ii≅C[[29]].

先行の同一P上でのcross-tangent rankは全a≠0で37、intersection tangent次元21。
全scheme intersectionをC=T/(I1+I2)とするとedim C=21、従ってdepth C≤21。
exact sequence0→A→C[[29]]⊕C[[29]]→C→0からdepth A≤22<29。
したがってreduced CAは全a≠0でnonnormal・singular・non-CMである。
射影branchはsmooth28、depth上限21。
whole CAのembedding dimensionを37と同定したり、intersectionのclean性を主張したりしない。

別の全CB旗packetでは、この同じ族でCBはsmoothである。
この非対称性は、個別CAの二旗と個別CBの一旗という実sourceの差であり、R19のacross-pairの転用ではない。
原点a=0の全conductor・全odd ideal・全completionは依然このpacketの範囲外である。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/ambient_CA_self_pair_all_nonzero_flags_20261103_v63.py --output <new-output-directory>

全係数とminor:
`reports/AMBIENT_CA_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63/`。
実行約1秒、制限60秒・2.2GB。新しいGBなし。
