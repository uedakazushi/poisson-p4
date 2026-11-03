# Q_b の全 CA/CB/N 旗：全非零parameterの独立監査

2026-11-03。`AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63` を独立に受理する。基本packetで得た一般bの二枝結論を、**全ての複素b≠0**へ広げるための、全proper旗の被覆証明である。先行CAのP_aとは別の族Q_bを使う。

## 1. 全端点と内在的な一意因子

独立コード `code/ambient_CB_self_pair_all_nonzero_flags_review_20261103_v63.py` は作者 helperを import せず、実際の `Q_b=∂y∧V_b` から次の全行列を作り直す。

* 定covector contractionの全75×5成分。選択5 minorは `−384b⁴/625`。
* `Euler∧v∧Q_b` の全175×5成分。ey列は零、選択4 minorは `2b`。

従って全b≠0で定covector Casimirはなく、定vector factorは厳密にL=eyの一つだけである。

CAとCBのproper模型のrank2端点の全被覆は受理済みのglobal模型を明示依存とする。独立計算では、全16個のHamiltonian係数を保った普遍式から、CAのf=0端点とCBのa_CB=0端点が二つの独立な定covector Casimirを持つことを確認する。したがってこれらの端点とF=0の共通境界は全て不可能である。残るCA g=0、CB b_CB=0には定vector factor Lがあり、上で一意にeyと定まる。よって候補は全てann(ey)の四変数(h,x,z,w)の旗であり、未検証の小chartだけを列挙していない。

## 2. 全 Darboux 一次式

`D=5V_b` とし、任意の非零一次式を

\[
\ell=Ah+Bx+Cz+Ew
\]

とする。`D(ell)=ell(r_hh+r_xx+r_zz+r_ww)` の全係数を検証する。E≠0ならw²、hw、zw係数から順にr_w=r_h=r_z=0となり、h²係数は5E=0を要求して矛盾する。従ってE=0である。

残るD(ell)はxの倍数である。ellがxの倍数でなければUFD上でellとxは互いに素なので、cofactorはλx、ellは三変数の線形作用のeigenformになる。ell=xの場合も同じeigenformに含まれる。全eigenbasisは

\[
x\ (-2b),\qquad K=x+bh-b^2z\ (3b),\qquad H=x+2bh\ (8b).
\]

eigenbasisの行列式は−2b³で、全b≠0で三eigenvalueは相異なる。このため他の不変一次式は存在しない。これはsamplingやDarboux方程式の有限degree探索ではなく、任意の一次式の完全な場合分けである。

## 3. 全 W の決定

普遍raw式をcanonical化し、全16 H係数を記号のまま残して次の恒等式を再構成する。marked lineをell、Wの別の生成元をk、binary sourceをfとすると

\[
\begin{array}{lll}
CA:&V(\ell)=f_k\ell/5,&V(k)=-f+f_kk/5,\\
CB:&V(\ell)=2f_k\ell/5,&V(k)=-f+2f_kk/5.
\end{array}
\]

ell=K,Hのときλ≠0なので、第一式のcofactorからx∈Wとなり、k=xを取れる。fをkで積分して第二式へ入れると、CBでは `f=λx²/4+αell²` と `λ=8b,α=0`、CAでは `f=λx²/2+αell²` と `λ=4b/3,α=0` が必要になる。ell²とx²の独立性も全b≠0で確認する。従ってCBではHだけ、CAではこの二候補はどちらも不可能である。

ell=xの場合はkの**全四係数**を保持する。x²係数以外の9係数を並べた9×4行列を直接構成し、保存した全成分と照合する。

* CB: `(D−3bx)k∈<x²>`。rank2のminorは−25b、独立なkernel二列はxとh−bz。従ってW=span(x,K)だけ。
* CA: `(D−8bx)k∈<x²>`。rank2のminorは−25、kernelはspan(h,x)。従ってW=span(h,x)だけ。

以上から全CB geometric fiberは既に実sourceで構成済みの二旗、全CA geometric fiberは元の一旗に限られる。固定旗の線形embeddingはinjectiveで、CBでは残るratioは[1:0]で一意なので、係数の異なる余分なsource点もない。

## 4. N、全scheme fiber、smoothness

独立に再構成したCA fixed20と、旧N fixed19がindex2を除く閉subbundleであることを全150係数で確認する。Q_bのこの追加係数は零である。Nの旗は同じCA旗への包含を持つので、その全fiberも元の一旗だけである。

CA全150×29、N全150×28のsource微分を元の物理GL作用から作り直す。固定fiberを除いた全150×9行列とminorを照合し、それぞれ

\[
-262144b^9/15625,\qquad 1572864b^9/78125
\]

を得る。従って全b≠0でrank29、28である。CBの両rank29、同じQにおける全58接方向のunion37 / intersection21は、独立受理済みの基本reviewをそのまま用いる。

いずれのsourceでもfiberの接空間は零であり、全geometric fiberは有限である。propernessによりその点の近傍で有限化し、有限local algebraのm=m²にNakayamaを使うと各fiber点は被約である。従って**全b≠0**で、CA/Nは一台長1、CBは二つの被約点を持つ。

CA/Nは唯一のsmooth sourceを全fiberとして持ち、source微分はfull rankである。ambient completed ringからそのsource completionへの全射と、全proper imageとの一致により、CAはsmooth cone29、Nはsmooth cone28となる。射影版では28、27次元である。

## 5. CBの全b≠0の二枝結論

基本reviewの全proper image argumentを、今回の全fiber二点へ適用できる。全b≠0でreduced CBは正確に二つのsmooth29 formal branchを持ち、その接空間の和は37、交わりは21となる。meeting schemeのedim21と完全列から、cone depth≤22<29、射影depth≤21となり、個別reduced CBは全b≠0でsingular・nonnormal・non-CMである。

whole imageのedimを37と同定せず、meetingのclean性や正確なidealも推論しない。b=0の全odd ideal・conductor・完成型、全global bad locus、actual CA/CB intersection Dは別問題として保持する。新しいPoisson成分は数えない。

## 再現と保存

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/ambient_CB_self_pair_all_nonzero_flags_review_20261103_v63.py \
  --output work/cb-Qb-all-flags-independent-review
```

保存先は `reports/AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63_REVIEW/report.json`。新しい全係数検算、基本Q_b review、global proper模型、N閉subbundleの既存証明書をhashで区別して結ぶ。GBは用いない。
