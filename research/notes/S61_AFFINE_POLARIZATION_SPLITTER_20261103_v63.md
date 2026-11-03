# 二次式のpolarizationによる既知五点の大域idempotent分離

2026-11-03。一般betaの全アフィンS61代数から、既知五点の因子を
**未知qによる局所化なしに、degree≤45の明示idempotentで大域的に分離**する。
元134式の全16betaを保ったQQ polarization、25組のmidpoint Jacobian、
分離linear formの判別式とJacobian normの単元性を有限計算で確認した。
一般係数のidempotentはdeterminant・有限逆行列・traceという有限recipeを指定するもので、
巨大な多項式係数表を展開保存した結果ではない。
未知残余の空性・一般次元零・全長5・全成分の尽くしは示さない。

## 1. 全アフィン代数と、実商としての既知五点代数

係数体Kを一般betaの標数零関数体とし、元45未知数をq、元134本の二次以下の式をFとする。
A=K[q]/I、I=(F)は全アフィンfiberの実代数である。
既存P29有限fiberの証明から、**実商写像**phi:A→Eがあり、Eは既知五点の有限étale5代数である。
これは単なる五点の数え上げから仮定した商ではない。
元134式の斉次有限17底模型でY=vを課し、t=1が成り立つ有限自由5代数として構成された。
有限族へのR^5全射と同じ既知五点によるrank下界の両方を使った証明である。
参照：S61_AFFINE_SEVENTEEN_BASE_LIFT_20261103_v62.md §4、
S61_AFFINE_QUINTIC_OPEN_20261103_v62.md §3。
今回その大きな有限模型を再計算せず、凍結証明と全input/evidence hashを再確認した。

45個の元行を選ぶselectorはS61_AFFINE_LINEAR_BLOCK_20261103_v63の
R10行とcritical-selected35行の連結である。row combinationではなく元行の選択そのものを使う。
専用producerは中心betaの既知五点を元134式へ代入してから、
選択45×45Jacobianが全五点で可逆であることを直接確認した。
従って一般beta上でd(p)=det J_selected(p)はEの単元である。
必要ならNorm_E/K(d)を反転したbetaだけの非空開で述べる。

E内の分離元はh=Σ_(j=0)^44 (j+1)q_jとする。
中心mod32003で五つのh値は相異なり、Vandermondeは非零。
Eの1,h,…,h^4基底行列式とその判別式を反転すれば、Eはmonogenic有限étale5代数になる。
これらもbetaだけの局所化であり、未知qの支持を切るものではない。
E内のp_j、d(p)^(-1)、traceはこの五元基底の有限行列演算で定義できる。

## 2. 元式の正確なpolarizationとadjugate

E係数の普遍既知点をp=(p1,…,p45)と書く。F(p)=0である。
全Fがq次数2以下なので、多項式恒等式

    F(q)−F(p)=J_F((q+p)/2)*(q−p)

が成り立つ。producerは元データの全16betaを記号のまま保ち、
QQ Fraction演算で全134行の全polarization係数を照合した。
これは中心betaだけの恒等式ではない。

選択45行についてM_p(q)=J_selected((q+p)/2)、D_p(q)=det M_p(q)と置く。
adjugateを掛けると、元134式のidealに対する正確な所属

    D_p(q)*(q_j−p_j) ∈ I*E[q]  (j=1,…,45)

を得る。実際、右辺のcofactorはadj(M_p)の各行と選択45式との積である。
D_pのq次数は45以下で、D_p(p)=d(p)はEの単元である。
行列式・adjugateは有限recipeとして完全に指定されており、係数展開は必要ない。

e_p(q)=D_p(q)/d(p)と置く。
通常の多項式差商により、任意f∈E[q]についてA⊗Eの中で

    e_p*f(q)=f(p)*e_p。

f=e_pとするとe_p(p)=1なのでe_p²=e_pを得る。

## 3. 五点への分裂、直交性、traceによる降下

Eを分裂させる忠実平坦な有限étale体拡大後でpを五点p_iへ分ける。
異なるi,jについて

    M_(p_i)(p_j)*(p_j−p_i)=0。

より、adjugateを掛ければD_(p_i)(p_j)は全座標差を消す。
分離元hの値の差が単元なので、D_(p_i)(p_j)=0となる。
よってe_i(p_j)=delta_ij、e_i*e_j=0 (i≠j)である。
producerは中心の全25組について、全134行のmidpoint kernel恒等式を検査し、
選択行列式が非零対角・零非対角になることも直接計算した。

ここで

    e(q)=Tr_(E/K)(D_p(q)/d(p))

と定義する。分裂後はe=Σ_i e_iなので降下し、A内の大域idempotentになる。
**degree_q e≤45**である。traceを5で割る正規化は行わない。

さらにf∈Eをsigma(f)=Σ_i f(p_i)e_iで送る写像が降下する。
これは単位元をeとするeAへの単位的環準同型で、
phi sigma=id_E、sigma phi(a)=eaを満たす。
従ってスキームとして正確に

    A ≅ E × A_res,   A_res=A/(e)

を得る。既知五点idealはA内で(1−e)、未知残余idealは(e)である。
元多項式環では

    I_known5=I+(1−e),  I_res=I+(e),
    I=I_known5 ∩ I_res,  I_known5+I_res=(1)

というcomaximalな大域分解になる。
cyclic determinantやquintic Hの零集合の外だけの主張ではない。
基礎beta環上でも、Eの有限étale性・dのnorm・h基底の行列式を反転して同じ証明が働く。
有限型A自身がbeta上有限だという仮定は用いない。

## 4. 残余の射程、保存、次の課題

このidempotentは未知の孤立点、正次元成分、非被約・埋め込み構造を
A_resにそのまま保持する。有限性や空性は新たに仮定しない。
producerの小さいQQ検査例では

    I=(x²−x,xy,y²), p=(1,0),
    e=x*((x+1)/2+y) ≡ x mod I,
    A_res=QQ[y]/(y²)

となり、非被約な未知残余も消さずに残る。

producer `code/s61_affine_polarization_splitter_20261103_v63.py`、
receipt `reports/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63/`。
全134QQ恒等式、25 midpoint determinant、全134 kernel検査、
45元row selector、hと五つの値、norm・判別式、構成recipeを保存した。
一般betaのdegree45多項式eを全展開したとの主張ではない。
専用巨大GB・1031乗算行列の再計算は行っていない。

次の尽くし課題は**この実残余ideal I+(e)=(1)**を示すことへ正確に移せる。
未解決の全支持はこの商に残り、既知五点を含む大きなquintic H補集合を
別に議論する必要なく、有限degreeの一つの大域式で既知因子と分離できる。
残余の空性は未証明である。
