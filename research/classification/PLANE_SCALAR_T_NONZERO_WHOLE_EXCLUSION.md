# rank零plane・scalar symbol T非零の全源は実既約quinticでは空

2026-11-08、PRIMARY E。純幾何 producer proof。新 operator／kernel／minor／GB／profile は無い。

**定理。** \(\pi\) が \(\mathbb P^4\) の rank 四 Poisson tensor、actual \(D\) が被約既約 quintic、
plane \(S\) 上で \(\pi=0\) とする。traceless co-Higgs が恒等零で scalar trace field \(T\ne0\) なら矛盾する。
全十五 Jordan source と全 parameter／kernel jump を保持する。
旧全像二十三以下を全 E と読み替えたのではなく、actual shear と全残核から新たに排除する。

## 1. 受理済み全 source と使用する exact 証明書

`PLANE_RANK_ZERO_SCALAR_FIRST_SYMBOL.md`、その
`PLANE_RANK_ZERO_SCALAR_SEPARATE_REVIEW.md` と固定／一parameter／二parameterの全 reports を入力とする。
全重複度の conormal 可換性から全35／30／9 source が有効で、全 scalar first symbol は
\[
 P_0=T\wedge E_n+K(z,w)\partial_z\wedge\partial_w,
 \quad T\in\mathfrak{sl}_3,\quad E_n=z\partial_z+w\partial_w.
\]
全 \(P_0\) は Poisson、projective rank 二以下。
\(P_1\) は canonical 全三十、\(P_2\) は全九方向。
既報の \(D_1,D_2\) は元350行に対する完全な linear Jacobi operators で、
全 kernel・nonzero maximal minor・全根は別作者 Fraction review 済みである。
今回その operators を再計算しない。

normal projection の w一 chart を \(t=z/w,x=(a,b,c)/w\) とすれば
\[
 \Pi_0=E_x\wedge(T_x+K(t)\partial_t).                            \tag{1}
\]
全 global fiber translations \(x\mapsto x+\theta(t)\)、
\(\theta\in H^0(O(1)^3)\simeq\mathbb C^6\)、は actual linear shears
\((a,b,c)\mapsto(a,b,c)+w\theta(z/w)\)。有限逆は符号反転である。

## 2. 全六 shear の stabilizer と rank 六の剛性補題

\(K\ne0\) のとき (1) の mixed anchor は \(\{x_i,t\}=x_iK(t)\)。
translation の \(P_1\) mixed 成分は \(\pm K(t)\theta_i(t)\) なので、
非零 polynomial \(K\) の乗法の単射性から六 shear の tangent map は全 \(T\) で injective。
K の単根／二重根又は T の固有値によりこの結論は跳躍しない。

\(K=0\) のとき各 normal 係数 \(v\in\mathbb C^3\) の shear kernel 条件は
\[
 v\wedge Tx+x\wedge Tv=0\qquad\hbox{for all }x.                  \tag{2}
\]
\(x=v\) より \(Tv=\mu v\)。従って (2) は \((T-\mu I)\mathbb C^3\subset\mathbb Cv\)。
\(T=\mu I+v\ell\)、\(\ell(v)=0\) と書け、trace零から \(\mu=0\)。
非零 \(T\) なら rank 一 nilpotent \(N_2\)、\(v\in\operatorname{im}T\) の場合だけが残る。
逆にその全 v は (2) を満たす。
したがって T非零では \(N_2/K0\) だけが shear stabilizer二方向を持ち、
他の全 model で六 shear は injective。

**剛性補題。** 六 shear が injective、\(\dim\ker D_1=6\)、\(\ker D_2=0\) の fiber は全 actual rank 四で空。
実際 shear の一次像は \(\ker D_1\) に含まれ dimension 六なので、全 \(P_1\) を有限 shear で零に出来る。
有限代入の \(P_2\) 補正を捨てない。変換後の weight 二 Jacobi は \(D_2P_2=0\)、従って \(P_2=0\)。
残る \(P_0\) は rank 二以下。
この論法は parameter の kernel basis の分母を反転せず、各全 fiber で次元と actual shear の単射性を使う。

## 3. 全十五 model と全複素壁の disposition

旧 complete certificates を戻すと剛性補題の例外は次だけである。

|全 symbol orbit|\(\dim\ker D_1,\dim\ker D_2\)|以下の処理|
|---|---|---|
|\(N_2/K0\)|(16,6)|旧 whole quadratic-factor theorem：空|
|\(N_3/K0\)|(8,3)|第4節|
|\(N_2/K\mathrm{nil}\)|(7,0)|第6節|
|\(T=\lambda\operatorname{diag}(1,1,-2),K=-2zw,\lambda=\pm2/3\)|(7,0)|第7節|
|regular ss／K0の \(T=\operatorname{diag}(1,-1,0)\)|(8,3)|第5節|

残る全 stratum は (6,0) であり空。
具体的には固定十型、四つの一parameter generic と全二十 root、二parameter generic、
全九 whole line とその全五十七 root／全交点を旧 accepted cover のまま使う。
二parameterの三つの excluded repeated-Jordan lines は \(\lambda D/H\) source へ戻す。
regular ss一parameterの excluded ratios一／−二／−一半も repeated sourceへ戻す。
\(\lambda=0\) はそれぞれ旧表の別 Jordan／T零 sourceへ戻し、T非零本定理に混入しない。

regular ss／K0の二 wall ratios零／−一は base permutation で一つの \(\operatorname{diag}(1,-1,0)\) 型。
\(\lambda=\pm2/3\) は normal swap \(z\leftrightarrow w\) と tensor scale \(-1\) で実際に同値。
全 wall／shear stabilizer jump を表に残し、generic (6,0) の結論をそこへ外挿していない。

## 4. N三／K零：全 residual二方向と全 P二

\(T=b\partial_a+c\partial_b\) とする。
保存された complete \(D_1\) kernel は各 normal \(z,w\) に対し次の四方向である。
\[
 R_0=c\partial_a\wedge\partial_b,\quad
 R_1=3a\partial_a\wedge\partial_b+b\partial_a\wedge\partial_c,
\]
\[
 R_2=2b\partial_a\wedge\partial_b+c\partial_a\wedge\partial_c,
 \qquad R_3=-2a\partial_a\wedge\partial_b+c\partial_b\wedge\partial_c.
\]
(2) の shear image は各 normal で \(R_0,-R_2,-R_1-2R_3\)。
これら三方向と \(R_3\) は独立なので、全 source は finite shear 後
\(P_1=\ell(z,w)R_3\)、\(\ell\) は全 binary linear form、と書ける。
\(R_3\) は strict linear Poisson tensor（\(\{a,b\}=-2a,\{b,c\}=c,\{a,c\}=0\)）だから
\([P_1,P_1]=0\)。weight 二 equation は \(D_2P_2=0\)。
保存された kernel は正確に
\(P_2=Q_2(z,w)\partial_a\wedge\partial_b\)、\(Q_2\) は全 binary quadratic。

この全 source では \(c\) に関わる raw brackets は全て cで割れる。
\(P_0\) の \(T_c=0\)、\(P_1\) の唯一の c-leg は \(c\ell\partial_b\wedge\partial_c\)、\(P_2\) に c-leg は無い。
従って hyperplane \((c=0)\) は Poisson invariant、actual Pfaffian は c因子を持つ。
Pf零なら rank 四でなく、非零なら被約既約 quintic ではない。
残る higher Jacobi equations はこの全 factor source を狭めるだけなので、それを解く必要は無い。

## 5. regular ss・零固有値／K零

\(T=a\partial_a-b\partial_b\)、zero-eigen coordinateを c とする。
六 shear は injective、保存された \(\ker D_1\) は dimension 八。
\(zR,wR\)、\(R=c\partial_a\wedge\partial_b\)、は実際に kernel に入る：
\([T,R]=0\)、\(T\wedge R=0\) から直接従う。
shear imageの ab係数は a,b だけなので、この二方向は六 shear imageから独立。
よって全 source は有限 shear 後 \(P_1=\ell R\)。
\(R\) は strict Poisson であり、\(D_2\) の全三次元 kernel は
\(Q_2(z,w)\partial_a\wedge\partial_b\)（同じ \([T,\partial_a\wedge\partial_b]=0\)、\(T\wedge\partial_a\wedge\partial_b=0\) と保存dimension三による）。
全 P は c-leg を持たず、hyperplane c零が invariant。前節と同じ actual c因子の矛盾を得る。
二つの零固有値 ratiosを一つの代表へ移す実 permutation の逆で全元 sourceへ戻る。

## 6. N二／nonzero nilpotent K

\(T=b\partial_a\)、\(K=-w^2\)。全 \(\ker D_1\) は dimension 七、shear imageは六。
保存された第3 kernel column（列を零から数える）は
\[
 R=-zb\partial_a\wedge\partial_b-zc\partial_a\wedge\partial_c
       +zw\partial_a\wedge\partial_w.                          \tag{3}
\]
これは旧三十basisの rows一／五が−一、row二十四が一、他零という literal 方向。
w一 chart の mixed anchor は \(-t^2\partial_a\)。
shear anchorは \(K\theta=-\theta\)、\(\theta\) は degree一以下なので、(3) はその像に入らない。
従って全 source は有限 shear後 \(P_1=\rho R\)、\(\rho\in\mathbb C\) 全体。

\(R=\partial_a\wedge(-zb\partial_b-zc\partial_c+zw\partial_w)\) は aに係数が依らず、strict Poisson。
\(D_2\) は injectiveなので変換後の \(P_2=0\)。
残る \(P_0+\rho R\) では非定数 ambient rational function \(b/c\) が Casimirである。
\(P_0\) の T は b,c をともに固定し、Rは b,c に同じ scalarで作用し、他の b,c brackets は零だからである。
従って全 \(\rho\) で projective rank 二以下。\(\rho=0\) を除いた主張ではない。

## 7. repeated semisimple／Hの全共鳴壁

\(T=(2/3)(a\partial_a+b\partial_b-2c\partial_c)\)、\(K=-2zw\) とする。
保存された \(\lambda=2/3\) wallの complete kernel は六shearと
\[
 R=z^2\partial_c\wedge\partial_w                               \tag{4}
\]
の直和である。実際旧三十basisの row二十だけ一の kernel columnをそのまま使う。
その w一 anchor は \(-t^3\partial_c\)、shear anchor \(-2t\theta\) の degree一／二では得られない。
全 source は \(P_1=\rho R\) へ移り、(4)は strict Poisson、\(D_2\) injectiveより \(P_2=0\)。
\(a/b\) は \(P_0+\rho R\) の非定数 rational Casimir（Tの a,b固有値は等しく、Rは両方を動かさない）。
従って全 rank 二以下で空。
\(\lambda=-2/3\) では旧 kernel の row二十九のみ一、\(R=w^2\partial_c\wedge\partial_z\) で同じ結論。
又は第3節の実 normal swap／scale とその逆で全 sourceを戻せる。

## 8. 全補集合と最終結論

旧 \(N_2/K0\) は full十六／六 kernelが raw \(\partial_a\wedge V+F_2\partial_b\wedge\partial_c\) を強制し、
actual Pfは quadratic×cubic：既受理の whole factor proofをそのまま使う。
これと第2節の全(6,0) stratum、第4–7節の全 jumpを合わせ、
**scalar \(T\ne0\) の全十五 sourceはactual被約既約quinticのrank四で空**となる。
actual finite shear と元の全 normalizer逆は PGL mapであり、rank／actual factor条件を保つ。
source係数の小さい開でしかない正規化を全源と同一視していない。

T零／K非零、weight零symbol全零、traceless非零 dependent global-V 等は別sourceである。
全rank零plane、全 E 又は総分類の完成は本定理から単独で主張しない。
新 component 数は零。旧 full reports／code／主稿／global／Gitは変更せず、同名 manifestへ全依存を固定する。
