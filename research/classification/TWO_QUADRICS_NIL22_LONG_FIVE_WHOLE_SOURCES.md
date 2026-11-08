# H＋相異二既約quadric：nil22 長 translation の全二十源を五つの実正常族へ

2026-11-08、PRIMARY E。純係数・因子幾何による全支持還元。
新しい operator、kernel、minor、profile、GB、jet 又は CAS は無い。
元二十係数の全 source を使い、既 H＋既約quartic 七源の定理を因子型外へ適用しない。
新 component の増分は零。以下の四つの未認識軌道を既知親と断定しない。

## 1. 全 potential と、因子条件からの二つの Hessian rank

固定 cofactor は \(Z=y\partial_x+\partial_y+w\partial_z\)。
`HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_SOURCES.md` とその独立 V2 review の
**全二十係数 potential 恒等式**を入力にする。
補助座標 \((t,u,v,w)=(y,x-y^2/2,z-wy,w)\) で
\(P=\partial_t\wedge V+\operatorname{Jac}\Phi\)、\(\operatorname{div}V=-1\)、
全 Jacobi は \(V\Phi+\Phi=k\) と同値で、actual affine Pfaffian は \(F=\Phi-k\)、\(VF=-F\)。
この多項式座標を PGL 同値とは呼ばない。下の正常化は全て元座標の affine 変換で行う。

\(D=hq_1q_2\)、相異二既約 quadric を仮定する。
nilpotent \(Z\) は各 factor を保ち、二次式上の nilpotent 作用に固有値は零しかないので \(Zq_i=0\)。
次数二以下の完全核を単項式係数で解けば
\[
 q_i=a_i u+b_i v+c_{i0}+c_{i1}w+c_{i2}w^2.
\]
斉次表示は \(\operatorname{span}\{h^2,hw,w^2,hx-y^2/2,hz-yw\}\)。
\((a_i,b_i)=0\) なら \(q_i\) は \((h,w)\) の binary quadratic で可約だから、両 row は非零である。

従って \(H=\operatorname{Hess}_{u,v}F\) は rank 一又は二。
rank 零は本因子型では不可能である。
独立 row なら rank 二、従属 row なら rank 一という積の Hessian の直接計算であり、
未知 source parameter の一般性は仮定しない。

\(Vq_i\) は fiber 次数一以下、\(q_i\mid Vq_i\) なので \(Vq_i=\rho_i(w)q_i\)。
fiber 係数行列 \(M\) は定数なので \(\rho_i\) は実際には定数で、\(\rho_1+\rho_2=-1\)。
この因子別 semi-invariance は \(VF=-F\) と互いに素性からの必要条件である。

## 2. Rank 一では q=−1/4 だけが残る

\(q=n_{14}\)、\(s=n_{13}\) と書く。
Hessian rank 一の場合、二 factor の fiber row は平行である。
active linear form \(U\) を共通 fiber 項に取ると \(q_i=U+p_i(w)\)。
conformal 恒等式 \(M^TH+HM+H=0\) は active eigenvalue を \(-1/2\) にするので
\(\rho_1=\rho_2=-1/2\)。相異性から \(D(w)=p_1-p_2\ne0\) であり、
\[
 (s+qw)D'=-D/2.
\]
\(D\) は定数ではあり得ず、degree は一又は二。
最高係数から \(q=-1/2\) 又は \(-1/4\) となる。
両値とも非零なので、既実 affine centralizer
\(z'=z+(s/q)y,\ w'=w+s/q\) により \(s'=0\) へ送れる。
逆も同じ係数の減算であり、空間 chart を除かない。

### q=−1/2 の全補集合は空

\(n_{17}\ne0\) なら既 universal actual action
\(x'=x+az,y'=y+aw,z'=rz,w'=rw\) を使って \(n_{17}=1,n_{18}=0\) と出来る。
Hessian rank 一と conformal 式は
\(n_{11}=-1/2,n_{19}=-5/4,n_{12}=-3n_8\) を強制する。
元 core 行25,12はそれぞれ \(n_2=-3/16\)、\(n_8^2=1/16\) を与える。
しかし二 factor が同じ active \(u\) だけに依存するには、全 \(F\) の passive \(v\) 係数
\(-(n_6+n_7w+n_8w^2)\) が零でなければならず矛盾する。

\(n_{17}=0\) なら Hessian rank 一から \(n_{18}=0\)、active 方向は \(v\)。
\(F\) の \(u\) 係数零は \(n_{19}=q\) を強制し、
active eigenvalue は \(n_{11}=q+1/2=0\) を与える。
すると \(H_{22}=-3n_{11}/2+n_{19}+1/2=0\) となって rank 一に反する。
従って q=−1/2 の全 source が本対象では空であり、以下は q=−1/4 だけ。

## 3. Rank 一・active u：一つの whole 族、全 T17

\(n_{17}\ne0\) を同じ実 action で \(n_{17}=1,n_{18}=0\) にする。
conformal／rank 一より
\[
 n_{11}=0,\quad n_{19}=-1/2,\quad n_{12}=-3n_8.
\]
元 core 行25,12から \(n_2=n_8=n_{12}=0\)。
passive \(v\) 係数零から \(n_6=n_7=0\)。
\(s=0\) にした全残 Jacobi は、\(L_0=n_{15},L_1=n_{16}\) として
\[
 n_0=-L_0/4,\ n_1=-L_1/8,\ n_3=-L_0L_1/2,
 \ n_4=-L_1^2/2,\ n_5=L_1/4
\]
だけを課す。\(n_9,n_{10}\) は任意に保持する。
actual \(x'=x-L_0/2-L_1w/2\)、\(z'=z-4n_9\) により、
\(B=n_{10}\) 任意の次の whole 族 \(U(B)\) となる。
\[
\begin{array}{lll}
 P_{xy}=R/2,&P_{xz}=y(Bw-z/4),&P_{xw}=-yw/4,\\
 P_{yz}=Bw-z/4,&P_{yw}=-w/4,&P_{zw}=-2R-w^2/2,
\end{array}\qquad R=x-y^2/2.
\]
全 \(F=-R(R+w^2/4)\) で、二 factor は rank 三と四の既約 quadric。
上の actual 変換は全て明示逆を持ち、自由 \(B\) の値を一つに選んでいない。

真の completed T17 helper を
\[
 F_T=pAB+uH^2,\quad G_T=qCD+vF_T,\quad
 \text{cross}=(A\partial_A-B\partial_B)\wedge(aC\partial_C+bD\partial_D),\quad a+b=p
\]
と書く。全五 parameter の named 同定は
`TWO_QUADRICS_CONFLUENCE_NAMED_MPS_PARENTS.md` の既独立入力である。
\(e\ne0\) で actual frame と parameter を次に取る。
\[
\begin{aligned}
 &H=h,\ A=h+ey+e^2x,\ B_T=h-ey+e^2x,\ C=w-ez,\ D=w+ez,\\
 &p=-e/2,\ u=e/2,\ q=e,\ v=-4/e^2,\ a=-e/4-Be^2,\ b=-e/4+Be^2.
\end{aligned}
\]
frame determinant は \(4e^4\)。逆は
\(x=(A+B_T-2H)/(2e^2),y=(A-B_T)/(2e),z=(D-C)/(2e),w=(C+D)/2\)。
\(R_e=hx-y^2/2+e^2x^2/2\) とすれば
\(F_T=-e^3R_e\)、\(G_T=e(4R_e+w^2-e^2z^2)\)。
period pole と inverse-frame pole を全六括弧で消すと、target tensor は正確に
\[
\begin{array}{ll}
 P_{xy}=R_e/2,&P_{xz}=y(Bw-z/4),\\
 P_{xw}=-yw/4+Be^2yz,&P_{yz}=(h+e^2x)(Bw-z/4),\\
 P_{yw}=(h+e^2x)(-w/4+Be^2z),&P_{zw}=-2R_e-w^2/2+e^2z^2/2.
\end{array}
\]
全 \(e,B\) で多項式であり、\(e=0\) は元 \(U(B)\)。
\(e\ne0\) では真の T17 source の actual PGL 像なので全 Jacobi は零であり、
多項式性から中心を含め全 \(e\) で零。
actual affine Pfaffian も
\(-R_e(R_e+w^2/4-e^2z^2/4)\) で中心は非零。
従って閉性から **whole \(U(B)\)、全 \(B\) は T17** に入る。
これは一点の parent membership を higher source へ延長したものではない。

## 4. Rank 一・active v：一つの whole 族、親は未認識

\(n_{17}=n_{18}=0\)、\(F\) の \(u\) 係数零と conformal 式から
\(n_{11}=1/4,n_{19}=q=-1/4,n_{15}=0,n_{16}=s\)。
\(s=0\) へ中心化すると元 core は \(n_{12}=0\)、\(n_2,n_8\) 任意。
平方完成した全 \(F\) は
\[
 F=-[v+8(n_6+n_7w+n_8w^2)]^2/16
       +(4n_8^2-n_2/4)w^4.
\]
相異性は \(c=n_2-16n_8^2\ne0\) と同値。
actual universal action の \(a=-8n_8,r=1\) は \(n_8'=0,n_2'=c\) とする。
これは \(v'=v-aw^2\) と full twenty-parameter action による実変換である。
続く actual \(z'=z+8n_6+8n_7w\) が \(n_6,n_7\) を消す。
全残 Jacobi は \(n_3=0,n_4=n_0,n_5=n_1\)。
actual \(x'=x+4n_0\) の後、自由 \(B=n_1\)、\(c\ne0\) の whole 族は
\[
\begin{array}{ll}
 P_{xy}=Bw+cw^2+R/4,&P_{xz}=wx/4-3yz/8+Bw^2,\\
 P_{xw}=(z-3wy)/8,&P_{yz}=wy/4-z/2,\\
 P_{yw}=-w/4,&P_{zw}=-w^2/4.
\end{array}
\]
\(F=-(v^2+4cw^4)/16\) は二つの相異 rank 四 quadric の積。
\(c=0\) の平方補集合は対象外として保持する。

actual \((x',y',z',w')=(t^2x,ty,trz,rw)\)、tensor scalar \(1/t\) では
\(c'=t^2c/r^2,B'=t^2B/r\)。
\(B=0\) は \(t=1,r^2=c\) で \((c',B')=(1,0)\)、
\(B\ne0\) は \(t^2=c/B^2,r=c/B\) で \((1,1)\) となる。
全逆は非零尺度を割るだけである。
\(c=1\) で \(B\) を動かす同じ実多項式族により、零軌道は非零軌道の閉包に入る。
従ってこの全枝の E に必要な未認識代表は \(c=B=1\) 一点の実軌道だけである。
それがどの既知親に入るかは本稿では未決。

## 5. Rank 二：全三 core と、各一 parameter の full actual source

conformal trace identity より \(q\det H=0\) なので \(q=0\)。
元 core 行19は
\(n_{17}(2n_{11}-1)=n_{18}^2\)。従って
\[
 \det H=n_{17}(n_{11}-2n_{19}).
\]
rank 二では \(n_{17}\ne0\)。既 actual action で \(n_{17}=1,n_{18}=0\) へ送り、
元 core から順に
\[
 n_{11}=1/2,\quad n_{12}=-3\xi,\quad n_2=-\eta/4,
 \qquad \eta=n_{19},\ \xi=n_8,
\]
\[
 (1-12\eta)\xi=0,\qquad 48\xi^2+4\eta^2+\eta=0.
\]
従って全許容 core は次の三つ。
\[
 (\eta,\xi)=(0,0),\quad(-1/4,0),\quad
 (1/12,\ \xi^2=-1/432).
\]
最後の二符号は actual \(r=-1\) で交換される。\(\eta=1/4\) は rank 二条件に反し、表に無い。
全非零／零補集合をこの二元方程式だけで解いたので、generic core のみを選んでいない。

二 factor の fiber row が独立なので、その constant／linear-\(w\) 係数を
actual \(x,z\) への affine-\(w\) 加算で同時に消せる。
係数行列はその二 row の可逆 \(2\times2\) matrix で、逆による四つの係数を使う。
この変更は \(Z\) を保ち、全 source tensor も同時に pushforward する。
残る固定 factor 二次部分は core のみで決まり、\(\kappa_0=\eta-1/4\ne0\)、
\[
 U=u-\eta w^2/2,\qquad W=v-(\xi/\kappa_0)w^2
\]
として
\[
 F=-U^2+(\kappa_0/2)W^2.
\]
各 factor の \(V\)-eigenvalue は \(-1/2\)。よって \(s=Vw=n_{13}\) を保った全 \(V\) は一意に
\[
\begin{aligned}
 Vu&=-u/2+\eta w^2/4+\eta s w,\\
 Vv&=-v/2+(\xi/(2\kappa_0))w^2+(2\xi/\kappa_0)sw,\\
 Vw&=s
\end{aligned}
\]
となる。表の全 core で \(\xi/(2\kappa_0)=-3\xi\)。
\(V\) と \(\Phi=F\) が元全 tensor を一意に定め、未記録の磁気項は残らない。
具体的な元二十係数の完全 substitution は

```
n0=0, n1=−eta*s, n2=−eta/4, n3=0, n4=0, n5=−eta*s,
n6=n7=0, n8=xi, n9=0, n10=−12xi*s, n11=1/2,
n12=−3xi, n13=s, n14=n15=0, n16=s, n17=1, n18=0, n19=eta.
```

全元 Jacobi は上の二元 core 方程式とその \(s\) 倍へ直接還元される。
全零 \(s\) を含み、各表の whole source は実 \(\mathbb A^1_s\) である。
二 factor は affine \(x\) について monic な degree 二式で、斉次化後も \(h\) で割れないから既約。
\(\kappa_0\ne0\) と独立 fiber row より相異なる。

actual \((x',y',z',w')=(t^2x,ty,t^2z,tw)\)、tensor scalar \(1/t\) は
各固定 core を保ち \(s'=ts\) とする。
従って \(s\ne0\) は \(t=1/s\) で \(s'=1\) へ入り、\(s=0\) は同じ全実 polynomial source の閉境界。
各 whole source の known-parent 認識に必要なのは \(s=1\) の一軌道だけである。
現時点では三つとも親を未認識として残す。

## 6. 全支持への結論と残る有限義務

元 full twenty source の HQQ は、rank 零 Hessian 不可、rank 一の \(q=-1/2\) 不可、
rank 一 \(q=-1/4\) の二族、rank 二 \(q=0\) の三族で尽くされる。
各正常化は actual affine action と明示逆を持ち、反転した元を零とする補集合は各節で処理済み。
補助 Ga 座標の非線形自己同型を PGL 変更として使っていない。

rank 一 active-u 全族は §3 の真 T17 実曲線で閉じた。
残る E は rank 一 active-v の \((c,B)=(1,1)\) と、rank 二三 core の \(s=1\) の
**四つの actual 軌道**の named parent／全成分性で十分である。
各零 parameter はその同じ全 source の閉境界なので、正の帰属が得られた後に閉性で戻せる。
本稿はこの四つを既知成分へ入れたとも、新成分を発見したとも主張しない。

入力は旧全二十 potential／universal action の独立受理と元二十六 Jacobi 係数であり、
本稿の新部分は上の有限手還元・全 source 再構成・一つの actual T17 曲線である。
新独立 receipt は別途必要。旧 G 上界、既 H＋既約quartic 七源又は別の因子型の theorem を
本結論の代わりに使わない。
