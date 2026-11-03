# generic W枝のCA・CB・N全完成環：独立検証

`AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63` を、generic W枝に限る結論として受理する。CA/CBはAc8、NはAc7と滑らか因子の積になる。全ての特定W≠0へのtransportは、別のorbit定理を必要とし、このreviewには含めない。

## actual曲線とminorのbinding

三つのproper sourceの全W曲線は同じphysical tensor

\[
P_W=\partial_y\wedge(x^2\partial_h+hx\partial_z+(h^2+z^2)\partial_w)
     +W x^2\partial_z\wedge\partial_w
\]

であることを、別の150係数抽出で確認した。raw row140は厳密にW、row106は1であり、actual algebraic曲線をformal tangent直線で置き換えていない。CA/CBの全monic critical表示と全150 Schur式、Nの全critical式もこのW曲線上で恒等的に零となる。

既受理のactual七一次odd functionをf、全firstjetが零になる二residualをhとする。七つのraw base方向のminorはCAが−8、CBが1/128である。これを逆行列で正規化し、critical tangentのU,V方向を加えた九方向を独立に再構成した。実R座標での全W曲線の微分が保存T_Wに一致し、T_W,T_U,T_Vはactual一次odd行列のkernel内にある。

二residualの実二次項は、CAでは

\[
C_A=\begin{pmatrix}4/25&-4/25\\9/25&1/25\end{pmatrix},
\]

CBでは

\[
C_B=\begin{pmatrix}7/45&-1/5\\-1/15&-1/5\end{pmatrix}
\]

をWU,WVへ掛けたもので、determinantは8/125、−2/45である。全actual coefficientとの照合とstaircaseのdeterminantを再計算した。

高次項の扱いも確認した。actual W曲線上で全odd functionは零であり、正規化七方向に関する上段blockはI+O(W)、T_U,V方向はO(W)である。hの**full R firstjet**が零なので、下段左もO(W)、下段右はWC+O(W²)となる。Schurの交差補正はO(W²)であるから、九minorの先頭係数は厳密にdet(C)W²となる。raw九方向では−64W²/125、−W²/2880である。未計算の高次blockを零と仮定していない。

## 全proper fiberとrank2 source

独立受理済みのP0での**whole proper fiber**の長さ2を使う。proper morphismのquasi-finite locusで有限化し、有限coherent moduleのfiber lengthの上半連続性により、P0を含むtarget開で全fiber長は2以下となる。曲線はW=0でこの開を通るので、generic Wも含まれる。

全critical表示の単元fiber blockと八flag blockは維持され、W曲線上でsource微分はCA/CBではrank28、Nではrank27、いずれもcorank1である。選択source点のfiber tangentは非零なので、その局所fiber長は2以上である。全fiber長2以下と合わせ、generic Wの**全fiberは一台長2**となり、他のsource点を捨てていない。

また、P0で検証したactual parameter底への非退化二次係数はgeneric Wでも非零である。選択linear blockを陰関数消去し、残るparameter idealの長さ2から、正則sourceはactual regular底R上finite free rank2となる。必要ならétale近傍で該当するlength2成分を分離し、1,pを基底に取り、traceを引いて \(B=R[t]/(t^2-\nu)\) とする。

全target imageは \(A=R\oplus It\) である。全proper fiberが一台なので、actual targetの完成環がこのsource完成環へ注入することも保持される。actual Iが非零と分かった後に \(t=(ct)/c\) を用いてfraction fieldが同じことを確認し、Bを全normalizationと同定できる。

## CA/CBのNakayamaとalgebraic generic W

全critical quotientはWをinvertするとsmooth20次元枝になる。従って \(J=(I,\nu)\) はregular29次元底で正則九生成idealである。上のactual九minorはgeneric Wで非零だから、Iの元がJ/mJの九classを張る。有限module J/IへのNakayamaによりI=J、特にν∈Iとなる。

formal Wで非零なminorを、特定のW数値へ外挿してはいない。finite sourceとactual parameter projectionは有限型代数から来る。étale局所のfinite free rank2 algebraでtrace、odd係数ideal、Fitting minorsはcoherentな代数的対象として定義できる。Jのregular locusとactual九minor非零は開条件である。もし元のalgebraic W曲線全体がそのbad locusに含まれれば、W=0でのformal曲線へのpullbackも零となり、上の非零先頭係数に矛盾する。従ってalgebraic曲線のgeneric point、さらに非空開の一般closed pointで同じ条件が成立する。

これによりCA/CBの全cone完成環は

\[
\mathbf C[[g_{20},s_1,\ldots,s_8,t^2,s_1t,\ldots,s_8t,t^3]]
\subset\mathbf C[[g_{20},s_1,\ldots,s_8,t]]
\]

となる。dim/edim/depth/multiplicityは29/38/21/2である。scalar方向を除くprojective版はgが19個となり28/37/20/2となる。

## NのAc7型

Nでは別途独立受理済みの全ideal \(I=J\) を使う。formal normal formは \((\nu,s_6,wv,v^2)\) で、physical W曲線はそのcritical w枝を通り、generic Wでcritical quotientはsmooth20次元となる。normal-formのwを無条件に全raw target座標と同一視する必要はない。critical quotient上のwはsource Wに一致するようliftでき、少なくともphysical曲線上でWの単位倍となる。従ってgeneric WでJは正則八生成idealとなる。

I=Jはactual coherent idealの等号であり、formal completionでの等号はfaithful flatnessにより必要なlocal algebraへ降りる。generic Wの全normalization fiberの議論と合わせ、Nの全完成環は

\[
\mathbf C[[g_{20},s_1,\ldots,s_7,t^2,s_1t,\ldots,s_7t,t^3]]
\subset\mathbf C[[g_{20},s_1,\ldots,s_7,t]]
\]

となる。coneでは28/36/21/2、projectiveでは27/35/20/2である。八正則生成元に対する28 Koszul関係と36 pair productsの全64式を独立に再構成し、pullbackが零となることを確認した。rank2基底によるnormal formと正則列のKoszul syzygiesから、これらは全kernelである。

三つのnormalizationはいずれもfinite elementary subintegralである。t²,t³がactual imageに入るためであり、二つの異なるgeometric preimageを持つ型ではない。

generic Vでの数値Hessianやodd-order計算は、この独立証明に使用していない。B18原点のCA/CB全ideal、actual D、全global不良locus、任意の特定Wへの延長も結論しない。

独立codeは `code/ambient_B18_generic_W_review_20261103_v63.py`、receiptは `reports/AMBIENT_B18_GENERIC_W_FULL_COMPLETIONS_20261103_v63_REVIEW.json` である。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/ambient_B18_generic_W_review_20261103_v63.py \
  --replay work/generic-W-producer-replay \
  --output work/generic-W-independent-review.json
```
