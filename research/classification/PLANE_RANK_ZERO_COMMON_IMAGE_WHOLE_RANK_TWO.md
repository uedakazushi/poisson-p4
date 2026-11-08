# rank零plane・common-image全四sourceは射影rank二以下

2026-11-08、PRIMARY E、担当 /root/log_manuscript。既受理の完全kernel・actual shear・全normalizer逆からの純証明。
新operator・rank・minor・点・GB・jet・数値計算は無い。main／frontier／globalと公開cutは変更しない。
別作者読解前の新producer証明として保存する。

**命題。** 全重複度abelian-conormal入口で、非零traceless symbolがglobal Vとconstant normal sl2行列の積、
trace TとVが線形独立・可換だがpointwise proportionalである全common-image sourceは射影rank二以下。
従って、同入口の適用対象であるrank四・被約actual退化因子の枝は空である。
旧G像二十上界を全点認識へ読み替えず、全tensorのrankを証明する。

## 1. 全入口と四normal forms

入力は PLANE_RANK_ZERO_COMMON_IMAGE_FIRST_JETS.md、PLANE_GLOBAL_COMMONIMAGE_SEPARATE_REVIEW.md、
不変producer CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMMONIMAGE_20261107.json、
不変review CLASSIFICATION_PLANE_GLOBAL_COMMONIMAGE_REVIEW_20261107.json。
全可換二plane分類とdet(Euler,T,V)=0によりcommon-image algebra Hom(C³/L,L)だけとなる。
image lineとordered independent pairの可逆GL2でT=b∂a,V=c∂aへ直す。pair determinant零は本仮定外。
全normal semisimple／nil、K零／非零の四型は

~~~
P0_ss = ∂a∧((b+c)z∂z+(b−c)w∂w)+k zw∂z∧∂w,
P0_nil = ∂a∧(bz∂z+(bw+cz)∂w)+k z²∂z∧∂w.
~~~

k非零の全actual逆は既受理である。ssは座標(a,b/k,c/k,z,w)＋tensor尺度1/k、
nilは(a,b,c/k,kz,w)＋尺度一でk一へ。逆はそれぞれ(a,kb,kc,z,w)、(a,b,kc,z/k,w)。
唯一の反転kの零補集合は後のK零二型で全て保持する。
全quadratic Euler商はP=P0+P1+P2、weight一raw36/Euler6の全三十方向と、
base pair ab/ac/bc×normal二次式z²,zw,w²の全九P2で尽くされる。
normalnil b非零affine型、dependent T=tV、traceless零、P0零は別入口であり本命題に含めない。

## 2. 実六shearの有限作用とP2

保存shear列の順序は**grouped normal order**であり、interleaved順ではない。

~~~
(s0,s1,s2,s3,s4,s5)=(a_z,b_z,c_z,a_w,b_w,c_w),
a'=a+s0*z+s3*w, b'=b+s1*z+s4*w, c'=c+s2*z+s5*w,
z'=z, w'=w.
~~~

全mapはdet一、逆は全sの符号反転。全quadratic pushのweight一Euler6を定数逆で除くと、
P1への作用は全sでP1→P1+Ss、Sは保存30×6 actual derivativeそのものである。
P0のweight一変化はsに一次。P0のs二次項、P1のs一次項、Euler補正の残りは全weight二へ入る。
全九P2を保持してこれらを吸収するので、inhomogeneous particularを捨てない。
これはnormal一次式をbaseへ加える有限三角代入からの厳密な次数論で、形式一次近似ではない。
P0は不変、逆shearで以下の全補方向へ実際に移せる。

## 3. SSの全非零K源

完全K1は六、actual Sもrank六で像はK1に入る。
元30座標rows (1,2,10,11,21,24)のS minorは

~~~
[[0,2,1,0,0,0], [0,1,0,0,0,0], [0,0,0,0,2,-1],
 [0,0,0,0,-1,0], [0,0,0,-1,0,-2], [1,0,2,0,0,0]]
~~~

でdet一。旧complete_kernel列順alpha0,…,alpha5に対し、全三十行の等号は

~~~
s=(alpha3−2alpha5,alpha4,alpha5,−alpha0+2alpha2,−alpha1,−alpha2),
P1=Ss.
~~~

逆shearでP1零とした後、weight二JacobiはD2(P2)=0。
旧全350×9 D2はrank九、保存minor−96なので全P2零。
実P1除去の後でforcingが零になるのであり、元P2をhomogeneous kernelへ狭めた議論ではない。
残P0_ssは三constant directions ∂a,∂z,∂wだけに属すから射影rank二以下。
SSについて旧producer BOUND §1が述べた短い空性の、有限作用と全P2の矢印を明確にした。

## 4. Nilの残る一方向も全rank二

旧K1は七次元、D2はrank九、保存minor−64。
完全kernel列順alpha0,…,alpha6は元free rows (18,19,20,21,24,25,26)そのもの。
同じrowsのSから、全三十行で厳密に

~~~
s=(alpha0+alpha2,alpha1,alpha2,alpha4+alpha6,alpha5,alpha6),
r=alpha3,   P1=Ss+rN,
N=∂a∧(−bw∂b−cw∂c+zw∂z).
~~~

Nは旧完全kernel第四列で、新chosen normalではない。非零raw座標はrow10=−1,row14=−1,row21=1。
source順はzL九、wM九、D三、B三、E三、C三であり、表示した三外積と正確に対応する。
全inverse shearでSsを消すとP=P0_nil+rN+P2。
Nはa非依存で共通∂aを持つためstrict [N,N]=0。旧D1N=0は元全350行で証明済み。
従ってweight二forcing零とD2単射性から、全P2零。rを反転せず、r零も含む。
残る全r tensorは

~~~
Pab=−r b w, Pac=−r c w, Paz=bz+rzw, Paw=bw+cz, Pzw=z²,
~~~

他の上三角entryは零。Eulを全radial場とすれば

~~~
P+Eul∧(−r w∂a)
 =∂a∧((bz+2rzw)∂z+(bw+cz+rw²)∂w)+z²∂z∧∂w.
~~~

右辺は三constant directionsだけなので全rでrank二以下。
Euler項はlinear fieldの外積であり、raw rankとprojective rankを混同しない。
同じ結論は実signed Euler Pfでも直接読める。omit bは−r c w z²、
omit cは−r b w z²、他三raw Pfは零だから

~~~
F=−b*(−r c w z²)+c*(−r b w z²)=0.
~~~

これは全rの全actual五次恒等式であり、一点の因子診断ではない。

## 5. K零と全scope

K零両型は旧独立全源証明で全十二P1＋全六P2を同時に戻して
P=∂a∧(W0+W1+W2)、全Wiはa非依存と完全に書け、全十八higherでrank二以下。
今回その証明を再実行しない。
これと§3–4、全k単位／零、全actual正常形の逆を合わせると、全四sourceにrank四点は無い。
新component増分零。dependent/scalar/zero-symbolなど別のrank零plane源、全E／総上界は残す。

本稿の新部分は旧完全S/K1の定数読解と有限weight作用、一つのdecomposable Nの手証明だけである。
全30行等号は別作者が既保存matrixからQQ比較できるが、新operatorを作る必要は無い。
旧全350・rank・minor・kernelを新再計算したとはしない。
全nonzeroKの旧像上界二十、旧all_passed=false、旧source bytesは不変である。
この純証明を別作者が受理した後にだけ新whole-E結論として採録する。
