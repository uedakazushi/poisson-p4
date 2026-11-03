# 修正された条件付き quad78 補題

2026-11-03。root・a_zero が提案した補正版の数学的段階を検証し、universal polynomial identity により復元を確認する。元のG1/G2だけによる未受理の推論は復活させない。実際のgate locusに対する新しい算術探索は行わない。

## 仮定

各幾何学的characteristic-zero点で aPJ≠0 とする。六変数は `(A,X,W,B,S,U)`、c=6075、

`Q1=AB/c+W*C1`, `Q2=40XB/c+W*C2`, `D=40X*C1−A*C2`

とする。C1,C2はA,X,W,Bの線形式だけを含む。actual source/frame、binary Dの既知unit、全五式の**定数**独立性は先行の独立binary/whole-pair/CI-open監査に固定される。三minorのbase C⊂P3はHilbert–BurchによりCMの純curveである。

追加する条件は以下である。

- **G1:** Q0,Q3 mod(B,W) が互いに素。
- **G2:** ternary pair Q0,Q3 mod(A,X,W) が互いに素。
- **G3:** vertex上のS²,SU,U²の2×3係数行列がrank2。あるいはG3の代わりに **Cが幾何学的にintegral** と仮定する。

strong2×5 theorem はG3を意味しない。これらの条件が全開で成立するとは主張しない。

## 五式のlinear syzygy

五式Q0,Q1,Q2,Q3,Dのlinear syzygyを考える。A=X=W=0でscroll三式は零であり、G2によりQ0,Q3のlinear cofactor K0,K3はspan(A,X,W)に入る。

G3を仮定すれば、S²,SU,U²の係数比較にはscroll三式のlinear cofactorは寄与しない。rank2によりK0=K3=0となる。B*S/B*Uは使わない。

残る三minorのheightは全aPJ-openで2であり、Hilbert–Burchの最小分解はlinear syzygy二つを与える。六変数への多項式拡大でもこの分解は保たれる。従ってlinear-syzygy空間は正確に

`W*D−40X*Q1+A*Q2=0`,

`B*D+c*C2*Q1−c*C1*Q2=0`

の二つで生成される。

### 代替条件：integral C

G3を仮定しない場合でも、Cが幾何学的にintegralなら同じ結論を得る。height2のCMなhomogeneous座標環はdepth2で飽和され、Cのintegralityと合わせて base coordinate ring はdomainとなる。S,Uを加えてもdomainである。

vertex binary係数行列のrank0は、ternary pairの両式がBで割れることからG2に反する。rank2は前述の通り。rank1の場合、非零の一定kernel vector `(b0,b3)` により `(K0,K3)=(b0*K,b3*K)` と書ける（一方のbinaryが零の場合も含む）。scroll idealを法として

`K*(b0*Q0+b3*Q3)=0`

となる。非零Kはlinearであり、scroll idealはquadricsで生成されるためそのclassは零でない。domain性によりQ0,Q3の非零定数結合がscroll idealに入り、次数2ではQ1,Q2,Dの定数結合となる。これは受理済みの五式の定数独立性に反する。従ってK=0である。

なお、smooth Cという条件を使う場合にはconnectednessも必要である。このCM cubic curveはHilbert–Burch分解からH⁰(O_C)=k、従って幾何学的にconnectedである。幾何学的smoothnessと合わせてintegralとなる。特定のE5条件がactual Cのsmoothnessを保証する算術証明は本packetでは検査せず、別の入力を要する。

## 四式のquadratic-cofactor syzygyを六Koszulへ縮約する

四式Q0,Q1,Q2,Q3のquadratic cofactorをF0,F1,F2,F3とする。mod(B,W)ではG1により `(F0,F3)=lambda*(Q3,−Q0)`、lambdaは定数。対応するKoszulを引くと

`F0=B*L0+W*K0`, `F3=B*L3+W*K3`

と書ける。全関係をmodWでB除算し、A=X=0へ制限するとG2よりL0,L3 modWは(A,X)に入る。W部分をKへ吸収し、L_i=alpha_i*A+beta_i*Xと取る。

各iについてQ1,Q2との二Koszulを引けば、

`B*(alpha*A+beta*X)+W*K − c*alpha*Q1−(c*beta/40)*Q2`

`=W*(K−c*alpha*C1−c*beta*C2/40)`

なのでB部分を消せる。計四Koszulを引いた後はF0,F3がW倍となる。残りをmodWに取るとA*F1+40X*F2=0。A,Xは互いに素なので、linear L,K1,K2により

`F1=40X*L+W*K1`, `F2=−A*L+W*K2`

と書ける。元の全関係をWで割ると、五式に対するlinear syzygy

`K0*Q0+K1*Q1+K2*Q2+K3*Q3+L*D=0`

を得る。上記の追加条件により二つのscroll syzygyだけである。第一は元の四cofactorへ復元すると零、第二はc倍のQ1,Q2間のKoszulとなる。これをuniversal polynomial ringで直接確認する。

従って元のquadratic-cofactor syzygyは六つのKoszulだけである。四式は五式の定数独立性から独立なので、この六つのKoszulも定数上独立である。六変数二次cofactorの次元は4×21=84であり、raw84×126のrankは84−6=78となる。

これは **G1∧G2∧(G3またはintegral C)** の下でのdegree4の補題である。全CIを仮定するshortcutは用いない。G1/G2のみの適用、全開quad78、全四次全射、長65、ambient除去は未証明のままである。

再現（新しい出力先）:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_quad78_corrected_conditional_20261103_v63.py \
  --output work/marked-square-v63/quad78-conditional-new
```

codeは二scroll恒等式・四Koszul除去・W除算・両scrollの復元と六Koszulの独立anchorを確認する。G1/G2/G3・integralityの実locusは計算しない。
