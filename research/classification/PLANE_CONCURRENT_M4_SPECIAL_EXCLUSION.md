# concurrent cubic plane・m4 squarefree special rank 二 pencil の全空性

2026-11-07。主義務 G。本ノートの結論は、以下の仮定を満たす rank 四 Poisson tensor が存在しないことである。

- restriction rank 二の plane S の zero cubic は三相異 concurrent lines。
- 実退化 quintic D の S 沿い位数は四、最初の normal binary quartic は generic squarefree。
- I≡0 又は J≡0、leading coefficient map の image rank は二。

generic j と rank 一は既受理の essential-three/Saito 矛盾により別途空である。
これらと本ノートを合わせると、concurrent reduced cubic の **全 m4 squarefree rank 四源が空**になる。
重根、m2/m3、非被約 R はこの結論に含めない。
全点 E 帰属、新成分又は他の支持の空性も結論しない。

## 1. 全四 source への入口

`PLANE_CONCURRENT_M4_SPECIAL_PENCILS.md` の全 normal-form 論証を使う。
基礎となる全 W 分類は既受理の `PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_PENCILS.md` §§2–3 であり、
I0 の二型、J0 の common-kernel/secant と compression の両方を保つ。
三相異 I0 pencil を捨てず、前三型で未使用の第三線も許す。
full Hamiltonian parallel 条件により接続は一意なので、log や diagonal を仮定しない。

f=xy(x−y)、normal u,v、残る base 座標 z とする。
raw P0=∂z∧W0 における W0 の normal linear part M(x,y)(u,v)^t は次の通り。

| 型 | M |
|---|---|
| J0 secant | diag((−x+5y)/8,(5x−y)/8) |
| J0 compression | diag(−x/2+y,x−y/2) |
| I0 mixed | diag(y/2,x−y/2) |
| I0 三相異 | [[−x/2+y/4,y/2],[3y/8,x/2−y/4]] |

W0 の base part は fy∂x−fx∂y。全 M は z に依存しない。
各 displayed leading generators G について

\[
 (f_y\partial_x-f_x\partial_y)G+(Mn)\partial_nG-\operatorname{tr}(M)G=0
\]

を全係数で確認し、八 unknown の complete linear system が rank 八であることを保存した。
M に z 係数を初めから自由に入れても、z 係数は stable quartic の twisted infinitesimal stabilizer に入り零となる。
したがってこの表示に未処理の z-dependent 接続を捨てる仮定は無い。
全十五 leading coefficients の完全核は順に 2,2,2,1 次元で、全 normal form generators に等しい。

actual D2=0 は raw {u,v}=0 を一意に強制する。
一方 ∂z∧W0 の全350 Jacobi と actual Pfaffian はともに恒等零。
この段階の rank 二は全高次 tensor の階数の仮定ではなく、normal filtration の weight 零部分の構造である。

## 2. 一回の full30/9 証明書と粗上界の失敗

有限契約 `jobs/plane_concurrent_m4_special.json` の gate 後、reusable runner の
`--concurrent` を **一回だけ**実行した。
report は `reports/CLASSIFICATION_PLANE_CONCURRENT_M4_SPECIAL_20261107.json`。
3.642640秒、全 Fraction 算術 PASS。
旧六 FunctionDef の AST は全一致で、isolated `concurrent_main` だけを追記した。
旧 producer/report と essential-three eight-source packet は一切変更していない。

全350 Jacobi weight 一と D3=0 全24行は、各型とも 374×30、rank18、完全核12であった。
全 weight 二 operator は rank3、P2 fiber6。
original-row left functional による全二次 compatibility は零なので、粗像上界は31,31,31,32に留まる。
report の `all_concurrent_special_rank2_G24_closed=false` は、この **粗次元ルートが失敗した**ことを正しく残す。
これを G 成功へ書き換えていない。

しかし同じ report の完全核は、はるかに弱く十分な「全源が共通外積因子を持つ」という構造を示す。
以下は既保存の matrix/kernel の読解であり、新しい GB、minor、jet 又は追加点の計算ではない。

## 3. 全 P1 が ∂z 因子を持つ

三十 source basis に戻した全十二 kernel vectors は、四型の全てで次の同じ空間を張る。

\[
\begin{aligned}
 &\{x,z\},\{y,z\}\text{ の各 coefficient に }xu,yu,xv,yv,\\
 &\{z,u\}\text{ に }uv,v^2,\qquad
   \{z,v\}\text{ に }u^2,uv.
\end{aligned}
\]

これ以外の bracket/coefficient は全て零である。
したがって任意の全 P1 は ∂z∧W1 で、W1 は z-independent。
この空間の線形独立性と全374行への代入は report の complete kernel に保持されている。
rank18 の nonzero minors は順に

\[
 1166067873304425/274877906944,\quad3486784401/16384,\quad
 492075/8,\quad-129437595/262144.
\]

よってこれは chosen slice でなく、全必要 P1 源の等式である。
この共通因子から Jac(P1)=0 と P1²=0 が全パラメータで従う。
実際、report の全 original-row quadratic matrix も全係数零である。
したがって weight 二式 L(P2)+Q(P1)=0 の右辺は、全 P1 で厳密に零になる。

## 4. 全 P2 も同じ外積因子を持つ

全九 P2 basis は bracket pair (x,y),(x,z),(y,z) ごとに normal monomials u²,uv,v² の順である。
report の全 weight 二 operator L の非零 columns は **0,1,2だけ**であり、
残る3〜8列は全係数零である。
非零三列の rank は三であり、保存された原行 minor は次の通り。

| 型 | original rows | determinant |
|---|---|---:|
| J0 secant | 31,32,33 | 81/512 |
| J0 compression | 31,33,52 | −27/8 |
| I0 mixed | 32,33,51 | −3 |
| I0 三相異 | 31,32,33 | −3/8 |

§3で Q(P1)=0 だから、L(P2)=0 は {x,y} の三つの normal-quadratic coefficients を全て零にする。
従って任意の全 P2 も、{x,z},{y,z} の normal-quadratic coefficients 六個だけからなり、
P2=∂z∧W2、W2 は z-independent。
ここでも generic fiber を全点へ外挿せず、同じ固定線形 operator を全 P1 上で使っている。

## 5. rank 四との矛盾

全 global filtration と Euler quotient を保持した完全源について、

\[
 P=P_0+P_1+P_2=\partial_z\wedge(W_0+W_1+W_2).
\]

したがって raw tensor の rank は高々二、projective quotient での rank も高々二であり、
P²=0、実 Pfaffian D≡0 である。
これは rank 四の実非零 quintic が normal order 四の generic-squarefree leading term を持つ仮定に反する。
四つの全 source が同時に空となる。

group fiber11/11/11/10や粗次元上界をこの最後の結論には使わない。
正しい normalizer は report に保存されるが、不要な次元改善を追加していない。
全高 weight Jacobi をさらに課しても、この必要全源の rank≤2 を変えない。

## 6. 固定入力と独立確認待ち

| artifact | SHA256 |
|---|---|
| code/classification_plane_reducible_cubic_m4_special.py | 9dc3295c6b09837f70e753cb911efed20e2e5fd59203c54223971ec1ab30be4c |
| reports/CLASSIFICATION_PLANE_CONCURRENT_M4_SPECIAL_20261107.json | 44e89fd84700a144fcaa9110ee085a8d1b604e955d9c68752df6a12fec60d357 |
| research/classification/jobs/plane_concurrent_m4_special.json | 321407eaa0f1027fd6766709b1b91927cdc5cb44e5cc368831b019449b1b20d0 |
| research/classification/PLANE_CONCURRENT_M4_SPECIAL_PENCILS.md | ab1dabb57fd62e8041c7762a8b57f3b94099321e7e9eb3b2b710ab857b36829a |

producer による全算術と以上の数学証明が完了した。
別作者の全 source・connection/minor・共通因子証明の再構成は独立 receipt を要する。
