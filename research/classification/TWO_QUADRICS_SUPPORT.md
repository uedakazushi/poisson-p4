# 超平面＋相異なる二既約二次因子：全入口とpair-zero枝の被覆

2026-11-07、課題E_reduced_nonnormal。実退化因子 `D=h q₁ q₂` を扱う。
q₁,q₂は相異なる既約二次式であり、各rankは3,4,5を全て残す。
全因子型の被覆完了ではない。既存CPLUS_EXTの個別成分性を全支持へ外挿しない。

## 1. 全配置・全90係数を保持する逆形式入口

h-Casimir gaugeで、h=1上のPを次数二以下のstrict Poisson bivectorとする。
この操作とprimitive βの橋は `REDUCED_NONNORMAL_DIVISOR_SCOPE_AUDIT.md` 第3節に従う。
全退化式の条件は `Pf P=κ q₁q₂`, κ≠0。
したがって逆形式を

`ω=α/(q₁q₂)`, `deg α≤2`

と書くと、全Jacobi・Pfaffian条件は正確に

`(q₁q₂)dα−d(q₁q₂)∧α=0`, `Pf α=k q₁q₂`, `k≠0`。       (1)

固定した二次式の対について第一条件は全90係数への線形条件である。
第二条件を全て課すことで非退化逆形式を得る。
逆向きにも、(1)の任意の源は二次以下のbivectorへ戻り、
そのPfaffianはq₁q₂の非零定数倍である。
全hyperplane配置・q₁,q₂の全rank・定数/一次/二次係数を保持し、
`[C,Q]=0` や一次項零を仮定しない。

## 2. 二つの可換なaffine線形normal場

第一添字縮約 `P#df=i_df P` を用い、`Z=div P` とする。
`P#d(Pf P)=(Pf P)Z` と因子の互いに素性から

`Xᵢ=P#dqᵢ/qᵢ`

は多項式のaffine線形場である。斉次化すればhを固定する五変数線形場となる。
全支持について次が成立する。

`Z=X₁+X₂`, `Xᵢ(qᵢ)=0`, `[Xᵢ,P]=0`,

`X₁(q₂)=λq₂`, `X₂(q₁)=−λq₁`,

`tr X₁=λ=−tr X₂`, `[X₁,X₂]=0`。                        (2)

実際、log-Hamiltonian場の発散は `div(P#df)=−Z(f)`。
よって `tr X₁=−X₂(log q₁)` は定数であり、反対称性から他方のeigen関係も得る。
log-Hamiltonian同士の括弧は定数λのHamiltonian場なので零となる。
非退化開上で得た等式は多項式として全域へ延びる。
q₁,q₂の相異なる既約性からX₁,X₂は定数係数では独立である
（定数線形関係をPで戻し、各既約因子上のdlogのresidueを取れば係数が零）。

rank(qᵢ)=5なら、その非退化二次形式の保存から `tr Xᵢ=0`、従ってλ=0。
rank四でHがそのvertexを避ける場合も、radical直線上の固有値はh固定条件から零、
四次元非退化商上のtraceも零なのでλ=0となる。
したがってλ≠0の補集合は、両quadricが特異であり、
rank四のvertexはHに含まれる配置だけに限られる。
λ=0の場合には二つの独立可換場がh,q₁,q₂を全て保存する。
この制限は次のpair/pencil幾何の入口であり、まだ全配置の有限分類ではない。

## 3. generic pair-zero半単純modular枝の全支持

`HYPERPLANE_QUARTIC_COFACTOR_ENTRANCE.md` 第3節(a)の**完全な**許容weight源を使う。
この源の完全性は、modular半単純部の固有値が(a,−a,b,−b)で、
追加の定数・一次・二次共鳴を持たない範囲について、全weight比較で独立確認済みである。
以下の被覆自体は表示された全源のどの点にも適用する。

`u=xy`, `v=zw`, `s=−p−q−r`, `v₁=p+q`, `v₂=p+r` と置く。
全括弧は

`Pxy=C+Au+Bv`, `Pzw=D+Eu+Fv`,

`Pxz=p xz`, `Pxw=q xw`, `Pyz=r yz`, `Pyw=s yw`。             (3)

交換項B,E、全六logcanonical係数、全定数項を残す。
全四Jacobiの係数idealは正確に

`BE`, `v₂C+BD`, `EC−v₁D`, `B(v₂+F)`, `E(A−v₁)`。           (4)

全Pfaffianは

`(C+Au+Bv)(D+Eu+Fv)+v₁v₂uv`。                            (5)

以下の三stratumは(4)の全支持を覆う。いかなるminor零点も捨てていない。

### B=E=0

最高Pfaffianは `(AF+v₁v₂)xyzw`。
hが単純因子なので、この係数κは非零である。
重み(h,x,y,z,w)=(1,0,0,0,0)は全定数項を消し、
最高部は非退化logcanonicalであり、その実退化式は `κhxyzw`。
従ってこの全stratumはMPS40へ入る。

### B≠0

(4)から `E=0`, `D=−v₂C/B`, `F=−v₂`、さらに `v₁v₂C=0`。
C=0なら(5)はv=zwで割れ、相異なる二既約二次因子だけを持つ型ではない。
v₂=0なら全Pfaffianは零でrank四でない。
従って対象型ではC,v₂が共に非零、よってv₁=0。
このとき

`Pf P=−(v₂/B)(C+Au+Bv)(C+Bv)`。                         (6)

A=0なら両quadricは同一なので、相異なる対象型ではA≠0。
重み(h,x,y,z,w)=(3,0,0,1,1)を取る。
Bzw∂x∧∂yの重みは2、全定数項は正重み、全logcanonical項は零重みである。
極限のPfaffianは `−A v₂ xyzw` で非零。
従って全B主開の対象点が実SNC軌道極限を持つ。

### E≠0

対称の全式を直接用いると
`B=0`, `C=v₁D/E`, `A=v₁`, `v₁v₂D=0`。
対象型ではD,v₁,Fが全て非零、v₂=0であり、

`Pf P=(v₁/E)(D+Eu)(D+Eu+Fv)`。                         (7)

重み(h,x,y,z,w)=(3,1,1,0,0)で交換項Exy∂z∧∂wと全定数が消え、
極限Pfaffianは `v₁F xyzw` で非零となる。

各実PGL軌道極限は五独立超平面のSNC点である。
PGL5は連結なので、元の点を含む全空間成分はそのSNC点も含む。
MPS Theorem1.4により、全点は既知MPS40の閉和に入る。
特定のMPS番号や成分の一意性は仮定していない。
任意の因子rankやhyperplane配置を表示源内で追加除外せず、
既約・相異なるという対象条件から必要な非零条件だけを導いた。

## 4. λ=0、joint半単純rank二にも同じ全源が適用する

λ=0の二normal場X₁,X₂が生成する可換代数のjoint半単純部のrankを考える。
rank二の場合は、hを含む五weightのうちhのweightは零、全traceは零であり、
両quadricはweight零である。少なくとも一つのquadricがrank三以上なので、
五weightは必ず

`0,u,−u,v,−v`、u,vは独立                                   (8)

となる。実際、零weight空間の次元が二以上なら、残りは高々三座標である。
rank三の対称形式は少なくとも一つの反対weight対を必要とし、trace零から
残るweightも零になるため、全weightのspanは高々一次元となって矛盾する。
零weight空間が一次元なら、非零の反対weight対があり、残る二weightはtraceから反対になる。
もし反対対が無ければquadricのrankは一以下である。

(8)は全て相異なるので、両normal場の可換冪零部は零である。
従ってX₁,X₂自身が二次元の対角torusを張り、両不変性からPの全許容単項式は
`Pxy=C+Axy+Bzw`, `Pzw=D+Exy+Fzw` と四crossのlogcanonical項だけになる。
一次項は無く、二つの独立characterを使うため、modular場単独の追加共鳴に依存しない。
さらにZ=X₁+X₂はこのtorusに属するので、そのx,y固有値和が零であり、
正確に `p+q+r+s=0` が従う。これは(3)の全源そのものである。
従って第3節の全Jacobi・実SNC証明は、λ=0の**全joint半単純rank二枝**を覆う。
特定のmodular場が追加wallに乗っていても、この二normal場の条件で全源を保つ。
joint半単純rank一又は零はこの議論に含めない。

## 5. 有限証明書と未処理範囲

`code/classification_two_quadrics.py` は(3)から全四Jacobiを疎Fraction多項式で作り、
全係数idealと(4)の一致を確認する。
両主開で `B*inv=1` 又は `E*inv=1` の実際の逆元関係を用い、
全Jacobi、(6),(7)、全単項式の非負重み、全極限Pfaffianを再構成した。
CAS/producerをimportせず、独立cofactor reviewerのRing演算だけを再利用する。
gate後の実行は0.0036秒、証明書は
`CLASSIFICATION_TWO_QUADRICS_PAIR_ZERO_20261107.json`、
契約は `jobs/two_quadrics_pair_zero.json`。
これは表示源の全支持定理であり、第4節によりλ=0の全joint半単純rank二を閉じる。

後続 `TWO_QUADRICS_NIL22_BASE_SUPPORT.md` はJordan(2,2),c∈im Aの全支持を扱う。
元27係数源のglobal binary-discriminant Jacobi式と、非零平方枝の全閉源・Pfaffianにより、
被約1+2+2は全て平面七配置となる。その全源次元は高々24なのでGの該当枝も優越しない。
特殊Eの二平面源（tan-simple、chord）の帰属は残る。

`TWO_QUADRICS_NONZERO_PAIR_COFACTOR_SUPPORT.md` はλ≠0の全支持を閉じる。
全joint-character matchingを列挙し、rank(4,4),(4,3),(3,3)の三配置だけへ完全還元した。
各配置の全480×90条件はrank89の一parameter完全源であり、全源の実SNC軌道極限から
全点が既知MPS40の閉和へ入る。独立Fraction再構成は全係数と6100 matching組を確認した。

`TWO_QUADRICS_JOINT_SS_RANK1_SUPPORT.md` はλ=0・joint半単純rank一の全支持について、
追加nil方向を保持した五つの完全200×92supersourceを作る。
二型は線形因子又は既済nil22-imageへ帰着し、残る三型は必要source次元6/7/5と
同次元の明示normalizerの自由作用から全PGL像次元≤23となる。
全五型の全行・kernel・minor・元Jacobi/Pfaffian・必要strataは独立Fractionレビュー済み。
これはGの未知優越成分を除く弱い矢印であり、Eの全成分帰属は別である。

Gについて残るのはλ=0内のjoint半単純rank零（二normal場ともnilpotent）である。
Eについてはjoint半単純rank一の全成分認識も別に残る。
他の冪零modular場、線形場の零点でC=L=0となる純斉次源、
一般の二quadric pencilとHの全配置は未処理である。
特にquartic既約型の有限24+10表を、この型へ無条件で転用しない。
その有限表の低次wall必要性は、元のquarticの**既約性**を使っているからである。
一般1+2+2支持の被覆、Gの完了、成分数上界は主張しない。
