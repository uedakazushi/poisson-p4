# 非被約 R・generic triple normal cubic：正しい curvature による全接続 V2

2026-11-08、producer 純証明。旧自由 r 分類の sign defect を
`PLANE_M3_FREE_CONNECTION_SIGN_ERRATUM.md` で明示撤回した後の独立した導出である。
旧九十源／十四 Ore 源の bytes とその入力に対する算術を変更しない。
本 note は全 first connection の有限被覆を証明する。
global K、全 higher source、G/E 閉鎖は別の義務であり、この note だけでは主張しない。

## 1. convention と正しい全八係数

rank 二 plane の actual reduced irreducible quintic が位数三で、normal leading cubic が
generic triple root とする。UFD と degree から primitive normal root は constant、
`G3=h z^3`、h は base quadratic で支持は R に含まれる。
normal rows の matrix を `β=[[p,0],[r,q]]` とし、その作用を
`Σ βij n_j ∂ni` とする。この linear-field 作用は matrix commutator の
anti-representation なので curvature は `dβ−β∧β` である。従って

`q=2p−dlog(f/h), dp=dq=0, dr−(q−p)∧r=0`。

normal Borel `w'=t w+s z`（t 非零）は
`r'=t r+s(p−q)`。inverse は `w=(w'−s z)/t` であり、
r の非零 scale と Borel image を実際の constant normal GL2 で戻せる。

別の検算：tracefree base field Z と
`r=−i_E i_Z vol/(3f)` を対応させると、`d i_E i_Z vol=3 i_Z vol` より
flatness は `(p−q−dlog f)(Z)=0`。
triple の q を代入すれば正確に **`(p+dlog h)(Z)=0`** となる。
この式は対角 p の正則特異性又は log 性を仮定しない。

R=2L+M では `f=a^2 b, h=a^e b^(2−e), e=0,1,2`、
`x=b/a,y=c/a` と置く。全 closed p は

`p=dlog f/3+α dlog(a/b)+γ dx+δ dy`。

これは全非 log 部分を保持した式である。
actual base change `c'=s c+t b+u a`、a,b の非零 scales により
δ 非零は `(γ,δ)=(0,1)`、δ 零／γ 非零は `(1,0)`、両零は別に取る。
δ 非零源は既存の全 r を残した universal obstruction
`PLANE_NONREDUCED_M3_TRIPLE_IRREGULAR_C_EXCLUSION.md` により全空。
以下はその排除に自由 r の flatness を使わない。

`η=4/3−e−α` とすると `q−p=η dlog x+γ dx+δ dy`。
全 global source `H0(Ω1(3))` は次の八係数である：

```
r=(P dx+Q dy)/x,
P=A+B x+C y+D x y+E y^2,
Q=F+G x+H y−D x^2−E x y.
```

従って正しい flatness の全 polynomial identity は

`x(Q_x−P_y)−(η+1)Q−γ x Q+δ x P=0`。

これが旧 plus curvature を用いた identity の代替である。
以下の各 kernel はこの八係数の全比較であり、選んだ r の部分空間ではない。

## 2. γ 一：二 whole line と η=−1 の一追加型

δ 零／γ 一では degree 三の係数から D=E 零、次いで G=H 零、
C=−F と `(η+1)F=0` を得る。
従って η≠−1 では `r=A dlog x+B dx`。
Borel image `η dlog x+dx` で B を消し、残 `A−η B` の零／非零を
normal ratio で正規化すると `r=0,dlog x` の二型である。
この二型は η=−1 にも存在するので、全 α line を保持してよい。

η=−1 の全 source は

`r=A dlog x+B dx+F(−y dx+dy)/x`。

Borel で B を消した後、F 非零なら y translation で残 dlog x を消し、
normal ratio と y scale で F 一へ戻せる。
従って追加は `r=(−y dx+dy)/x` の一型だけ。
これは `α=7/3−e` であり、旧 `α=1/3−e` ではない。
全 e で六 whole line と三 fixed となる。

## 3. log-only：generic 三型と四つの全特殊 fiber

γ=δ=0 の全係数条件は

```
C=−η G,
(η+1)F=(η+1)H=0,
(η−2)D=(η−2)E=0.
```

η≠−1,2 の全 source は
`r=A dlog x+B dx+G(−η y dlog x+dy)`。
η 非零なら Borel で A を消せる。
actual base stabilizer `x'=λ x, y'=s y+t x+u` の作用では
G 項から A に `η G u/s`、B に `(η−1)G t/s` が加わる。
従って η≠0,1 かつ G 非零なら B を消せる。
normal ratio と λ,s の非零 scales により generic 全三型は

`0, dx, −η y dlog x+dy`。

この三 whole α line は η∉{−1,0,1,2} の範囲で使い、
除いた全 fiber は次の有限正常形へ戻す。

* η=1：A を Borel で消すが B/G の両零非零を残し、
  `0, dx, −y dlog x+dy, −y dlog x+dy+dx` の四型。
* η=0：Borel image は零。G 非零なら y shear で B を消す。
  `0, dx, dlog x, dlog x+dx, dy, dlog x+dy` の六型。
* η=−1：`d(xr)=0` なので `r=dψ/x`。
  global ψ は全 degree 二以下 polynomial、constant は微分で消え、
  Borel は ψ の x 項を変える。actual affine triangular group と normal ratio で
  ψ は `y²,y²+x²,xy,xy+y,0,x²,y,x²+y` の八型。
  y² 非零では y shear/translation で xy,y を消し残 x² の零非零を取る。
  y² 零／xy 非零では x² を消し y の零非零を取る。
  両者零では x² と y の各零非零の四型。
* η=2：`d(r/x²)=0` なので `r=x²d(ψ/x²)`。
  global ψ は degree 二以下、x² は写像の kernel、constant は Borel image。
  全 ψ は `y²,y²+x,xy,xy+y,y,0,x` の七型。
  y² 非零では xy,y を消し x の零非零、y² 零／xy 非零では
  x を消し y の零非零、残 linear source は y,0,x へ入る。

ここで全 group element は λ,s,normal ratio 非零の actual base/normal change で、
その inverse も同じ triangular affine group にある。
η と α の対応は常に `α=4/3−e−η`。
η の零、normalizer resonance 一、kernel jumps −一／二を同一視しない。

## 4. canonical necessary equation による十二 fixed の全排除

前節の全八 r source に
`PLANE_M3_CANONICAL_TRANSVERSE_OFFDIAGONAL_OBSTRUCTION.md` を適用する。
log-only と γ 一では全 p/α で Cbad は非零 c 方向、κ 零となり、
canonical strict uv² equation は E=H 零を要求する。
これは free-r sign のどちらにも依存しない全 higher-source 必要条件である。

η=−1 の `ψ=y²,y²+x²` は H 非零、
η=2 の `ψ=y²,y²+x` は E 非零なので、各 e で四型を全空として排除する。
従って log-only の一 e あたり特殊 fiber は

```
η=1: 4型、η=0: 6型、η=−1: 6型、η=2: 5型
```

の計二十一 fixed。
排除後に r を勝手に gauge 固定したのではなく、全元 source への必要式である。

## 5. R=3L の全 offdiagonal source

`f=a³,h=a²` では `p=dlog a+dφ,q=dlog a+2dφ`。
φ は degree 二以下 polynomial modulo constant、actual affine base GL2 で
`0,x,x²,x²+y,xy` の五型。
全 r は

```
r=P dx+Q dy,
P=A+B x+C y+D x y+E y²,
Q=F+G x+H y−D x²−E x y.
```

flatness は **`Q_x−P_y−φ_x Q+φ_y P=0`**。

* φ 零：全 r=dψ。全 ψ は 0,x,x²,x²+y,xy の五型。
* φ=x：D=E=G=H 零、F=−C、
  `r=(A+B x+C y)dx−Cdy`。
  Borel で A を消す。C 非零なら y'=s y+t x+u で
  B は `B−Ct/s` へ動くので B を消し、全三型は
  `0,x dx,y dx−dy`。C 零の B 非零は normal ratio で一へ。
* φ=x²：全 `r=(A+B x)dx`。Borel で B を消し 0,dx の二型。
* φ=x²+y と xy：全 r は constant multiple dφ、Borel で零へ。

最後の二 potential は q(Cbad) の constant-vector contraction が injective なので
全 higher source 空。φ 零／ψ=xy は canonical transverse lemma により全空。
残 source は φ 零四型、φ=x 三型、φ=x² 二型、計九 fixed。

## 6. 全 source dimension と正確な有限計算対象

以上から全 generic triple first symbols は、排除済み全源を除き、

* R2L γ 一：六 whole α line＋三 fixed。
* R2L log-only：九 whole α line＋六十三 fixed。
* R3L：九 fixed。

従って **十五 whole α source＋七十五 fixed＝九十 families** で覆われる。
これは正しい curvature での全 r8／actual normalizer 分類であり、旧九十入力とは異なる。

固定 plane の source dimension は gauge 消去前で数える。
R2L では R の ordered lines 四、πS amplitude 一、constant root 一の六に、
p と r の元 solution stratum を加える。
γ 一 generic は二＋二、η−一 fixed も一＋三なので B≤10。
log generic は一＋三で B≤10、η零／一 fixed は B≤9、
η−一／二 fixed は r 五次元を保持して B≤11。
R3L の source 上界は φ 零四型で B=4,6,7,8、
φ=x 三型で B=7,8,9、φ=x² 二型で B=8,9。

全 normalizer の image だけを数えて消去した fiber を落とさない。
この上界から full thirty/nine higher kernel と合わせて
`B+k1+k2+5≤24` なら G 非優越が従う。
full D3 image 零又は full bad normal rows 零による actual divisor の z 因子も
さらに弱い十分条件である。

次の V2 計算は旧 literal source と β/P0 が同一な場合だけその exact 証明書を再利用する。
whole-line の旧除外 fiber は新 corrected table 又は original matrix の exact specialization で戻す。
旧 source と同一でないものだけ一 bounded batch に含める。
算術 partial PASS を全 G/E 閉鎖又は旧 coverage の復権と呼ばない。
