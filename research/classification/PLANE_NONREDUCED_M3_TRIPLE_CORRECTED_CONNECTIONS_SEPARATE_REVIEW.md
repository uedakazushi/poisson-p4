# 非被約 R・m三 generic triple：corrected V2 全接続の独立監査

2026-11-08、PRIMARY G。producer と別作者 `/root/quartic_sources` が
`PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_CONNECTIONS_V2.md`
（SHA256 `2fd5596aba1129b976e6336082379bfc19c34a615749a2767f9b76f27a65171b`）を
全八係数・全 actual normalizer と零補集合から独立に手監査した。
**十五 whole line＋七十五 fixed＝九十 family の必要被覆を受理する。**
旧 plus-sign 九十／十四 source の whole coverage は撤回されたまま。
その literal 入力での算術と、新しい corrected 全被覆を同一視しない。
本 receipt は純証明であり、新 CAS、operator、matrix、minor、jet 又は GB の実行は無い。
fullhigher/G/E の閉鎖は、この有限入口だけの結論ではない。

## 1. curvature と全 r八係数

primitive triple normal root は degree 零、G3=h z³。
β=[[p,0],[r,q]] の linear-field 作用は matrix bracket の反表現であり、
curvature は dβ−β∧β。
従って `q=2p−dlog(f/h), dr−(q−p)∧r=0`。
normal Borel w'=tw+sz は `r'=tr+s(p−q)`、t非零、逆は w=(w'−sz)/t。
別に r=−i_E i_Zvol/(3f)、tracefree Z の微分恒等式を使うと
flatness は `(p+dlogh)(Z)=0` となり、corrected sign を確認出来る。

R2Lでは f=a²b,h=a^e b^(2−e)、x=b/a,y=c/a、e=0,1,2。
全closed p は `dlogf/3+αdlog(a/b)+γdx+δdy`。
全global off-formは

`r=(Pdx+Qdy)/x`,
`P=A+Bx+Cy+Dxy+Ey²`,
`Q=F+Gx+Hy−Dx²−Exy`。

η=4/3−e−α と置き、全 flat identity を手で展開すると
`x(Qx−Py)−(η+1)Q−γxQ+δxP=0`。
三次／二次／一次／定数の各係数を全て比較し、root の V2 表と一致した。
δ非零の whole source は、全 r を固定せず original badvector 式で既に空となる。
その proof の q(C)=0 は Ca=Cc=0、α=5/6−e/2を与え、
元projective coefficientはκ=−7,α=4/3を必要とする。
e0/1/2では両者が一致しない。
この排除に自由rの誤sign分類は不要である。

## 2. R2L 全 kernel と全 resonance の actual 正常形

δ=0,γ=1では D=E=G=H=0,C=−F,(η+1)F=0。
η≠−1の全sourceは `A dlogx+Bdx`。
Borel image ηdlogx+dx でBを除き、残A−ηBの零／非零をnormal ratioで戻す。
全αを保持した二line r=0,dlogx となる。
η=−1ではさらに `F(−y dx+dy)/x`。
F非零は y translationで残dlogxを除け、y/normal scaleでFを一へ戻せる。
追加一fixedのαは **7/3−e**。
旧1/3−eを再利用していない。

log-onlyの全条件は
`C=−ηG, (η+1)F=(η+1)H=0, (η−2)D=(η−2)E=0`。
η∉{−1,0,1,2}では、BorelでA、y shearでBを消し、全三型
`0,dx,−ηy dlogx+dy` となる。
y'=sy+tx+u のG成分からAにはηGu/s、Bには(η−1)Gt/sの項が入り、
η0とη1の両失敗条件を区別する。
全λ,s,normal ratioの非零逆は同じ actual triangular affine group にある。

特殊 fibers は次の通り。

| η | canonical全型、canonical必要排除の前 |
|---|---|
| 1 | 0, dx, −y dlogx+dy, −y dlogx+dy+dx |
| 0 | 0, dx, dlogx, dlogx+dx, dy, dlogx+dy |
| −1 | r=dψ/x、ψ=y²,y²+x²,xy,xy+y,0,x²,y,x²+y |
| 2 | r=x²d(ψ/x²)、ψ=y²,y²+x,xy,xy+y,y,0,x |

η−1はd(xr)=0で、ψのx項だけがBorel image。
全二次係数のy²、xyの零／非零を順に分け、残linearをactual affine triangular群で戻せば表の八型を尽くす。
η2はd(r/x²)=0、ψのx²はkernel、constantはBorel image。
y²／xyの同じ二分とlinearの全零補集合で七型を尽くす。
一般 kernel を特殊値へ単に代入した分類ではない。

### badvector 方向の residue-zero 壁への補強

原note§4のlog-only全αでCbad∥∂cという箇所は、q(C)=0単独では一residue零の壁に不十分である。
以下の元projective必要式が、追加の局所化無しにその全域を補う。
log-onlyの全raw fieldsは
`X=3α a∂c, Y=(3e−4+6α)a∂c, κ=0`。
f=a²b の全 badvector 式を比較すると、Ca非零はXの係数−8,Yの係数−2を必要とし、
これらを代入するとe=6。
Cb非零はXの係数4,Yの係数1を必要とし、e=−1。
どちらもe0/1/2に反する。
従って全log-onlyでCa=Cb=0、既約性からCc非零。
γ1ではq(C)のb²／ab比較だけでCa=Cb=0が出る。
producerにもこの補強を共有し、frozen原noteは上書きしていない。

全r8のliteral curl inverseを展開すると ∂c(Za,Zb)=(3E,3H)。
canonical strict uv²式のa/b成分は(DZ C)ab=0なので、全r8で E=H=0を必要とする。
η−1のψ=y²,y²+x²、η2のψ=y²,y²+xが全higher源で空となる。
各eで四つ、計十二fixedを排除し、log-only特殊fiberは4+6+6+5=21固定源になる。

## 3. R3L 全 off-source と全 potential 正常形

f=a³,h=a² では p=dloga+dφ,q=dloga+2dφ、
全φはactual affine base GL2で `0,x,x²,x²+y,xy`。
全rは同じP,Q八係数だが、ここでは `r=Pdx+Qdy`。
flatnessは `Qx−Py−φx Q+φy P=0`。
各monomialを比較すると、

- φ0：r=dψ、全ψは0,x,x²,x²+y,xy。
- φ=x：D=E=G=H=0,F=−C、r=(A+Bx+Cy)dx−Cdy。
  BorelでA、C非零のy shearでBを除き、全三型0,x dx,y dx−dy。
- φ=x²：r=(A+Bx)dx、BorelでBを除き0,dx。
- φ=x²+y又はxy：rはconstant倍dφで、Borelで零へ。

φ0の全ψ分類はquadratic partのrank零／一／二と、rank一時のkernel-linear零／非零を保持する。
φ=xでは正しい非log型は **y dx−dy**。
旧 y dx+dy 型をcorrected whole源と扱わない。

最後の二φのq(C)は全constant-vector contractionがinjectiveなので全source空。
φ0,ψxyではq=dloga、C_a=0。
この時Z=diag(0,3,−3)、X=Y=K20=0。
canonical **全vector** uv²式は DZ C=0を要求し、Cb=Cc=0も強制する。
よってこのψxy whole sourceも空。
c方向だけにbadvectorを固定してから除いたのではない。
残りはφ0四型、φ=x三型、φ=x²二型の九fixedである。

## 4. 有限被覆数・安全な次元・正確な次工程

R2L γ1はe三通りの二wholeline＋一fixed、log-onlyはe三通りの三wholeline＋二十一fixed。
R3Lは九fixed。
合計十五wholeline＋七十五fixed＝九十family。
全marked lines、constant normal root direction、anchor amplitude、normalizer imageの全逆を戻す。
primitive rootのconstant性はdegreeから来ており、source parameterの一般点だけで被覆していない。

固定planeのsource次元はgauge前のmarked parameterで安全に数える。
R2L baseline六、γ1ではp/r全stratumを戻してB≤10。
log genericもB≤10、η0/1はB≤9、η−1/2はB≤11。
R3Lはφ0の四型でB4,6,7,8、φ=xの三型でB7,8,9、φ=x²の二型でB8,9。
これらはsource上界であり、canonical点のorbit次元を特殊fiberへ外挿した数ではない。

この純証明は全corrected firstconnection入口を受理するだけで、
九十の全K/P0/30/9とα壁、旧literal有効Oreへの全橋、十八affine modules、
新fixedの元Jacobi ideal因子は、それぞれ独立算術を要する。
旧九十reportのliteral入力を再利用する場合はrawP0と全targetのexact同一性だけを根拠とし、
撤回された旧whole coverageを依存に戻さない。
全G閉鎖又は全pointwise E帰属は、このreceipt単独では主張しない。
