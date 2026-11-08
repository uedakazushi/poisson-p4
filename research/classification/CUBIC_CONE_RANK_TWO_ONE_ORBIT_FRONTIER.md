# Cubic-cone rank二の残一軌道：真正ambientと一意二重零点

2026-11-08、PRIMARY E。全sourceの還元は
`CUBIC_CONE_RANK_TWO_VERTEX_SHEAR_REDUCTION.md` と
`CUBIC_CONE_RANK_TWO_CRITICAL_ORBITS.md` の凍結済み証明に従う。
残る全unit源は一つのactual PGL軌道O（二十四次元）へ戻り、nil全fiberも
その閉包Cへ入る。ここではCのcomponenthood又は新componentを主張しない。

## 1. 一点を選ぶ理由と真正全ambient

λ=3/5、保存六kernelのk0=1,k5=2、他零は、全unit源を実shear／Cartanで
一軌道へ戻した後のcanonical代表である。標本から全sourceへ推測していない。
座標を(x0,x1,x2,x3,c)とし、整数分母消去scaleは5。
全十括弧は `CLASSIFICATION_CUBIC_CONE_RANK_TWO_ONE_ORBIT_AMBIENT_20261108.json`
にliteralに保存した。全350projective Jacobi零、全126Pf、全150原係数を保持した。

exact QQの結果は次の通り。

* full350×150のDはrank100、全kernel次元50。
* Euler25はrank24、Euler+GL25+scalarの全51列はrank49。
* projective orbitは24、projective tangentは25、真正normalは一本。
* 全maximal minor、全kernel、free-coordinate identity、D*M=0を保存した。

active gate後1.30534秒、RSS51,656KiB。既helperのexact bytesも報告へ埋め込んだ。
したがってsmooth24という結論は得ていない。normal一本だけからcomponenthood、
親成分又は新成分を推測せず、追加point／rank profile／jet反復も行っていない。

## 2. 実quinticと全reduced singular locus

q0=x0*x2−x1²、q1=x0*x3−x1*x2、q2=x1*x3−x2²、
Cminus=x2*q2−x3*q1、Cplus=x0*q1−x1*q0とおく。
canonical代表のquinticは非零定数を除き

`F=c*q0*q2 + q2*Cminus + q0*Cplus`。

より一般に二係数a,b非零の
`F=c*q0*q2+a*q2*Cminus+b*q0*Cplus`について同じ議論が通る。
∂cF=q0*q2なのでSing(F)はq0*q2=0に含まれる。
q0=0,q2≠0ならF=0からCminus=0。

`x0*Cminus=−q1²−x2²*q0`

よりx0≠0の場合q1=0、さらにsyzygy
`x0*q2−x1*q1+x2*q0=0`からq2=0となり矛盾。
x0=0の場合はx1=0、q2=−x2²≠0、Cminus=−x2³≠0となり矛盾。
q2=0,q0≠0も、

`x3*Cplus=q1²+x1²*q2`

と対称なsyzygyで同様に矛盾する。従ってSing(F)はq0=q2=0に含まれる。
q1=0ならS=V(q0,q1,q2)であり、FはI(S)²内なので全SがSing(F)に入る。
q1≠0なら二つのsyzygyからx1=x2=0、x0*x3≠0。
このとき

`∂x1F=−a*x0*x3³, ∂x2F=b*x0³*x3`

は非零なのでsingularでない。以上から**Sing(F)_redは正確にtwisted-cubic cone S**。
特にsingular planeは無い。q1≠0のsecant planeをSに混同していない。

## 3. projective bivectorのorder≥2零点はc-vertexだけ

projective bivectorが二次以上で零になる点では、局所Pfaffianは四次以上で零になる。
従ってその点は上のSing(F)=S内である。
Sのvertex以外では、まずx0=1のchartを取り

`(x0,x1,x2,x3,c)=(1,t,t²,t³,eta)`

と書く。保存済み整数tensorのdehomogenizationから、S上の制限は正確に

`5*(1−t^5+t*eta) * (∂t wedge ∂eta)`。

t=0では零でない。t≠0の零ではeta微分が5t≠0なので一階で既に非零。
x0=0の非vertex点はx3≠0であり、x3=1、s=x2/x3=0のchartでは
`{s,c/x3}=5`なのでそこでも零でない。
Sはvertex外でsmoothなので、制限の一階微分が非零ならambient bivectorも
order≥2では零になれない。
一方、c=1のvertexでは全affine bracketの最小次数が二であり、実際order二零。
従ってこの点は**全P4で一意のorder≥2零点**である。

## 4. displayed T36 vertexとのexact jet比較

これは別のanchor計算でなく、二つのliteral sourceからの手による係数比較である。
completed T36は既存 `MPS_SMOOTHING_COMPLETION_20261103_v62/T36.json` を使う。
そのx0-vertexのaffine quadratic jetをQ_T、normal変数を(y1,y2,y3,y4)とする。
column divergenceは

`(−9p0*y1/5, +9p0*y2/5, +3p0*y3/5, −3p0*y4/5)`。

これはp0≠0でregular semisimpleであり、全t0,…,t3に依存しない。
我々のc-vertexのjetはPbase+Ebase wedge dで、column divergenceは

`(−15x0,−5x1,+5x2,+15x3)`。

scalarを除く固有値照合は二通りだけで、正の対応では
`(y1,y2,y3,y4)=(x0,x3,x2,x1)`、負の対応は逆順となる。
固有値が相異なるため可能なlinear eigenframeは各固有線のscaleのみ。
正対応では我々のjetは

`Q12=5*y1*y2+10*y3*y4, Q14=5*y1*y4`。

T36では

`Q_T12=p0*y1*y2+10*t2*t3*y3*y4/p0, Q_T14=(p0/5)*y1*y4`。

log monomial二係数の比は我々で1、T36で5。対角frameでもtensor scalarでも
この比は変わらない。負対応でも両係数が同時に符号反転するため同じ矛盾。
vertexを保つprojective mapの四vertex shearはquadratic jetを変えない。
従って**T36のdisplayed order二vertexを我々のvertexへ送るactual PGL同値は無い**。

§3により我々のorder二vertexは一意なので、我々の軌道はcompleted T36の
generic全nonzero sourceの軌道と異なる。既受理のT36二十四次元成分と我々の
二十四次元軌道閉包が等しいなら両者のdense open orbitは一致するはずであり、矛盾。
これはhigher-dimensional親からの非包含を述べるものではない。

## 5. 残る最弱義務

全sourceを一軌道閉包へ戻す証明とD100/M49は保存したが、Cが全Poisson schemeの
componentであることは未証明である。別のknown parentへの全actual curve、又は
真正一本normalの有限支持判定が必要。新component数は零のまま維持する。
他のconormal/source branchを、このnormal一本やT36非同値だけで消去しない。
