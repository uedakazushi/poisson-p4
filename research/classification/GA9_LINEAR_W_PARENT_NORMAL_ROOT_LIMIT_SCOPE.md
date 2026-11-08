# Ga9親認識：m二からm四へのnormal-root limitで保持すべき補集合

2026-11-08、純数学。T36／T37への所属・非所属を主張しない。新しいoperator、point、profile、jet又はCAS無し。

`GA9_LINEAR_W_FIXED_DIVISOR_MAGNETIC_NORMAL_FORM.md` は、正則なactual endpointに対する実Dとunique singular plane上のrank二からGa9全軌道を認識する定理である。本稿は、そのendpointへ至る全sourceを根拠なく狭めないための必要な区別を与える。

## 1. T36のproper planeはspecialでも一意にalign出来る

completed T36の一般源ではSε={u=v=0}がplaneで、実Dのnormal二次項は `(3p²/25)ABC uv`。pABC≠0でnormal二次式はnodeとなり、plane上の一般Poisson階数は二である。normal三次以上の全periodを捨ててこの事実を得たのではない。

Ga9 endpointの実Dは−4h(4h³w−8h²y²+y⁴)で、reduced singular supportは唯一のplane S={h=y=0}。proper Grassmannian上のincidence `Sε⊂Sing Dε` を閉包すると、special planeもSに一致する。DVR上のtautological plane bundleとそのquotientはfreeなので、unit determinantのregular GL5変更によりこのplane familyをSへalign出来る。

plane上のzero cubicもproper factor incidenceで扱える。T36の三linear factors A,B,Cが与えるtriangleからGa9のw³へspecializeする場合、各factorのprojective limitはwである。generic factorの一時的なscalar発散は、各lineをprojectiveにmarkすることで除ける。

## 2. generic nodeの二rootはtargetのquartic rootへ自動的には送れない

上のproper plane alignmentだけでは、normal二次式uvの二rootをtargetのnormal四次式h⁴のroot hへ送ったことにはならない。

S×Δをblow upすると、generic DはI_S²に属するので、全familyではexceptional factorを二回だけ除く。special DがI_S⁴に入ると、同じ全familyのstrict-expressionのspecial fiberには、なおexceptional divisorの二重成分が残る。

したがってgeneric node rootsをmarkしたproper exceptional flagsのlimitは、このvertical成分へ入ることが出来る。special Dそのものを四回割って得るstrict transformとの交点へ、同じflagsを移したと仮定してはならない。この移行にはno-vertical等の別の仮定又は証明が必要である。

## 3. Poisson条件を保つ明示反例

C⁴の座標(x,z,h,y)、tをparameterとし、任意の独立linear forms ℓ1(h,y),ℓ2(h,y)について

```
π_t=∂x∧∂z+f_t(h,y)∂h∧∂y,
f_t=t ℓ1ℓ2+h⁴+hy⁴
```

と置く。二つの直和blockは別の座標だけに依存するので、全tでstrict Jacobiは零。全familyでS={h=y=0}上の制限は∂x∧∂z、従ってrank二である。

t≠0でSの一般normal二次式はtℓ1ℓ2となり、二つの異なるrootは任意の固定ℓ1,ℓ2である。一方、t=0のnormal四次式はh⁴で、選んだ二rootがhと無関係であってもよい。例えばℓ1=y,ℓ2=h+yが使える。

y=r,h=rsのblowup chartでは

```
f_t/r²=t ℓ1(s,1)ℓ2(s,1)+r² s(s³+r).
```

special fiberのr²は正確に上のvertical exceptional成分である。generic二rootのlimitを、special strict transformのroot s=0へ送る議論はこの式で直接破れる。

これはglobal degree二P4 sourceの新しい例ではなく、**rank二Poisson planeとm二→m四だけからnormal flagsのalignmentを導く一般論の反例**である。whole T36への所属又は非所属へ拡張しない。

## 4. 残す正確な親認識義務

T36をS固定の全sourceで扱う場合、base三rootがwへcollideする全相対orderと、normal二rootがdistinctに残る場合／coalesceする場合、normal frameのdeterminantのvaluation、全normal-to-base shearを保持する。generic normal coordinatesの線形独立性をspecialにも仮定しない。

Ga9の固定D定理を用いた成功証明には、正則で非零なactual tensor curve、そのendpointのsigned Pfaffianが所望のhF、unique S上rank二、の三点で足りる。ただしD-only limitからtensorの正則性を補ったり、normal-root flagsを未証明で固定したりしない。

現段階でこのwhole actual curve又はwhole parentを排除するclosed invariantは未構成。Ga9 E所属は未解決のまま、既closed sourceや成分総数を変更しない。
