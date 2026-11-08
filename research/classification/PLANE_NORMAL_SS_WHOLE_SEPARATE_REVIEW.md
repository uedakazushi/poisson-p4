# Normal SS二parameter三源とmixed三parameter：全source独立査読

2026-11-08、quartic_review、PRIMARY E。**両producer全源定理をPASS**。
新D/operator/M1/kernel/minor/profile/jet又はCASを実行していない。
旧whole QQ二源の全line／全根／全kernel／全sliceは既別作者reviewを固定入力とし、
今回のactual shear、元compatibilityの支持帰結、全FULL九P2、三parameterの普遍恒等式を独立に検査した。

## 二parameter三源の全shearと非剛性点

P0=T∧En+V∧(z∂z−w∂w)+kzw∂z∧∂w、T,V可換・trace零とする。
w=1,u=z/w,x=(a,b,c)/wではP0=Ex∧Ux+u(Bx)∧∂u、U=T−V,B=2V+kI。
actual shear x→x+v+u rのkernel条件はproducer (1)と一致した。
非零trace零Rに対するv∧Rx+x∧Rv=0を全xで課すと、
Rv=mu v、im(R−mu I)⊂Cv、traceとv上の固有値からmu=0、
Rはrank一nilpotentでv∈imRとなる。
これよりMIXHのshear rank低下は(t,k)=(±1,−2)、NILREGは(±1,0)、
MIXNは全域rank六となる。rank五点で六自由度を引いていない。

有限shearは元homogeneous base-from-normalの線形PGLで、逆は全六係数の符号反転。
weight一の全kernelがshear像なら全P1を零へ戻せる。
そのときFULL九P2はD2核に入り、D2単射ならP2=0。
有限shearで生じる二次補正も同時に移しているため、inhomogeneous P2を失わない。

旧全line／root証明書から戻した非剛性点は、MIXNの(3,3),(−3,−3)と、
MIXHの十二点(±1,4),(±3,−2),(±1,−8),(±1,1),(0,1),(0,−2),(±1,−2)。
producer metadataのrank／全点表と旧complete sourceを照合した。
前六MIXH点及びMIXN点のextra Rは全て共通constant方向を持つstrict Poisson tensor。
shear補空間としての独立性は、保存列／mixed anchor又は零anchorと可逆Bから確認できる。
D2単射で変換後P2=0。全一parameter残余のb又はc-rowがその座標で割れる。
四MIXHのD2核一点源は全kernelを保持しても同じb又はc因子を持つ。

負点への実normal swapとbase共役を検査した。
MIXNはN→−Nとの合成で(t,k)→(−t,−k)。
MIXHはさらにtensor scale −1を入れて(t,k)→(−t,k)。
NILREGはJ→−J,J²→J²との合成で(t,k)→(−t,−k)。
全てglobal可逆であり、正点の一部だけからの外挿ではない。

## MIXH criticalでFULL九を直接保持する

(t,k)=(1,−2)で、保存核八のshear像は
span(k0,k1,2k2+k3,k6,k7)、残三方向はk2,k4,k5。
保存列から、全残P1は
w(rb−2sa)∂a∧∂b+scw∂b∧∂c+ell zw∂a∧∂zと一致する。
FULL九P2=A2(z,w)∂a∧∂b+B2(z,w)∂a∧∂c+C2(z,w)∂b∧∂cを保持する。

元cyclic Jacobiをc=0で直接手展開した結果は

```
Jacz=z(−2B2+C2−2w∂wB2), Jacw=w(C2+2z∂zB2), Jczw=0,
Jbcz=−2zC2−2zw∂wC2, Jbcw=2zw∂zC2.
```

従ってdegree二のEuler恒等式から、元projective四成分は正確に
(E∧J)aczw|c0=6zwB2、(E∧J)bczw|c0=6zwC2となる。
strict Jを零と仮定したのではない。
polynomial域でB2=C2=0が必要、全r,s,ell,A2はそのまま保持してc不変。
負criticalも上の実逆で同じ因子結論となる。

## NILREGのK零whole line

k非零は旧全wall後の剛性に戻る。
k零ではgeneric二slice変数の両square compatibilityと全例外±3,±1,±1/3を保持する。
generic及び±1/3でsliceは原点となり、FULL P2はQ2(z,w)∂a∧∂bの全三核方向。
t=3の旧rows132/134=u0²/4、u0u2/2+3u1²/4からu0=u1=0。
残P1=u2w²∂a∧∂zはstrict Poisson。
t=1はshear五・slice四を保持し、rows33/132=u3²/2、u0²からu3=u0=0。
残P1=w[(u1c+u2b)∂a∧∂b+u2c∂a∧∂c]もa非依存の∂a-wedgeでstrict Poisson。
両方とも**残った全P1のself-Jacobiが零**なので、全inhomogeneous forcingが消え、
P2は全三核Q2∂a∧∂bへ戻る。先にhomogeneous核だけを仮定していない。
全tensorのc-rowがcで割れ、負の三例外も実逆で同じ結論を得る。

これでMIXH／MIXN／NILREGの全二parameter源は、全higher・零境界・shear jumpを含め、
actual quinticにproper hyperplane因子又はPf零を持つ。

## Mixed三parameterの普遍bad-vector条件

T=tH+sN、V=H+N、H=diag(1,1,−2)、N=b∂a、t≠sを全て保持する。
FULL P1=zL+wM+z²D∧∂w+zw(B∧∂z+E∧∂w)+w²C∧∂z、
FULL九P2はbase-pair binary quadraticである。
C=0なら全z-rowがz倍、D=0なら全w-rowがw倍なので、対象は両方非零を必要とする。

w=1,u=0において、真のprojective affine括弧のP0は
{xi,xj}=xi(Ux)j−xj(Ux)i、{xi,u}=u(Rx)i、U=T−V,R=2V+kI。
P1のmixed anchorはC、base bracketはM(x)でlinearである。
元Jac(xi,xj,u)の**linear-in-x**係数を全sourceから取ると

```
C∧Rx−∂C[x∧Ux] = C∧(−T+3V+kI)x−x∧(T−V)C = 0.
```

他のP1自己積はこのnormal rootでbase degree零、
未知base-linear MのP0との項はmixed P0のu因子で消える。
FULL P2はbaseのみでmixed anchorを持たず、このlinear係数へ寄与しない。
従って未知base P1又は非斉次P2を零と置いて得た式ではない。
別chartのnormal swapはV→−V,k→−k,C→Dである。

C=v≠0について商C³/Cvへ写すとUv=beta v、
im(−T+3V+kI+beta I)⊂Cv。
v上の作用からVv=nu v、traceからbeta=(2nu−2k)/3、
商上でTbar=3Vbar+(k+beta)Iを得る。
V=H+Nの固有lineはa方向とc方向だけで、b方向はJordan generalized vectorである。
全必要支持は

```
Ca: t=3,k=−2;          Cc: s=3,k=3t−5;
Da: t=−3,k=−2;         Dc: s=−3,k=−3t−5.
```

四つのC/D交差は全て空：tの±3矛盾、sの±3矛盾、又はk=−2と−14の矛盾。
全parameter値でC=0又はD=0が強制される。
FULL P2がnormal rowを持たないことから、そのhyperplane因子をhigherで修復できない。

## 結論と射程

旧全normalizer正逆を戻すと、normalSS MIXH／MIXN／NILREG全二parameterと
mixed全三parameterの指定rank四・被約既約quinticは全空。
これは旧G次元上界を全点帰属へ外挿した結果ではない。
Cartan、common-kernel/low21、normalnilその他のsourceは既別義務のままで、
このreceiptだけで全rank零plane又は全Eの完了を主張しない。新成分増分は零。
今回の認証は純証明／既証明書の読解であり、新operator算術再実行とは区別する。
