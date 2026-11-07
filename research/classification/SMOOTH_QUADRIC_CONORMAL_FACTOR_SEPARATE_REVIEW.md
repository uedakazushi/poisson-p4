# 滑らかなquadric全支持のhyperplane因子補題：独立監査

2026-11-07。reviewer `root`、producer `plane_collision_review`。
**全R、全normal-cone型を含む命題の幾何証明を独立確認し、PASS。**
新CAS、有限係数標本、接空間の次元推測は使わない。

## 1. 正確な命題と全Euler lift

smooth quadric S=(h,q)⊂P4 がPoisson invariant、conormal bracketが零なら、
そのspan H=(h=0)もPoisson invariantとなる。
これはrestriction rankの仮定を要しない補題である。
非正規支持への帰結ではrank二の既受理conormal補題と、被約Dに対するrank零のmodular補題を別々に用いる。

任意のquadratic raw liftのnormal contractionがS上のEuler lineにあることは、
Euler列をそれぞれO(1)、O(2)でtwistしたkernel条件である。
従って係数λ、μは局所有理函数ではなく、正確にH0(O_S(1))とH0(O_S(2))のglobal section。
S⊂P3のquadric完全交叉列、又はSegre monomialで両restriction mapの全射を確認した。
q=x^TQx、μ=x^TBxをQ,B symmetricで書き、
W(h)=λ、W(x)=Q^-1Bx/2を同時に選べることを直接計算した。
W(q)=μであり、W(h)の指定はqがhに依らないので干渉しない。

## 2. intrinsic bracketの補正を落としていない

Π'=Π+E∧Wの全normal contractionはordinary I=(h,q)に入る。
chart Xj非零でu=h/Xj、v=q/Xj²と置いた実bracketは

`Xj³{u,v}={h,q}−(2q/Xj){h,Xj}−(h/Xj){Xj,q}`。

両補正項はI²、左辺もintrinsic conormal可換性からI²である。
ゆえにstrict homogeneous {h,q}がsheaf I²に入るという矢印は成立する。
raw liftがstrict Poissonであるとは仮定しない。degree零函数上のEuler gaugeは零であるため、
上の計算は元のprojective bracketそのものである。

prime regular sequenceのgr_I各pieceがR/I上freeなので、I²はI-primary。
g∉IならR/IとI/I²でのg乗法が単射であることから直接確認出来る。
I外のlinear formを一つ選べば、同じprimary性からI²のirrelevant saturationもI²。
従ってsheaf membershipからordinary homogeneous idealへの移行に穴は無い。

## 3. 全rowの次数とactual Pfaffian

新liftの全h-rowはdegree二かつIに入り、必ずh*linear+q*constantで尽くされる。
{h,q}∈I²をmod hへ送るとqΣa_i q_i∈q²。
degree三対四なのでΣa_i q_i=0、Qのrank四から全constant a_i=0となる。
chosen submatrixのrankではなく、全h-rowを戻した結論である。

Hのtangent dimension三を使ってrank≤二を得ることも出来るが、
actual Euler Pfaffianの各項を直接見ると、Eのh成分又はh-rowを必ず一つ使うので全項がhで割れる。
E∧(Π+E∧W)²=E∧Π²により、同じ実五次因子が元のπに強制される。
rank四・既約被約Dに矛盾する。

## 4. 全適用範囲

非正規曲面成分SはHamiltonian flowで保存される。
rank二restrictionではnormal一次Pfaffian比較とbundle写像の延長がconormal可換性を与え、
rank零ではlocal modular vector fieldによる全重複度補題が与える。
どちらも既受理の独立証明を参照した。
零曲線Rを一切用いないので、可約・非被約R、vertical因子、square normal symbolも全て含む。
従来のintegral-R排除と四root-template還元を成功理由として混ぜる必要は無い。

quadric coneではqのradicalがあり、Q^-1と最後のa_i零の二段階が使えない。
coneを含めた全quadric、他の全G/E又は成分数上界はこの命題から主張しない。
旧producer/primary packetsのbytesは全て不変である。
