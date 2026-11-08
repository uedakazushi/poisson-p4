# constant TRIPLE：同一軌道graph入口と全T40排除の独立手査読

2026-11-08、PRIMARY E_multiple_hyperplane、担当 `/root/exceptional_coverage`。
**独立手読解・手展開及び既exact証明書のreadback。新算術実行は行っていない。**
producerの唯一の実行結果はT39/T40とも陰性であり、T39帰属の陽性証明とはしない。

## 1. 全四源と同一actual tensor

全有限還元とrootの別作者査読を読んだ。残源は全p,d,l,Lに対する

```
f=x³+3xy+p x−3y/2+r+d,
Yx=x/3+y−l,
Yy=−xy+l x+y/6−r/3+v,
Yr=−3y²+xr+L x+(3l−p−5/4)y+r/2+w,
v=2p/9+l/2−L/3, w=d+lp+p/3+3l/4−L/2,
P=∂t wedge Y+Jac(f).
```

r'=r+px+d、t'=t+(p/3)yはdet一の実affine変更で、l'=l+p/3、L'=L−d−p/6。
続くf0保存写像とt-shearは三角形の実逆を持つ。
A=3l/2、alpha=(L−2A³−2A²−A)/3を使いl,Lを消す際にもsource係数を反転しない。
したがって四parameter域の全点が同じactual orbitに属する。
元のbB-unitを正規化した非零対角尺度以外に新しい開集合は無い。
この事実はsource点の局所smoothnessを必要としない。
元のcube、triangular、rank零等の補集合は既有限還元の別known枝であり、この一軌道と混同しない。

## 2. 固定r変更とfaithful coframe

canonical六括弧へr_old=r_new+1を代入すると、Ptyに−h²/3、Ptrにhx+h²/2が加わるだけ。
これは同じtensorの実affine表示変更であり、別source点又はrank用の標本ではない。
新原点の六skew値は(0,−1/3,1/2,1,3/2,0)、Pfは1、detは1。

```
S=((1,0,0,0),(0,0,−3,0),(3,1,0,0),(0,0,1/2,1/3)).
```

手でdetS=1、SKSᵗ=Jを確認した。
具体的にはa=t,b=−3y,c=3t+x,d=y/2+r_new/3で、
{a,b}(0)={c,d}(0)=1、他のcross値は零。
逆t=a,y=−b/3,x=c−3a,r_new=3d+b/2、r_old=r_new+1も整合する。
これらの分母は固定非零有理数であり、元四parameterの特殊値を除かない。

旧caller/S61 adapterを読んだ。Euler24を線形にquotientし、実Darbouxと
perspective Ga4を施し、全150の復元、元134式、全射影Jacobiを検査してから
全45個のQ51/Q52差を評価する実装である。
単なる線形部分の一致又は134式への代入だけを成分所属とはしていない。

## 3. 完成した唯一のgraph結果と独立した一残差検算

root GO後のproducer reportは2.63919871699909秒、全150復元／全134式PASS。
実coframeは上のS、perspective=(0,2/9,0,3/2)。
全45残差中T39は31本、T40は32本が非零。
producer結果noteのT40も31本という記述は転記上の誤りで、元execution reportの全値は正しい。
executionのall_passed=trueは計算完走を意味し、graph所属の真偽とは別fieldである。

保存された全61座標とQ51公式から、u2残差を手で独立検算した。
u2=6、u13=1/3、u27=−9、u35=1、u51=−9/2。
当該公式の他の積は全て零なので

```
Q51_u2=(u13*u27)/4−(u35*u51)/2=−3/4+9/4=3/2,
u2−Q51_u2=9/2 ≠ 0.
```

全45式を別算術エンジンで再実行したとは主張しない。
実frame・その逆・完全な旧adapterのscopeを手監査し、producer全証明書の入力hashと
全残差をreadbackした結果と、上の一式の独立手検算を区別する。

Q51の全A16閉graphが被約T39のS61閉包に等しいことは既普遍切断定理と独立査読を用いる。
さらに既受理の全評価fiber定理は、connected Sp4により全J-fiberが既約、
perspective四pivotの定数det36によって全fiber=A4×S61 sectionとなることを証明する。
T39/T40は各dim25なので全sectionが既約dim16、含まれる閉A16 graphと全一致する。
従って今回のfaithful非退化入口の陰性は、当該実tensorが**全global T39/T40の外**であることを与える。
別評価点又はframeの探索は不要。全四parameterの一actual orbitも同じく両成分外となる。
軌道閉包の全境界点まで外だとはしない。

## 4. T40の別のwhole幾何排除

GENUS_ONE_RECOGNITION §6の引用を、その依存IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE §4と
GENUS_ONE_SEPARATE_REVIEW §§5–6まで戻って点検した。
必要な入力は「非正規支持がscroll」という一般論だけではない。
既known q52の一次資料MP Lemma3.2.3は、elliptic curve上のstable rank二degree五bundle Wに対して
**projective rankzero locusの被約台がsmooth embedded scroll P(W∨)**と述べる。
これは二次元で、full q52 construction/PGL像はT40に稠密である。

projective Poisson空間上のuniversal projective Zero incidenceはclosedで、
P4因子を忘れる射はproper。そのfiber dimension上半連続性によりdimZero≥2はclosed。
従って全T40がこの条件を満たす。unstable特殊点から逆向きに一般幾何を推論していない。

TRIPLEはh=1でPxy=1だから零点無し。
h=0のEuler三次式はx²y、x(3y²−xr)、3tx²−2xyr+3y³、3x³。
手でこの四式と逆向きの消滅を確認すると、全零台は厳密にh=x=y=0の一本のP1。
ゆえにT40の全closed条件に反する。零schemeの厚み又はT40の特殊parameter resolutionは必要無い。

## 5. 判定

同一軌道化、固定translation、Darboux／逆、faithful全graphへの含意、
一非零Q51残差及びwhole T40零曲面閉性を独立手査読PASSとする。
producerが実行した全graph証明書は陰性のまま受理する。
**T39帰属は得られず、既一軌道は全T39/T40外。新成分又は全E完了は主張しない。**
この固定graph経路を停止し、別origin／point／normal／jetを追加しない。
主稿、共有コード、元report/job、global台帳及びGitは変更していない。
