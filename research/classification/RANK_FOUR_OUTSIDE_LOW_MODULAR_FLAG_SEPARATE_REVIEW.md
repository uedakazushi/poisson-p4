# 低rank marked modular旗外の全25次元上界：独立純査読

2026-11-08、別作者 `/root/plane_collision_review`、PRIMARY E。
対象は `RANK_FOUR_OUTSIDE_LOW_MODULAR_FLAG_DIMENSION_BOUND.md`
（SHA256 `ab9bd9f3272575b5f8915821b0e1c6436adebc24bccd4dda559dfaac1f029bf8`）。
**全case合成と新しい nil22-image fiber 証明を独立に手読解して PASS とする。**
新CAS・点・Jacobi行列・minor・kernel・jet・GBは実行していない。
既source定理の全式による独立証明書は再実行せず、元の適用仮定を保持して使った。
入力35ファイルのhashは全一致した。文書のhash検査は数学的証明の代用ではない。

## 1. 正則な閉旗

rank四開集合上ではactual Dは非零である。canonical divergence-free liftの
Poisson hyperplane ℓに対するcolumn-sharpは ℓB と一意に書け、
Bはtrace零、B(ℓ)=0である。これはprojective tangencyから
Y=ℓA−A(ℓ)E、divY=ℓ trA−5A(ℓ)=0と計算すれば従う。
よって Nℓ=5B はℓ-Casimir liftのdivergenceで、係数にregularに依存する。
ℓのunit座標chartにおけるmonic exact divisionを使うため、未証明の有界な
normalizerやfactorの極限を仮定しない。

したがって ℓ²|D かつ rankNℓ≤1 はmarked factor incidence上で閉じ、
双対P4へのproper性により存在条件 L はrank四開内で閉じる。
Pf零点に任意の因子を付けたincidenceは使わない。
既CA/CBの全磁気係数を保持したNℓの一行表示から、それらのrank四閉包全体はL内である。
固定prior129のrank四次元>25はこの二28成分だけであり、後続の新採録成分も25以下である。
rank二以下の閉成分はrank四開と交わらない。

## 2. 重複Hのない全sector

五の因子分割を全て保持する。normal actual Dは全T39、integral nonnormalは
修正版全pointwise合成へ入る。被約可約1+4、1+1+3、1+2+2では、全normalizer、
surface/plane/frame、source振幅を戻した既知閉和又は大域像≤24という結論を使う。
これはprimitive βの非優越そのものではない。
特に1+1+3の五非平面残源は、それぞれ全sourceと実preserving群による
23/24/23/22/22上界であり、平面十三配置は全六磁気を含め23以下である。
1+4の全nil十型にも固定cofactorの全源と実field orbitを戻した上界を使う。
2+3、1+1+1+2、五Hは既受理の全pointwise閉和へ入る。
非被約で重複Hがなければ正確にH+2Q、Q prime rank≥3で、同じく全被覆済みである。

これらの既知閉和のうち25を超える可能性はCA/CBだけだが、その全rank四部分はL内。
よってL外には全像≤25が得られる。この行でquartic-openの全源を使ってよいことと、
次節のmultiple-H入口へその表を転用してはいけないことを区別した。

## 3. marked Hからの全90係数入口

H-Casimir liftで P=C+L+Q、Z=c+Ax、trA=0、PfQ=0 を保ち、
全Jacobi、[Z,P]=0、Q∧Ax=0を使う。Cを非退化に固定しない。
N_Hのmatrix rankは rank[c A] であるため、A=0とrankA1/c-imageは全体がL内。
その大きな全sourceを小さい部分sourceへ置き換えていない。
残るrankAと全Jordan/coker型はproducer §§4–7に尽くされる。
非nilpotent rank≥2の全period、zero-eigenvalue、translation、linear resonanceは、
全低次項を保ったactual polynomial familyでMPSへ送る既定理を使う。
SNC初期形だけの一般点から全sourceを推測していない。

nilpotentの全像上界は次の算式であり、固定symplectic値やβを引かない：

| affine modular型 | 全固定source上界 | field orbit | H4−scale1を戻した上界 |
|---|---:|---:|---:|
| N4/image、五源 | 4 | 15 | 22 |
| N31/image、三源 | 10 | 12 | 25 |
| N31/long-top、六源 | 5 | 14 | 22 |
| N31/short-top、六源 | 6 | 13 | 22 |
| N22/image、十一canonical源 | 12 | 10 | 25 |
| N22/top、八源 | 7 | 12 | 22 |
| N211/long-top、十七源 | 12 | 10 | 25 |

N4/nonimageはPf零。N211/shortは対応する全long源への実曲線で戻す。
N31/N22五型の根拠は最高Pf零を保持する元90係数reportであり、
quartic-openの型5/6/7の表ではない。n7-unitの全十源とn7零half/thirdの
全十源を分け、unit条件を零境界へ課していない。
field orbitはGL centralizerとtranslation image、全coker orbitで数える。
非零tensor尺度は各nilpotent affine型へactual chain-weight共役で戻せるので、
全affine飽和はconeでありprojectivizationの−1は正当である。

## 4. 新nil22-image上界の独立確認

元全必要十分表示は M=−adj(S+HessD3/2)+p⊗V、F0は全二次六係数である。
actual同時GL2作用(q,p)→(Aq,Ap)から
D3'=det(A)²D3(A⁻¹p)、S'=det(A)²A⁻ᵀSA⁻¹、V'=AV(A⁻¹p)、
F0'=det(A)F0(A⁻¹p)を得る。
十一canonical源の自由次元は6,6,6,6,6,8,6,6,9,7,9。
GL2四次元を戻すだけで、二つの九変数源以外は12以下となる。

triple_CではD3=x³,S=0、V0=(0,n)、V1の下行(c,d)。
A_r=diag(r²,r)はD3とSを厳密に保ち、(n,c,d)→(rn,c/r,d)である。
quadratic_BではD3=0,S=diag(1,0)であり、A_r=diag(r,1)は
それを保ち、(n,c,d)→(n,c/r,d)。両方とも全F0六係数を可逆に動かす。
この二式はn,c,dのどれも反転せず、全零境界にも同じactual作用が存在する。

tensor写像GL2×sourceのfiberは
(G,s)→(G A_r⁻¹,A_r s)を含む。第一factor上の作用がfreeなのでfiberは必ず一次元。
特定tensorの安定化群が一次元であるとの仮定は要らない。
従って両行とも4+9−1=12以下。D3=S=0はdetM=0でrank四外であり、
Hessianの全rankとbinary cubicの全root型は旧独立必要十分分類が保持する。

## 5. 全上界から成分性への矢印

全因子分割は有限である。非nilpotentの連続spectrumは全MPS閉和へ送り、
nilpotentの全有限sourceを上表で押さえたから、L外は有限個の像閉包dimension≤25で覆われる。
L外と交わる全空間既約成分Cではその交わりが稠密開なので、有限閉和の一つにC自体が入る。
ゆえにdimC≤25。一般βや未知成分dim25–27の下界は一切使わない。

別作者が既に認証したZ_M1,Z_M11は全係数逆、全StabX fiber、scalar stabilityと
一般markedH一意性により正確な既約25次元閉集合である。
それぞれの一般D=4hL、3hL1L2ではhだけが重複Hで、
N_h=h∂x+x∂yはrank二。従って一般点はL外にある。
これを含む成分CはdimC≤25=dimZ、よってC=Zとなる。
二者相互およびpriorからの区別は既入口のproper因子条件/marked-rank監査を使用する。
**両whole像閉包の全空間成分性を受理する。** 全境界もその閉成分内である。

## 6. 正確なscope

これは被約支持の成分性であり、smooth、局所成分一意性、横断multiplicity、
全schemeの非被約構造を主張しない。一点の接次元28/27を変更しない。
L内全sourceの次元上界25、全Eの完成、成分数の上界も主張しない。
本査読者自身が以前作成した全G/非正規合成は既別作者の正式受理依存として使い、
その旧算術を独立再実行したと称していない。新しい独立部分は上の全仮定接続と
二つのfree Cstar fiber、およびそれからの全component-support矢印である。
