# concurrent cubic・m4 squarefree special-rank2全源の独立排除

**独立PASS。** rank二restriction planeの零cubicが三相異concurrent lines、
実退化五次のnormal order四、leading binary quarticがgeneric squarefree、
I恒等零又はJ恒等零かつcoefficient map rank二の場合、全rank四sourceは空である。
旧粗上界31,31,31,32のfalseを変更せず、全必要sourceの共通外積因子により排除した。

既受理のWitt/Hankel theoremのcommon-kernel/secant、compression、I零mixed、
I零三相異の全四型を保持した。
前三型のdiscriminantは二つの線だけを使用してよく、第三線を捨てていない。
各stable pencilは少なくとも二つの相異discriminant rootsを持つため、その二つの
pullback線の交点がcoefficient mapのkernel点である。両線はRのconcurrent線なので、
その点はRの共通点に等しい。従ってleading係数をx,yだけに取る操作は
z係数を任意に捨てる操作ではない。
rank二coefficient mapの全base線形対は、全concurrent三線を送るPGL3、normal GL2及び
非零振幅により公開normal formへ入る。I零三相異では全三線を使用する。
normal formは全leading pencilであり、個別のhigher点への特殊化ではない。

f=xy(x-y)、normal u,vとして、各全leading generatorに対する
(f_y dx-f_x dy)G+(Mn)partial_n G-tr(M)G=0の全係数を独立再構成した。
normal matrix Mの全八未知のrank八、全equation/RHS、唯一解を確認した。
解を全十五base-linear*normal-quartic monomialへ作用させた全210*15行列の
完全kernelは2,2,2,1次元で、非零minorと独立kernelを照合した。

offdiagonal又は非log接続を初めから除外したのではない。
generic squarefree binary quarticのtwisted GL2 infinitesimal stabilizerは零である。
trace scalarは四次式とdet twistの差二で非自明に作用し、sl2部分は四相異根を固定するため零。
したがって一つのparallel connectionを全leading tensorで確認すれば、別の接続との差は零。
これがMz係数も零とする理由であり、z依存の未知接続を捨てる仮定ではない。
前三型の全二非零振幅、第四型の全共通振幅でも同じ接続がparallelなので、全振幅に適用できる。

全P0をpartial_z wedge W0として手式から作り、全150係数、全350 Jacobi、実Pf零を照合した。
normal-normalの任意u2,uv,v2項を加えた実D2はそれぞれ3f倍になることも確認し、
raw normal-normal bracket零がnormal order四により強制されることを確かめた。

全Euler quotient sourceはweight一三十、二九であり、その基底も先行の独立実装から構成した。
全350のweight一JacobiとD3全二十四行の374*30行列を再構成すると、
全四型でrank十八、完全kernel十二。
全十二fieldはzと対を作るbracketだけを持ち、coefficientは全てz非依存である。
従って任意P1はpartial_z wedge W1、[partial_z,W1]=0、Jac(P1)=0となる。
これは全kernelの等式であり、選んだsliceから全点へ外挿したものではない。

同じ固定P0に対する全350*9 weight二Jacobi operatorを再構成した。
非零columnは(x,y)のnormal-quadratic三列0,1,2だけであり、保存された原Jacobi行minorがrank三を与える。
したがってJac(P1)恒等零から全P2は残りの六standard columnだけとなる。
これらもpartial_z因子とz非依存性を持つ。
追加leading-shapeや高weight Jacobiを課す前の全必要源で既に

    P=P0+P1+P2=partial_z wedge (W0+W1+W2),  P squared=0, actual D=0

である。従ってrank四、非零generic-squarefree leading termと矛盾する。
全original quadratic compatibility matrix零という保存値とも一致した。

normal simultaneous centralizerの全8*4行列と完全kernelも再構成した。
前三型は二次元diagonal、第四型は一次元scalarであり、正しいgroup fiberは11,11,11,10。
存在しない第四型の二次元torusを引いていない。
最終排除はrank二の構造であり、このgroup次元や粗次元上界の改善を必要としない。

checker `code/classification_plane_concurrent_review.py`、SHA256
`6f400f9e18355aff275ec3b409ea0bff684a8de605541166711769998bbc7703`。
job `research/classification/jobs/plane_concurrent_special_review.json`、SHA256
`c7579cf24667f26be16ce148a2ddeacafef9447d581f9e78d1c14ea939efc97a`。
gate PASS後1.135395776秒、最大RSS21000KiBで全四源PASS。
report `reports/CLASSIFICATION_PLANE_CONCURRENT_SPECIAL_REVIEW_20261107.json`、SHA256
`ee45a4ab177f58380b95ca52208780f49e2825b74eeccdda3129fa97602f4e9a`。
producer/CASをimportせず、以前の独立Fraction実装だけを再利用した。
producerの全旧false、111関連packet、旧reviewerは一切変更していない。

generic j及びrank一の全concurrent squarefree枝は、既受理essential3/Saitoの
ordinary四重点mu九と全tau十二の矛盾で別途空である。
それと合わせて全concurrent m4 squarefree rank四sourceが空となる。
重根、m2/m3、非被約R、他の支持面、全Eはこの結論の補集合として保持する。
