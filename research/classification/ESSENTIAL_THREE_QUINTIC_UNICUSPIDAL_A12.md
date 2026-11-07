# 最大重複度二・一尖点五次のA12代表への全還元

2026-11-07、PRIMARY E。全本質三変数既約五次Poisson定理から、
対象curveは有理cuspidal quintic、全δ=6、全μ=τ=12である。
最大重複度二で尖点が一つなら、その型はA12となる。
ここでは全degree5 normalization sectionsからPGL3代表を一意にする。
Nambaの全normal form表やMoeの存在構成だけから一意性を推論しない。
新しいSaito/周囲/障害計算は行わない。

## 1. 全degree5 linear seriesの正規化

尖点の逆像をs=0とする。
接線との交点重複度は四である。
実際、三又は五ならfirst Puiseux pairは(2,3)又は(2,5)でδは一又は二となり、
六以上ならdegree五のBézoutに反する。
接線の残交点は一つのsmooth点で、その逆像をt=0へ移す。
接線のsectionをX=s⁴tに取れる。

他の二sectionのうち、Yを尖点でorder二、Zを尖点で非零とする。
YからXの倍を引き、ZからX,Yの倍を引き、それぞれのleading coefficientをscaleすると、
全sectionは必ず

    X=s⁴t,
    Y=s²t³+C s³t²+A s⁵,
    Z=t⁵+B st⁴+D s³t²+E s⁵                     (1)

となる。ここでは五parameterを全て残した。
残交点の位置を固定した後にsourceのpoleを持つ座標変更は使っていない。
degree5曲線の正規化mapなので、三sectionに共通basepointが無いことも保持する。

## 2. 有限Puiseux係数比較

pでr=s/t、den=1+B r+D r³+E r⁵、P=1+C r+A r³と置く。
affine座標は

    x=X/Z=r⁴/den,   y=Y/Z=r² P/den.

den(0)=P(0)=1である。
A12なら、xをyの解析的関数で補正したとき、最初の奇数次は十三である。
従って以下の奇数次五、七、九、十一は全て消えなければならない。
これは固定degree五curveの型を決める有限条件であり、Poisson変形の高次jet反復ではない。

まず

    x−y²=r⁴(den−P²)/den².

その奇数次五の係数はB−2Cだから

    B=2C.

この後の偶数次六をC²y³で消す。
式 x−y²+C²y³ の奇数次七の係数はD−2A+C³なので

    D=2A−C³.

ここでL2=C⁴−2ACと置くと、偶数次八をL2*y⁴で消した後の奇数次九の係数はE+AC²である。
従って

    E=−AC².

最後にM2=−A²+6AC³−C⁶と置くと、

    x−y²+C²y³−L2*y⁴−M2*y⁵

の奇数次十一の係数は正確に−A²Cである。
従ってA²C=0となる。

上の係数を直接照合するための中間恒等式も記す。
B=2C,D=2A−C³,E=−AC²の後に

    x−y²+C²y³ = r⁶ N(r)/den³

と書くと、Nのr²,r³,r⁴,r⁵係数は順に

    L2=C⁴−2AC,
    L3=2C L2,
    L4=−A²+2AC³+C⁶,
    L5=−6A²C+6AC⁴.

L2*y⁴の補正後は

    M2=L4−2C²L2,
    M3=L5+2C L4+(D−4C³−4A)L2.

M2*y⁵の補正後の奇数次十一は

    M3−3C M2
      = L5−C L4+(C³−2A)L2
      = −A²C.

であり、全て有限多項式恒等式として手で検算できる。

## 3. basepoint補集合を除くと一候補だけ

A=0ならB=2C,D=−C³,E=0となり、(1)の三sectionは全てtを因子に持つ。
これは正規化のdegree5 line bundleにbasepointがあることになり、対象曲線ではない。
従ってA≠0で、C=0が強制される。
残る全sourceは

    X=s⁴t,
    Y=s²t³+A s⁵,
    Z=t⁵+2A s³t²,   A≠0.                            (2)

sourceのsをλsへ変え、targetのscaleを取るとAはAλ³へ移る。
C上でλ³=−1/Aを選ぶと、全非零Aは

    (s⁴t : s²t³−s⁵ : t⁵−2s³t²)                      (3)

という一つのPGL3型になる。
A=0補集合を捨てたのではなく、上でbasepointとして排除した。

(2)はA≠0ならbasepointを持たず、三sectionは線形独立である。
map degreeとimage degreeの積は素数五なので、imageはdegree五、mapはbirationalである。
消去した四つの奇数次からpのδは少なくとも六となる。
算術種数六からδは正確に六で、pはA12、他の特異点は無い。
従って候補の存在も単なる形式germでなく実際のdegree5曲線で確認できる。

## 文献と残義務

(3)はMoe, arXiv:1511.02691v1, §6.3,p.82のC2 parametrizationと一致する。
同論文のsecondaryな構成を全族の一意性に読み替えてはいない。
Dimca–Sticlaru arXiv:1504.01242v4, Proposition3.6(ii)のk=0、Example4.4も
degree5・一Puiseux pair(2,13)・free(2,2)の曲線を与えると報告する。
そのhomogeneous定義式の一つは

    z(yz−x²)²−2xy²(yz−x²)+y⁵.

全十八係数Saito核と具体parametrizationの全係数一致は次の有限算術義務に残す。
固定Fの全Poisson源は既済の全磁気inverseで一つのactual23軌道へ入るが、
その軌道閉包の成分性・既知帰属や新しい成分数を本稿だけから主張しない。
