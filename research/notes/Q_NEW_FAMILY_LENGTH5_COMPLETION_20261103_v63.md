# 新25次元成分の一般長さ5と n111 における全完成環

2026-11-03。generic な二法線 quadratic rank1の全関数体証明と独立監査が完了したため、
先行する条件付き補題の仮定はすべて満たされた。
新成分の一般横断 Artin 長さは**正確に5**、幾何的一般 Hilbert functionは(1,2,2)である。
n111 の実 full projective Poisson 完成局所環は

    C[[s1,...,s25,u,v]]/(uv,u³,v³)

に同型である。被約成分はこの点で滑らか、全環は Cohen–Macaulay・type2でGorensteinではない。
nilradicalの三乗は零で、二乗は非零である。
この点の全完成環と一般長さの結果であり、成分の全点で同じ完成環を持つとは主張しない。

## 実有限加群と rank の一致

actual family・群・scalarの26次元像を横断する cut 座標を使い、
R=C[[s1,...,s26]]→A→R という恒等合成を持つ全 canonical Poisson 完成環 A を取る。
独立監査済みの全126方向・全350式の閉 fiber は長さ5であり、完備 Nakayama により
R^5→A が全射となる。ここで26方向はすべて群方向ではなく、形式積やflatnessを先に仮定していない。

generic Jacobi rank98により、新しい26次元 affine 成分の generic 横断埋め込み次元は2である。
さらに generic augmented rank99の exact 証明は、intrinsic 二法線の二次関係空間の rank1を与える。
普遍 tangent を actual family方向で微分する降下補題も独立監査されている。
従って generic branch の Hilbert 値は少なくとも(1,2,2)で、長さは5以上となる。

K=Frac(R)とすれば augmentation が選ぶ branch の residue field はKで、
この長さだけで dim_K(A⊗K)≥5となる。一方R^5→Aから上界5を得る。
全射の kernel は整域R上の自由加群の rank0部分加群なので零であり、AはR自由 rank5となる。
これは閉 fiber の長さだけを用いた誤った kernel basechange/Nakayama 論法ではない。
追加の generic rank 下界が Tor や base torsion の可能性を除いている。
全generic rankを同じ branch が使うため、一般長さは厳密に5となる。

## N/N²を先に扱って乗法を決める

N=ker(A→R)と置く。split augmentation と自由性から N はR自由 rank4である。
generic には長さ5、埋め込み次元2、二次部分の次元2なので N_K³=0である。
AのR-torsion-freenessにより N³=0を全完成環へ戻せる。

N/N²は閉 fiber で二生成、generic rank2である。R²→N/N²の全射の kernel は
自由加群の rank0部分加群なので零となり、**N/N²が先にR自由 rank2**と分かる。
したがって0→N²→N→N/N²→0はR加群としてsplitし、N²も自由 rank2となる。
N²の特殊化を未確認のまま閉 fiber の二乗と同一視していない。

u,vをN/N²の基底のliftとする。N³=0なのでN²はu²,uv,v²で生成される。
Sym²(R²)→N²のkernelは自由 rank1で、二次形式qが生成元となる。
閉 fiber のbinary quadratic discriminantは非零なので、qのdiscriminantはRの単元である。
Rは複素数上の完備環なので、Henselによる平方根と可逆な二変数線形変換でqをuvへ移せる。

R[u,v]/(uv,u³,v³)→Aは全射で、両方R自由 rank5、閉 fiberでも同型なので同型である。
これが全 canonical 完成環の同定である。
非零 cone の自由scalarを一つのbase座標に取り、射影化で固定すれば、冒頭の25変数baseの完成環となる。

## 帰結と保存

nilradicalは(u,v)、被約環はRであるため、この完成環には極小素が一つしかなく、
被約成分はn111で滑らかである。R上の自由性はCohen–Macaulay性を与える。
Artin因子のsocleはu²,v²で張られる二次元空間なのでtype2となる。
一般横断長さ5は成分の一般多重度であり、適当な非空開の幾何的一般 Artin環は同じ長さ5を持つ。
大域的に自明な積や全境界点での正規形は主張しない。

元の条件付き原稿は変更せず保存する。新しい無条件 receipt は
`reports/Q_NEW_FAMILY_LENGTH5_COMPLETION_20261103_v63.json`、依存検証は同名小文字Pythonである。
実閉fiber・普遍kernel・generic二次rankの各独立PASSと、有限-rank論法を区別して固定する。
