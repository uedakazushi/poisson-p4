# Poisson部分多様体のincidence：現行研究ノート

2026-11-10 / v5。研究方針はRESEARCH_POLICY.md、証明全文はpapers/poisson_subvariety_incidence_ja.texとその共有TeXにある。

## 共鳴Logの直交性：何が機構か

V=U1⊕L⊕W (2,1,2)、q=H(u)+2zb(u)+cz²+w0²+w1²、HおよびqU非退化とする。重みξ=(1,1,-2)とW回転RからP0=2ξ∧R+ξ(qU)Jを作る。affine斉次P0の零スキームはP(W)という直線とP(U1⊕L)内の平面二次曲線C=(ξ(qU)=0)の非交和である。Cは非被約でもスキーム論的な線形包がP(U1⊕L)となる。したがってstrict so(q) stabilizerは二ブロックを保存する。

U1⊕L上で[ B, ξ ]=0を課すとB=diag(D,d)。直交性はDᵗH+HD=0, Dᵗb+db=0, 2cd=0。Hが非退化な2次元のためso(H)は一次元、非零生成元は可逆。c≠0ならb≠0でD=0。c=0でもqU非退化性はbᵗH⁻¹b≠0を強制し、bがskewなDᵗの非零固有値ベクトルになることは不可能。このため全非退化域でstrict stabilizerがb≠0ならCR、b=0ならCR⊕so(H)。

発散がintrinsic対称性へ入る線形Jacobi制約を使うと、b≠0では全延長がP0+tKRでrank2。b=0では重みと可換な座標変換で標準二直線模型に帰り、全延長Spec C[t,v]/(v(t+2))、rank4の実因子はQ+V(z)+V(cz²+w0²+w1²)。元の一変数q_λのminor 13824λ/(λ²−1)^4も保持する。今回の定理を特別なLogの全分類と同一視しない。

## 次に調べる特別な境界

Hが退化する場合、あるいはqUが退化する場合には、零スキームの成分と直交代数の層別を見直す必要がある。特にHの核方向とbの位置を不変条件とし、c=0や非等方・等方の分岐を調べる。追加のLie対称性があるだけでPoisson rank4延長が生じると判断しない。実JacobiとPfaffianの厳密計算で確認し、概念的説明へ戻す。

## 他の残課題と文献

一般4次del Pezzoの延長剛性は証明済みだが真の例外locusは未分類。P3の6 intrinsic成分の横断対称性も未分類。滑らかなQ3の一般Rat/Log/Affは一点/A1/二点の延長を持つ。LPT (DOI 10.1002/mana.201100354) はintrinsic分類、Mayanskiy (DOI 10.1016/j.jalgebra.2014.08.049) はsmooth cubic。新定理の文献上の優先権は未確認。

再現：code/subvariety_incidence/check_log_orthogonality.pyは五実例のstrict stabilizer rankとPoisson恒等式を実行し、research/subvariety_incidence/certificates/log_orthogonality.jsonと照合する。これは一般定理の証明の代わりではない。正本はGitHub main、生成PDFは非追跡、Actionsなし、Libraryは実push失敗時のみ緊急避難。
