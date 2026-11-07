# rank零plane・m二・係数rank二・M(p)=0の全first-symbol正規化

2026-11-07、主義務G。producer純幾何還元。
全normal-weight Jacobiの有限jobに先立つsource及びactual逆変換の証明。

`PLANE_RANK_ZERO_MULTIPLICITY_TWO_FIRST_JET.md` のO(1)型で
M:C³→sl2がrank二、かつinclusion Uのpoint pについてM(p)=0の場合を扱う。
rankM三、M(p)≠0、O(0)型、m_D≥3は別補集合として保持する。

## 1. 全Lax式と二つのimage-plane型

U=∂aを選ぶと、M=bM1+cM2、M1,M2は独立で、W=imM⊂sl2は二次元。
weight零sourceは

`Π0=T∧E_n+U∧(M n·∂n)+q_R(n)∂z∧∂w`。

全projective Jacobiは前掲Lax式

`T(M)−t M+[R,M]=0`、`[T,U]=−tU`

と同値である。
trace-free liftを用い、T(a)=t a+r_b b+r_c c、T(b,c)=A(b,c)と書く。
Wへのad_RのrestrictionをJとすれば、

`R∈N_sl2(W)`, `A=tI−J`, `3t=traceJ`。              (1)

r=(r_b,r_c)は初めは全二係数を保持する。
actual base shear a'=a+ℓ(b,c)は

`r→r+ℓ(A−tI)=r−ℓJ`                               (2)

と作用する。逆はℓの符号を戻すだけであり、base/normal gradingを保つ。

sl2のnondegenerate Killing pairingでWのorthogonal lineを見ると、
WのAdPGL2 orbitはnonisotropic lineの場合とisotropic lineの場合の二つだけ。
それぞれ

`W_nd=<F,E>`、`N(W_nd)=C*H`、
`W_B=<H,E>`、`N(W_B)=<H,E>`

と取れる。H=diag(1,−1)、E=[[0,1],[0,0]]、F=[[0,0],[1,0]]。
W_BはBorel subalgebraである。
nondegenerate orthogonal basis又はWitt basisを選ぶことで全可逆normal transformationが得られ、
M1,M2のbase basis changeは可逆なので、この二whole orbitはarbitrary sampleではない。

R=0ならJ=t=A=0、TはUのmultipleとなってleading normal quadratic Pfaffianが零、
m_D=2に反する。従ってR≠0を保持する。

## 2. nondegenerate Wは一つのwhole symbol

M=bF+cE、R=αH、α≠0と取る。
Jはinvertibleなので(2)で全rが消える。
tensor scalar σ=1/αを掛け、actual base change
`a'=a, b'=σb, c'=σc`
を行えばMのcoefficientを再び一へ直し、R=Hとなる。
全変換は可逆である。
全symbolは

`U=∂a`, `M=bF+cE`, `T=2b∂b−2c∂c`, `K=−2zw`。     (ND)

これが全nondegenerate W、M(p)=0、m二の一orbitを代表する。

## 3. Borel W・semisimple Rではrow parameterを全て保持する

R=αH+βE、α≠0の場合、normal unipotent conjugationでβを消し、
同時にM1,M2のbase basisを変える。
例えばAd(exp(sE))R=αH+(β−2αs)Eなのでs=β/(2α)である。
M=bH+cEと直すと、J=diag(0,2α)、t=2α/3。
(2)でr_cを消せるが、r_bは消せない。

tensor σ=3/α及び `a'=a,b'=σb,c'=σc` でR=3Hとなる。
この変換でr_bは**変化しない**。
従って残るscalarをρと書き、全ρ∈Cを保持して

`U=∂a`, `M=bH+cE`,
`T=(2a+ρb)∂a+2b∂b−4c∂c`, `K=−6zw`              (BSSρ)

という一parameter whole familyを得る。
ρを任意に0又は1へ取った点から全sourceを推論しない。
この全familyについてQQ(ρ)上のoperator、selected-maximal-minorの全例外根を扱う。

## 4. Borel W・nilpotent Rの二whole orbit

Rはnonzero nilpotentなのでnormal conjugationでEへ直し、M=bH+cEとする。
J(H)=−2E、J(E)=0、t=0。
(2)でr_bを消し、r_c=εを保持すると

`T=2b∂c+εc∂a`, `K=−w²`。

ここではε≠0を一へ直す実変換がある。
tensor scalar σを掛け、normal diagonal conjugationでE→σ^−1E、H→Hとし、
baseを `a'=a,b'=σb,c'=c` と変える。
するとM=b'H+c'E、R=E、2b∂cはそのまま、ε→σεとなる。
σ=1/εを選べばε=1、逆はσ^−1とnormal/base逆変換。
normal diagonal matrixには必要なsquare rootを複素数体上で選べる。

従って全nilpotent branchは

`U=∂a`, `M=bH+cE`, `T=2b∂c+εc∂a`, `K=−w²`, `ε=0,1`。 (BN0/BN1)

零と非零の補集合を明示したwhole-orbit分類である。
semisimpleのρへ同じrescalingを流用していない。

## 5. 全parameter次元と最弱fiber目標

fixed planeでpを二parameter、Mのpをannihilateする二columnを六parameter取ると、
rank二coefficient mapのbaseは八次元。
nondegenerate Wはopen、Rのnormalizer一、row二なので(ND)全affine symbol sourceは十一。

Borel WはGrassmannのdiscriminant一式でbase七次元。
semisimple Rはnormalizer内open二次元、row二なので全(BSSρ)源は十一、
ρを固定した層は十。
nilpotent Rはnormalizer内一次元なので、row全二で十。
ε=0の閉条件で一つ下がり(BN0)は九、ε≠0の(BN1)は十。

全30 weight一及び全9 weight二operator kernelの和をdとする。
normalizationはgradingを保存するactual変換とtensor scalarなので、kernel dimensions及び
full nonlinear higher-fiberの必要上界は全orbitで保存される。
fixed plane sourceへplane choice六を加え、tensor projectivization一を引くと

| 全symbol層 | affine symbol base | G像upper | G≤24に十分なd |
|---|---:|---:|---:|
| ND | 11 | 16+d | 8 |
| BSS一般ρ | 11 | 16+d | 8 |
| BSS固定ρ | 10 | 15+d | 9 |
| BN1 | 10 | 15+d | 9 |
| BN0 | 9 | 14+d | 10 |

を得る。
全350 coefficientを使う有限operator jobがこの弱い上界を満たせば、
非線形higher Jacobiを追加せずこのwhole first-symbol stratumのGを閉じられる。
満たさない場合はexact kernelと上界の限界を保存し、full80-variable GBへ拡大しない。

このnoteは上界をまだ計算しておらず、全Poisson支持又は全Gの完了とはしない。
