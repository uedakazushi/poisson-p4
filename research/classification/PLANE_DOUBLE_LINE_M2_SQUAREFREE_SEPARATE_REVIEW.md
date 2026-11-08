# R=2L+M・m二 SF：全 flags／全 higher 空性の別作者監査

2026-11-08、PRIMARY G。`/root/quartic_sources` が root の次の二純証明を独立に手監査した。

- `PLANE_DOUBLE_LINE_M2_SQUAREFREE_WHOLE_FIRST_SYMBOL_REDUCTION.md`
  SHA256 `d49078dd973bcea6965943f0fb46031cbb85938473e3da9fff6443870b7f65ab`。
- `PLANE_DOUBLE_LINE_M2_SF00_RATIO_EXCLUSION.md`
  SHA256 `9a30c720a02ac1736edc8082b00500d48de397695e06561316d19850da329da1`。

全 pole-budget 入口、全 I/II・unit 補集合、全 literal Z、全 actualu² 比較、
全 normalizer inverse、全 ratio spectrum／Gauss 橋を再導出した。
**R=a²b・normal order二・generic SF の全被約既約 quintic source は空**を受理する。
新 CAS、operator、matrix、minor、GB、jet、parameter sample は実行しない。
先行の一fixed P0だけを全higherの存在又は帰属へ外挿していない。

## 1. 全 pole-budget と flags の補集合

`PLANE_NONREDUCED_M2_SPLIT_ROOT_ENTRANCE.md` を必要部分まで独立に読解した。
accepted discriminant 恒等式から G2=h ell1 ell2、primitive rows の M と h detM=f を得る。
flags は 00,01,02,03,11,12 の全六通り。
各 prime で primitive row の unit 成分を先に用いて theta_i の pole bound を得る。
M mod prime が rank一のとき、両 primitive rows が非零なので
adj(M)diag(1,0)M の少なくとも一つの off-diagonal が generic unit。
eta の pole order k は middle term で m+k、M^-1dM は高々 m、
scalar theta2 は off-diagonal に寄与しないので、全条件は pole(eta)≤div(h)。
f の非被約性や h と detM の衝突を落としていない。

h定数の 03/12 は eta が regular homogeneous degree−1 なのに Euler 値が非零で矛盾。
11 は h一次、closedness と Euler 値零から eta=0。
detM=ab の全 rank三／rank二 matrix 正常形では、全leftUを保持しても
actualu²係数は 3f k20 で、target h x²又はh U11U21 x²が f で割れず矛盾。
detM=a² では全 M=U[[a,t],[0,a]]V、primitive rows により U の第一列二成分 p,r は共に非零。
eta=0 の同一 diagonal character により actualuv係数は3f k11。
target h a[2pr t+(ps+qr)a] は a²h で割れず矛盾する。
t は a と独立、h の重複も cancel 後に同じ矛盾である。
従って全 03/12/11 が空、残 00/01/02 の irregular closed pool は原表の全 source を保持する。

## 2. full bad-vector 必要式と全 I/II

f=a²b、H=(a,−2b,c)、全 closed diagonal fields は
X=λH+(p a+q b)∂c、Y=μH+(r a+s b)∂c。
全 irregular d(c/a)、d(b/a) 成分を含み、log source への制限ではない。
全 actualuv 係数から λs=μq、γ=3(κ+λr−μp)非零。
Z と k20 はこの uv 係数へ入らない。

constant u root の全 P1 bad term v²C は、既約 D のため C非零が必要。
元 projective coefficient 式
`E∧((DY+κI)C)−(2Y−X)∧C+P_(C(f))=0`
を独立に三成分へ展開した。
特に

`Cb(κ−4μ+λ)=0, Ca(2λ−5μ−κ)=0, Ca(μ−λ−κ)=0, Cb(2r−p+2)=0`。

Ca非零なら λ=2μ,κ=−μ。
Ca=Cb=0,Cc非零でも `(κ−μ+λ)Cc=(5μ+κ−2λ)Cc=0` から同じ結論。
残 Ca=0,Cb非零では κ=4μ−λ,p=2r+2。
よって全非零 C を尽くす必要 union は
I:(λ=2μ,κ=−μ)、II:(κ=4μ−λ,p=2r+2)。
unit は Iで μ(2r−p−1)非零、IIで (λ−2μ)(r−1)非零。
Ca,Cb,λ,μの零を初めから除外していない。

00 の v root にも同式を使い、X/Y交換と κ→−κ を行う。
四つの交差のうち I/I、I/II、II/I は全て γ=0。
II/II のみ μ=−λ,κ=−5λ,p=r=−2、flatness は q=−s、γ=−27λ非零。
new c=c+q b/(3λ) で q,s が同時に消え、new c=λc,tensor 1/λ で λ=1となる。
anchor P_f と p=r=−2 は保存される。
全 inverse は minus shear と非零 scale の逆で、λ=0はunit違反として別に排除済み。
従って00全first-symbolは原noteの一fixed P0の全fiberへ入る。

## 3. 01全 literal Z と全 actualu² の監査

以下の signs は収縮順序ではなく、原noteが明示した cross の literal 成分で固定する。
`E cross Z=(bZc−cZb,cZa−aZc,aZb−bZa)=−3(f/B)(dA+A eta)`。
tracefree inverse は Z=−curl(E cross Z)/3。
この式から全四 Z を手で再導出した。初期の収縮記号の解釈差はこの明示 cross で解消され、
rootの literal Z や後続係数は変更していない。

B=a では eta=ρdloga−(1+ρ)dlogb、λ=μ=ℓ,p−r=3(1+ρ),q=s。
I は γunitに反し、II は κ=3ℓ,p=6ρ+4,r=3ρ+1,γ=−9ℓρ非零。
A=b では Z=−3ρb∂c、actualu² target比較は k20=0。
A=c では Z=diag(ρ+2,ρ−1,−2ρ−1)、actualu² の abc 係数は −9ℓ(ρ+1)、
target係数は −9ℓρ。ℓ非零により矛盾する。
A の両配置は actual Aut(R) と normal Borel の全逆を持ち、片配置の chosen point ではない。

B=b では eta=−dloga+d((ρb+σc)/a)、
λ=μ−σ,p=r−3,q=s−3ρ,σs=3μρ。
A=a は Z=(0,−3σa,3ρa) で、全u²比較は
`k20+3λρ−qσ=0, γ=3σ(1−p)`。
I では σ=−μ,s=−3ρ,q=−6ρ,k20=0、
target γ=3μ(r−4) と unit 式 γ=3μ(r+2) が一致せず矛盾する。
II は r=−5,p=−8,κ=3μ+σ,γ=18(μ+σ)非零。
flatnessで3λρ−qσ=0、従ってk20=0。
γ比較からσ=2μ、μ非零、κ=5μ。
A=c は Z=(−ρa,−3a−ρb−3σc,2ρc)。
fullu²係数は
`3a²b(k20+ρ−q)+3(1−p)a³+3(1−p)σa²c`。
flatnessにより他のabc係数は零。
target γa²c の a³ 比較で p=1、続いて γ=0となるので全空。

以上で01の残りは二つの全higher ratio source だけとなる。

## 4. 02は全六 quadratic 係数で矛盾

pole-budgetの全B2はa²又はab、etaは−2dlogb又は−2dloga。
両場合λ=μ=ℓ,q=s、Iはunit違反、IIのみ。
B2=a² は p10,r4,κ3ℓ,γ−9ℓ、B2=ab は p−8,r−5,κ3ℓ,γ18ℓ。
両ℓ非零を actual c-shear／scale で一へ戻しても、A2 は任意primitive quadraticのまま残る。

B2=a² は Z=3(Ac,0,−Aa)、β11=4dloga−3dlogb−d(c/a)。
actualu²比較は
`(a∂a+(4a+c)∂c−1)A2=(k20/3)a²`。
A2=αa²+βab+χac+δb²+εbc+ζc² とすると、
c²,ac,b²,ab の係数は順に ζ=0,χ=0,δ=0,ε=0を強制し、
A2=αa²+βab。a共通因子により primitive性に反する。

B2=ab は Z=3(0,−Ac,Ab)、β11=−2dloga+3dlogb−d(c/a)。
actualu²比較は
`(3a∂c+b∂b−2)A2=−(k20/3)ab`。
c²,ac,a²,bc の係数から ζ=χ=α=ε=0、
A2=βab+δb²。b共通因子でprimitive性に反する。
非primitive境界は元 UFD の他flagへ戻り、本flagの全 source に含めて誤解していない。

## 5. 全 raw36/P2九を保つ ratio 定理の再導出

constant u rootのnormal bracketはK=u²k20+uvκ。
ζ=v/u、k_eff=k20+κζ、base A(ζ)=Z−ζ(X−Y) とする。
全 higher mixed contractionの追加を B(u,v) と書くと、B は normal homogeneous degree一の rational vector。
raw36 の tracefree base derivative、全9P2のnormal column零により row divergence は
Ru=κu,Rv=−2k20u−κv。higherはこの二成分に寄与しない。

canonical strict Casimir lift Π=P+E∧R/5 の M=Π#dζ は

`Mnormal=(3/5)k_eff E_normal, Mbase=(A−(2/5)k_eff I)base+B`。

この全source恒等式から M(F)=0 を読み、homogeneous degree五を使えば
`[(A base+B)·∂base+k_eff(E_normal−2)]F=0`。
base'=base+T、T=−(k_eff I−A)^−1B で全translationが消える。
canonical M においても translation係数は
`3k_eff/5−(A−2k_eff/5)=k_eff I−A` なので同じ式である。
T は normal degree一、逆は base=base'−T。
有理 normal field K=C(u,v)で読む操作であり、これを global PGL とは呼ばない。

normal degree m、base degree5−m の作用素は
`Sym^(5−m)(A)+k_eff(m−2)I`。
原normal order二はこのtranslation後も保たれる。
m3/4/5が可逆なら全higher summandsは零、Fは固定leading cubicの有理translationに等しい。
P2九と任意のinhomogeneous particularを捨てた議論ではない。

01のB=a,A=bでは A(ζ) は全ρでnilpotent、k_eff=3ℓζ はKのunit。
translationと全m3/4/5作用素が可逆で、leading cubic ab(b+ζa)はpositive base factorsを持つ。
01のB=b,A=aでは A(ζ) は全ρでtriangular、固有値 σζ,−2σζ,σζ、
σ=2μ,κ=5μなのでk_eff/(σζ)=5/2。
Jordan化の仮定は不要で、triangular matrixの固有値だけで可逆性が決まる。

## 6. SF00一fixedの全higherも空、全R2L/m二SFの結論

00のfixed P0は X=H−2a∂c,Y=−H−2a∂c,Z=k20=0,κ=−5。
A(ζ)=−2ζH、k_eff=−5ζ、固有値係数s0=−2ζに対しk_eff/s0=5/2。
全 spectrum を入れた失敗比は次の通り。

| 条件 | 全 k_eff/s0 |
|---|---|
| translation | 1,−2 |
| m3、Sym² | −2,1,4 |
| m4、Sym¹ | −1/2,1 |
| m5、Sym⁰ | 0 |

5/2はいずれにも入らず、00及び前節semisimple-ratio01の全sourceで可逆。
ζをKの元として扱うため normal coordinate divisor は仮定として除外しない。
有理恒等式で得た因子分解は Gauss lemma で元 polynomial ring へ戻る。
00では F=−27uv(a+Ta)²(b+Tb)、01でも leading の h(次元二)×linear が残る。
両方ともK[a,b,c]に正のbase degreeの因子が二つ以上ある。
元Fがirreducibleでpositive base degreeなら C[u,v]上primitiveなので、Gauss lemmaが適用され矛盾する。

したがって00/01の残全higher源も空、02/03/11/12のglobal矛盾と合わせ
**全 R=a²b・m二 genericSFの被約既約 quintic source は空**。
normal double、R=a³、他normal order、全Poisson分類は本receiptからは結論しない。
全G/Eのaggregateは、それぞれ別のaccepted packetをrootが接続する義務である。
