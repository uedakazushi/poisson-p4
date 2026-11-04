# 42法線によるexact残余表示の独立監査

2026-11-04。凍結された `S61_SPARSE_GENERIC5_NORMAL42_EXACT_PRESENTATION_20261104_v63` の入力、typed operator、全chartの恒等式とschemeの降下を監査する。producer、Gram展開、126 Gaussian、新しい点・素数・GBは実行しない。旧 `S61_AFFINE_NORMAL_FINITE5_20261103_v63` と同じ抽象定理を、今回の実K5・130残差へsource-boundに実装した成果として扱う。

## 入力と型の点検

24501個の全ノードについて、参照順序、scalar/行列/多項式の型、行列shape、次数上界、normal原点でのorder下界を独立に再導出する。134式・129残差・G6の既受理source回路のunknown次数も読み、parameter scalarとして輸入する全根がz,nを含まないことを確かめる。

実Xの16380成分は、元129残差のcross129×126とG6残差のcross126の正確な参照である。b130は同じ130残差のz=0であり、各残差の既受理factor identityによりb(0)=0である。Gram、唯一のparameter線形系、P=−Lb、H=b+XPの入力と符号を確認する。新しいscalar係数を展開しない。

130本のH、2583本のcompatibility、210本のquarticを全て自由多項式環で照合する。Pの126成分とsource scalarは独立記号として扱うため、この検査は特殊化ではない。選択五式のconstant-in-z部分には全normal二次項を含み、linear-in-z部分には全cross成分を含む。2923根の列挙と次数≤2/3/4を確認する。

全42×1081=45402個の元次数≤2単項式について、n_j²による分母払いを独立な指数counterで検査する。さらに全3×42×42=5292個のchart recoveryの符号・反対称性を確認する。生成operatorにはnormal変数の逆元はなく、逆行列に相当する操作はparameter行列の一意線形求解だけである。

## Gramと四次元のpure-normal空間

K=QQ(H,L,s,q,t)は形式的実体で、Xは既受理のrank126を持つ130×126行列である。Cauchy–Binetの平方和に少なくとも一つ非零項があるためdet(XᵗX)≠0。従って(XᵗX)L=XᵗはK上で唯一の解を持ち、LX=Iとなる。左逆一般が唯一だという主張ではなく、このnormal equationによるcanonicalな選択である。

(I−XL)X=0なのでHの各成分はpure-normal行空間に属する。逆にvX=0ならvH=vbなので、H130は全pure-normal空間を張る。全135行の独立性、plane rank5、cross rank126からその次元は4である。四本の明示基底やGram係数は未計算である。

形式的実性を用いるのはK上のGram非零の証明だけである。generic complex係数体への拡大では非零性は保たれる。特定の複素parameterへの特殊化では、Xがfull rankでもGramが消える場合がある。従ってこのrecipeをその点へ移すには、既存adapter/selectorの全parameter分母とdet Gramが単元であることが別途必要である。全rank126特殊化への一様なGram recipeは主張しない。Σn_j²を反転する操作は一切ない。

## 全42chartのscheme同型

u_hk=z_h n_k、R130=Xu+bとする。AではLX=Iからu=P、H=0となり、compatibilityとquarticも零になる。逆にR/JのD(n_j)上でz_h=P_hj/n_jと定めると、compatibilityから全P_hk=z_h n_kとなる。従ってR130=X(u−P)+H=0であり、quarticと可逆なn_j²から選択五式も零となる。両方向の写像はnでは恒等、zも同じ式で回復する。全nilpotentを保つ商環の同型であり、radicalだけの一致ではない。

## 有限性と閉埋込み

既受理G6 reviewのstandard triangle証明から、選択五式の最高z二次部分はHF(1,3,1)、その一次倍は全三次を張る。最高z部分はparameter係数であり、normalを含む項はそのz次数を増やさない。よってAはR=K[n42]上で高々5生成である。自由性やK有限長は含まない。

actual known-five splitter A=E5×Bと、同じplane配置のA/(n)=E5を使用する。従ってBは有限R加群でB/mB=0、m=(n)である。局所化した中山の補題によりB_m=0となる。各D(n_j)ではE5が消え、前節からR_(n_j)→B_(n_j)は全射。R→Bのcokernelは有限で、mでは零、他のprimeではいずれかのn_jが単元なので零である。従ってR→Bは全射、Spec Bは原点を含まない閉部分schemeとなる。quasi-finiteだけからfiniteやclosed immersionを推論していない。

I_B=ker(R→B)はJと全D(n_j)で一致する。原点を避けることからI_B+m=Rであり、I_Bはm飽和している。有限個のchartの局所所属を分母払いすれば

`I_B = intersection_j (J:n_j^infinity) = J:m^infinity`

となる。逆包含もmB=Bから従う。従ってB≅R/(J:m^infinity)、B=0とあるdでm^d⊂Jは同値である。

この等式はexact scheme presentationである。飽和を実際に計算した、冪所属cofactorを得た、残余が空である、generic五点尽くしを証明した、という主張ではない。全2923式のscalar展開、四本の基底、minimal生成数も未計算である。

専用reviewコード・証拠は同名REVIEWに保存する。全input/evidence hashと実行時間・メモリをreport.jsonへ記録する。
