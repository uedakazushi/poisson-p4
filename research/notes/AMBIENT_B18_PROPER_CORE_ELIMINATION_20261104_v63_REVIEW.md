# 実 proper core の五逐次消去・七残余式の独立監査

PROPER_CORE_ELIMINATION の凍結 packet を、独立検証済み reviewer の sparse QQ 演算だけで
監査する。作者の worker、SymPy、CAS cancellation、GB、新 jet を使わない。

元の H12 のうち H0,H4 の (v3,v4) 係数行列を別構成し、det=-16(1-d²)(d²+5)、
原点値-80を確認する。五 pivot v4,v3,v5,v7,v8 の raw source 係数が実際に一次で、
未解決の後続 pivot を含まないことを検査する。各係数は非零 QQ 定数と D=1-d²、
P=d²+5 の冪だけであり、原点値は10,8,80,16,192である。

先行根を代入した五 source equations と保存された五根の恒等式を、分母を払う全係数
交差積で確認する。全五 Hsource の零性に加え、残り七つの元 Hrow について、保存
された rational numerator、D/P denominator powers、非零 QQ normalization と実 K を
独立に照合する。残余七多項式の原点 Jacobian は零である。

五つの最終根は全て base(a,d,v2,v6) のみの有理関数である。pivot_i=u_i+root_i への
ambient coordinate change は D,P が単元の完成環で可逆で、d を固定する。検証した
三角 unit 関係と残余恒等式から、全 H ideal は (u1,...,u5,K1,...,K7) に移る。
従って全 E4 は正確に

    QQ[[u1,...,u5,a,d,v2,v6]] / d*(u1,...,u5,K1,...,K7)

となる。divided core の商 QQ[[a,d,v2,v6]]/(K7) と全 E4 を混同せず、五つの d*u_i
関係も残す。actual D0、CM、depth、graph membership、次段の CAS 結果は本監査の
結論ではない。結果と入力 hash は report と independent-proof に保存する。
